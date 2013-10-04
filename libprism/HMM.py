# Copyright (c) 2012-2013 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/downloads/sequencing/licenses/>.
#

import os
import sys
import numpy as np

np.seterr(divide='ignore')

from math import log
from math import exp
from numpy import logaddexp
from scipy.weave import inline, converters

from libprism.math import logsumexp, H

class HMM(object):
    def __init__(self, panel, subject, genetic_map, positions):
        self.panel = panel
        self.subject = subject
        self.genetic_map = genetic_map
        self.positions = positions

        self.back_pointers = None
        self.haplotype = None # optimal viterbi path
        self.fwd_values = None
        self.bwd_values = None

        self.dir = os.path.dirname(os.path.realpath(__file__))

        self.collect_common_parameters()

    def collect_common_parameters(self):
        panel = self.panel
        common_arguments = dict()

        common_arguments['K'] = panel.get_haplotype_num()
        common_arguments['S'] = panel.get_haplotype_num()**2
        common_arguments['N_e'] = panel.get_effective_population_size()

        # variable lam (for lambda) is used to compute emission probabilities:
        K = panel.get_haplotype_num()
        theta = 1.0/H(K-1)
        common_arguments['lam'] = theta/(2*(theta+K)) * 10 # NOTE: times ten optional here

        common_arguments['hap_counts'] = panel.get_haplotype_counts()
        common_arguments['T'] = sum(panel.get_haplotype_counts())

        self.common_arguments = common_arguments

    def run_inner_loop(self, pos, code, algorithm_arguments):
        # panel = self.panel
        # subject = panel.get_subject()

        # merge common, position-specific, and algorithm-specific requirements in one dict
        arguments = dict()
        arguments.update(self.common_arguments)
        arguments.update(algorithm_arguments)

        # run inner loop
        # for key, value in arguments.items():
        #     print key, value

        return inline(code,
               arg_names = arguments.keys(),
               local_dict = arguments,
               headers = ["<math.h>", "<stdio.h>",
                          "\"" + self.dir + "/loops/indices.h\"",
                          "\"" + self.dir + "/loops/logaddexp.h\""],
               type_converters = converters.blitz,
               compiler = "gcc")
        #extra_compile_args=["-pg", "-O3", "-finline-functions",
        #"-ffast-math"]

    def run_viterbi_factorial(self):
        positions = self.positions
        subject = self.subject
        is_start_position = positions.is_start_position
        K = self.panel.get_haplotype_num()

        M = np.zeros((K**2,2))
        M_fact = np.zeros((K,K))
        M_fact_pointer = np.zeros((K,K))
        back_pointers = dict()

        # load inner loop C code:
        code_file = open(self.dir + "/loops/viterbi_factorial.c")
        viterbi_code = '\n'.join(code_file.readlines())
        code_file.close()

        # start Viterbi:
        for j, pos in enumerate(positions):
            if j % 200 == 0:
                print '%d/%d' % (positions.get_position_index(pos), len(positions))
            sys.stdout.flush()

            back_pointers[pos] = np.zeros((K**2,2))
            bpointers = back_pointers[pos]
            M_prev = M.copy()

            log_prob_equal, log_prob_nequal = self.get_transition_liks(pos)
            prob_equal, prob_nequal = exp(log_prob_equal), exp(log_prob_nequal)

            block_start_position = int(subject.is_block_start(pos))
            start_position = int(is_start_position(pos))
            panel = self.panel.get_matrix_at_position(pos)
            sample = subject.haplotype_at_position(pos)

            self.run_inner_loop(pos, viterbi_code, locals())

        # recover haplotype by following backpointers:

        # functions to convert between 2d and 1d indices
        sub2ind = lambda x: (int(x) / 2, int(x) % 2)
        y2hap = lambda x: (x / K, x % K)

        y_max_index, z_max_index = sub2ind(M.argmax())

        haplotype = list()
        path = list()
        haplotype.append(z_max_index)
        path.append(y2hap(y_max_index))

        pos = positions.last_position
        while not positions.is_start_position(pos):
            idx = back_pointers[pos][y_max_index, z_max_index]
            y_max_index, z_max_index = sub2ind(idx)
            haplotype.append(z_max_index)
            path.append(y2hap(y_max_index))
            pos = positions.get_previous_position(pos)
            
        haplotype.reverse()
        path.reverse()

        self.haplotype = haplotype
        self.path = path

    def run_forwards_factorial(self):
        positions = self.positions
        subject = self.subject
        is_start_position = positions.is_start_position
        K = self.panel.get_haplotype_num()

        fwd_values = dict()
        fwd_prev = np.zeros((K**2,2))

        # load c code:
        code_file = open(self.dir + "/loops/forwards_factorial.c")
        forwards_code = '\n'.join(code_file.readlines())
        code_file.close()

        # actual forwards algorithm starts here:
        for j, pos in enumerate(positions):
            if j % 200 == 0:
                print '%d/%d' % (positions.get_position_index(pos), len(positions))
            sys.stdout.flush()

            if not is_start_position(pos):
                prev_pos = positions.get_previous_position(pos)
            else:
                # the C code doesn't use any variables, we just need to fill them with something
                prev_pos = pos

            fwd = np.zeros((K**2,2))
            fwd_values[pos] = fwd

            log_norm_const = logsumexp(fwd_values[prev_pos])
            fwd_prev = np.exp(fwd_values[prev_pos] - log_norm_const)

            log_prob_equal, log_prob_nequal = self.get_transition_liks(pos)
            prob_equal, prob_nequal = exp(log_prob_equal), exp(log_prob_nequal)

            block_start_position = int(subject.is_block_start(pos))
            start_position = int(positions.is_start_position(pos))
            panel = self.panel.get_matrix_at_position(pos)
            sample = subject.haplotype_at_position(pos)

            self.run_inner_loop(pos, forwards_code, locals())

            if not is_start_position(pos):
                fwd_values[pos] = np.log(fwd_values[pos]) + log_norm_const

        last_pos = positions.last_position
        log_data_prob = logsumexp(fwd_values[last_pos])

        self.fwd_values = fwd_values
        self.fwd_log_data_prob = log_data_prob

    def run_backwards_factorial(self):
        positions = self.positions
        subject = self.subject
        is_end_position = positions.is_end_position
        K = self.panel.get_haplotype_num()

        bwd_values = dict()
        for pos in positions:
            bwd_values[pos] = np.zeros((K**2,2))

        # load c code:
        code_file = open(self.dir + "/loops/backwards_factorial.c")
        backwards_code = '\n'.join(code_file.readlines())
        code_file.close()

        # actual backwards algorithm starts here:
        for j, pos in enumerate(reversed(positions)):
            if j % 200 == 0:
                print '%d/%d' % (positions.get_position_index(pos), len(positions))
            sys.stdout.flush()

            if not is_end_position(pos):
                next_pos = positions.get_next_position(pos)
            else:
                # the C code doesn't use any variables, we just need to fill them with something
                next_pos = pos

            block_start_position = int(subject.is_block_start(next_pos))
            start_position = int(positions.is_start_position(next_pos))
            panel = self.panel.get_matrix_at_position(next_pos) #TODO: rename panel to matrix here and in C code
            sample = subject.haplotype_at_position(next_pos)

            log_prob_equal, log_prob_nequal = self.get_transition_liks(next_pos)
            prob_equal, prob_nequal = exp(log_prob_equal), exp(log_prob_nequal)

            log_norm_const = logsumexp(bwd_values[next_pos])
            bwd_next = np.exp(bwd_values[next_pos] - log_norm_const)

            bwd = bwd_values[pos]
            end_position = int(is_end_position(pos))

            self.run_inner_loop(pos, backwards_code, locals())

            if not is_end_position(pos):
                bwd_values[pos] = np.log(bwd_values[pos]) + log_norm_const

        # log_data_prob:
        log_data_prob = -np.inf
        y_states = [(i,j) for i in xrange(K) for j in xrange(K)]
        z_states = (0,1), (1,0)

        start_pos = positions.first_position
        sample = subject.haplotype_at_position(start_pos)
        panel = self.panel.get_matrix_at_position(start_pos)
        lam = self.common_arguments['lam']

        for y_index, y in enumerate(y_states):
            for z_index, z in enumerate(z_states):

                # compute emission probability:
                h = panel[y[0]], panel[y[1]]
                log_emission_prob = 0

                for x in (0,1):
                    if sample[x] == h[z[x]]:
                        log_emission_prob += log(1-lam)
                    else: 
                        log_emission_prob += log(lam)

                log_data_prob = logaddexp(log_data_prob,
                                          bwd_values[start_pos][y_index, z_index]
                                          + log_emission_prob
                                          + log(0.5) + log(1.0/(K**2))
                                          # + log(hap_counts[y[0]]/T)
                                          # + log(hap_counts[y[1]]/T)
                                          )

                if np.isnan(log_data_prob): exit("ERROR: log_data_prob=NaN")

        self.bwd_values = bwd_values
        self.bwd_log_data_prob = log_data_prob
        self.log_data_prob = log_data_prob

    def get_posterior_code(self,pos,y_index,z_index):
        return exp(
            self.fwd_values[pos][z_index][y_index] + self.bwd_values[pos][z_index][y_index] \
            - self.log_data_prob)

    def get_posterior_lik_z_only(self,pos,z_index):
        return logsumexp(self.fwd_values[pos][:,z_index] + self.bwd_values[pos][:,z_index] - self.log_data_prob)

    def get_posterior_joint_lik_z_only(self,pos,z,z_prev):
        subject = self.subject
        K = self.panel.get_haplotype_num()

        # load c code:
        code_file = open(self.dir + "/loops/joint_probability.c")
        probability_code = '\n'.join(code_file.readlines())
        code_file.close()

        prev_pos = self.positions.get_previous_position(pos)

        fwd_const = logsumexp(self.fwd_values[prev_pos][:,z_prev])
        bwd_const = logsumexp(self.bwd_values[pos][:,z])

        fwd_prev = np.exp(self.fwd_values[prev_pos] - fwd_const)
        bwd = np.exp(self.bwd_values[pos] - bwd_const)

        panel = self.panel.get_matrix_at_position(pos)
        sample = subject.haplotype_at_position(pos)

        log_prob_equal, log_prob_nequal = self.get_transition_liks(pos)
        prob_equal, prob_nequal = exp(log_prob_equal), exp(log_prob_nequal)

        joint_prob = self.run_inner_loop(pos, probability_code, locals())

        log_joint_prob = log(joint_prob) + fwd_const + bwd_const
        return log_joint_prob - self.log_data_prob

    def get_transition_lik_z_only(self,pos,z,z_prev):
        prev_pos = self.positions.get_previous_position(pos)
        joint_lik = self.get_posterior_joint_lik_z_only(pos,z,z_prev)
        z_lik = self.get_posterior_lik_z_only(prev_pos,z_prev)

        return joint_lik - z_lik

    def get_transition_liks(self, pos):
        """ Returns transition likelihoods from position pos to position prev(pos),
        or 0.0, 0.0 if pos is the first position. """
        if self.positions.is_start_position(pos):
            return 0.0, 0.0

        panel = self.panel
        subject = self.subject
        genetic_map = self.genetic_map
        prev_pos = self.positions.get_previous_position(pos)

        N_e = panel.get_effective_population_size()
        N = panel.get_haplotype_num()

        if not subject.is_block_start(pos): N_e = 10
        # try: N_e = 0; N_e = old_N_e

        rho = 4*N_e*(genetic_map[pos]-genetic_map[prev_pos])/100.0
        if rho > 0.0:
            NR = exp(-rho/N)
            log_prob_equal = log(NR + (1-NR)/N)
            log_prob_nequal = log((1-NR)/N)
        else:
            log_prob_equal = 0.0
            log_prob_nequal = -np.inf

        return log_prob_equal, log_prob_nequal
