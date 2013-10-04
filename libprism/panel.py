# Copyright (c) 2012-2013 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/downloads/sequencing/licenses/>.
#

""" Classes that express the panel that is to be phased.
"""

import numpy as np

###################################################################################################

class Panel:
    """ Panel to be phased. This object will be used later for phasing with
    different types of panels (1000G + HapMap), including ones that are unphased.

    CURRENTLY, IT IS NOT USED. """

    def __init__(self, positions, members, N_e):
        self.positions = positions

        self.members = members

        haplotypes = list()
        for member in self.members:
            self.haplotypes.extend(member.get_haplotypes())
        self.haplotypes = haplotypes

        self.N_e = N_e  # effective population size

    def get_full_matrix(self):
        matrix = np.zeros(len(self.positions), len(self.haplotypes))
        for k, h in enumerate(self.haplotypes):
            matrix[:,k] = h.get_vector()

        return matrix

###################################################################################################

class PhasedPanel:
    """ Panel of phased reference haplotypes. """

    def __init__(self, positions, haplotypes, haplotype_counts, N_e, haplotype_distances=None):
        self.positions = positions
        self.haplotypes = haplotypes

        matrix = np.zeros((len(self.positions), len(self.haplotypes)), dtype=int)
        for k, h in enumerate(self.haplotypes):
            matrix[:,k] = h.get_vector().squeeze()
        self.haplotype_matrix = matrix

        self.haplotype_counts = haplotype_counts
        self.N_e = N_e  # effective population size

        if haplotype_distances is None:
            self.haplotype_distances = [0 for h in haplotypes]
        else:
            self.haplotype_distances = haplotype_distances

    def get_haplotype_num(self):
        return len(self.haplotypes)

    def get_effective_population_size(self):
        return self.N_e

    def get_haplotype_counts(self):
        return self.haplotype_counts

    def get_matrix_at_position(self, pos):
        verify_position(pos, self.positions)
        j = self.positions.get_position_index(pos)
        return self.haplotype_matrix[j,:]

    def get_full_matrix(self):
        return self.haplotype_matrix

    def visualize(self, subject):
        print "%s" % ' '.join(["(%f,%d)" % (d,c) for d,c in zip(self.haplotype_distances, self.haplotype_counts)])
        for pos in self.positions:
            if subject.is_block_start(pos):
                block_start = "X"
            else:
                block_start = " "
            h = subject.haplotype_at_position(pos)
            M = self.get_matrix_at_position(pos)
            print "%d\t%d\t%d %d\t%s %s" % (pos,self.positions.get_position_index(pos),
                                            h[0], h[1], ''.join([str(x) for x in M]),block_start)

###################################################################################################

class GeneticMap:
    def __init__(self, positions):
        self.positions = positions
        self.genetic_map = dict()

    def __getitem__(self, pos):
        verify_position(pos, self.positions)
        return self.genetic_map[pos]

    def get_distance(self, pos):
        verify_position(pos, self.positions)
        return self.genetic_map[pos]

    def set_distance(self, pos, d):
        verify_position(pos, self.positions)
        self.genetic_map[pos] = d

###################################################################################################

class Member(object):
    """ Member of the reference panel. """

    def __init__(self, positions, h0, h1):
        self.positions = positions
        self.h0 = h0
        self.h1 = h1

    def get_haplotypes(self):
        return h0, h1

    def get_haplotype_vectors(self):
        return self.h0.get_vector(), self.h1.get_vector()

    def get_haplotype_matrix(self):
        matrix = np.zeros((len(self.h0),2), dtype=int)
        matrix[:,0] = self.h0.get_vector().squeeze()
        matrix[:,1] = self.h1.get_vector().squeeze()
        return matrix

    def haplotype_at_position(self, pos):
        verify_position(pos, self.positions)
        return np.array((self.h0[pos], self.h1[pos]))

    def genotype_at_position(self, pos):
        verify_position(pos, self.positions)

###################################################################################################

class Subject(Member):
    def __init__(self, positions, h0, h1, blocks):
        super(Subject, self).__init__(positions, h0, h1)
        self._blocks = sorted(blocks, key=lambda b: b.start)

        self._block_starts = set()
        self._blocks_at_start = dict()
        for block in self._blocks:
            self._block_starts.add(block.start)
            self._blocks_at_start[block.start] = block

    @property
    def blocks(self):
        return self._blocks

    def is_block_start(self,pos):
        verify_position(pos, self.positions)
        return (pos in self._block_starts)

    def block_at_start_position(self, pos):
        return self._blocks_at_start[pos]

###################################################################################################

class Haplotype:
    def __init__(self, positions):
        self.positions = positions
        self.seq = np.zeros((len(positions),1), dtype=int)

    def __setitem__(self, pos, x):
        verify_position(pos, self.positions)
        idx = self.positions.get_position_index(pos)
        self.seq[idx] = x

    def __getitem__(self, pos):
        verify_position(pos, self.positions)
        idx = self.positions.get_position_index(pos)
        return self.seq[idx]

    def __len__(self):
        return len(self.seq)

    def set_vector(self, v):
        assert np.shape(v) == np.shape(self.seq)
        self.seq = v

    def get_vector(self):
        return self.seq

    def __hash__(self):
        l = tuple([int(x) for x in self.seq])
        return l.__hash__()

    def __eq__(self, other):
        return all(other.seq == self.seq)

###################################################################################################

class Block:
    def __init__(self, positions, name, start, end, h0, h1, num_snps):
        self.__name = name
        self.__start = start
        self.__end = end
        self.h0 = h0 #TODO: do we need haplotypes here?
        self.h1 = h1
        self.positions = positions
        self.phase = -1
        self.__num_snps = num_snps

    @property
    def start(self):
        return self.__start

    @property
    def end(self):
        return self.__end

    @property
    def name(self):
        return self.__name

    @property
    def num_snps(self):
        return self.__num_snps

## HELPER FUNCTIONS ###############################################################################

def verify_position(pos, positions):
    if pos not in positions:
        exit("Invalid position:" + (pos))
