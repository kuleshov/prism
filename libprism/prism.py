#!/usr/bin/python
#
# Copyright (c) 2012-2013 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/downloads/sequencing/licenses/>.
#

import sys
import os
import argparse

from libprism.HMM import HMM
from libprism.panel import PhasedPanel, Subject
from libprism.positions import Positions
from libprism import loader

from math import log, exp

def main():
    ###################################################################################################
    # ARGUMENTS

    parser = argparse.ArgumentParser(description="Prism statistical phaser")

    parser.add_argument('--blocks', help="locally phased .bed",required=True)
    parser.add_argument('--positions', help='reference file',required=True)
    parser.add_argument('--start', help='starting position in b37 coordinates', type=int, required=True)
    parser.add_argument('--end', help='ending position in b37 coordinates', type=int, required=True)
    parser.add_argument('--phase', help='output file containg the phase')
    parser.add_argument('--path', help='output file containing the Viterbi path through the reference panel')

    parser.add_argument('--K', help='panel size', type=int, default=100)
    parser.add_argument('--N', help='effective population size', type=float, default=15000)
    parser.add_argument('--visualize-panel', action='store_true')

    args = parser.parse_args()

    ###################################################################################################
    ## LOAD POSITIONS, PANEL, AND SUBJECT

    print '[prism] Starting new phasing run.'
    positions = Positions(args.positions, args.start, args.end)
    print '[prism] Loaded set of positions'

    haplotypes, genetic_map = loader.load_haplotypes_and_genetic_map_from_panel_gz(args.positions, positions)
    print '[prism] Loaded locally-phased haplotypes.'

    subject_blocks, subject_hap0, subject_hap1 = loader.load_blocks_and_haplotypes(args.blocks, args.positions, positions)
    subject = Subject(positions, subject_hap0, subject_hap1, subject_blocks)

    ###################################################################################################
    ## LOAD PHASED REFERENCE PANEL

    block_scores = loader.compute_block_scores(haplotypes, subject, positions)
    min_haplotype_distances = loader.compute_haplotype_distances(block_scores)
    haplotype_counts, min_indices, min_distances = loader.choose_smallest_haplotypes(haplotypes, min_haplotype_distances,
                                                                            args.K)

    k_haplotypes = [haplotypes[i] for i in min_indices]
    hap_counts = [haplotype_counts[h] for h in k_haplotypes]

    panel = PhasedPanel(positions, k_haplotypes, hap_counts, args.N, haplotype_distances=min_distances)
    print '[prism] Loaded reference panel.'

    if args.visualize_panel:
        print '[prism] Panel visualization:'
        panel.visualize(subject)

    ###################################################################################################
    ## RUN VITERBI, FORWARDS, BACKWARDS

    hmm = HMM(panel, subject, genetic_map, positions)
    print '[prism] Loaded HMM model.'

    # Run Viterbi, forwards, and backwards
    print '[prism] Starting the Viterbi algorithm...'
    hmm.run_viterbi_factorial()
    print '[prism] Viterbi completed.'
    print '[prism] Starting the forwards algorithm...'
    hmm.run_forwards_factorial()
    print '[prism] Forwards completed.'
    print '[prism] Starting the backwards algorithm...'
    hmm.run_backwards_factorial()
    print '[prism] Backwards completed.'

    if abs(hmm.fwd_log_data_prob - hmm.bwd_log_data_prob) > 1e-4:
        exit("An error has occured: the probabilities obtained from the HMM are inconsistent.")
    print '[prism] The forwards and backwards data probabilities are matching.'

    ###################################################################################################
    ## SAVE THE PHASE OF THE BLOCKS

    out_file = open(args.phase,"w")
    out_file.write('#ID\tPHASE\tSCORE\n')
    tr_scores = dict()

    for pos, z_index in zip(positions, hmm.haplotype):
        if subject.is_block_start(pos):
            if not positions.is_start_position(pos):
              qscore_tr = exp(hmm.get_transition_lik_z_only(pos, z_index, z_index_prev))
            else:
              qscore_tr = 1
            block = subject.block_at_start_position(pos)
            out_file.write("%s\t%d\t%f\n" % (block.name, z_index, qscore_tr))

            tr_scores[pos] = qscore_tr  # to avoid re-computing this at the next step (costly)

        z_index_prev = z_index

    out_file.close()

    print '[prism] Phase of local blocks saved.'

    ###################################################################################################
    ## SAVE VITERBI PATH, DISPLAY VISUALIZATION

    ref_file = open(args.path, "w")
    ref_file.write('#POS\tREF1\tREF2\tREFALLELE1\tREFALLELE2\n')
    short_to_full_ind = min_indices     # for translating between indices in the panel of size K and the full panel

    for pos, (y0, y1), z_index in zip(positions, hmm.path, hmm.haplotype):
        qscore_tr = tr_scores.get(pos,1.0)

        panel_at_pos = panel.get_matrix_at_position(pos)
        ref_file.write("%d\t%d\t%d\t%d\t%d\n" % (pos, short_to_full_ind[y0], short_to_full_ind[y1],
                                                 panel_at_pos[y0], panel_at_pos[y1]))

    ref_file.close()

    print '[prism] Path through reference panel saved.'
    print '[prism] Prism finished successfully.'

    # for pos, (y0, y1), z_index in zip(positions, hmm.path, hmm.haplotype):
    #     block_start = " "
    #     qscore_tr = tr_scores.get(pos,1.0)

    #     if subject.is_block_start(pos): block_start = "X"

    #     panel_at_pos = panel.get_matrix_at_position(pos)
    #     subject_at_pos = subject.haplotype_at_position(pos)

    #     print "%d\t%d\t%d %d\t%d %d\t%d %d\t%.4f\t%d %s" % \
    #           (pos, positions.get_position_index(pos)+1,
    #            subject_at_pos[0], subject_at_pos[1],
    #            panel_at_pos[y0], panel_at_pos[y1],
    #            y0, y1,
    #            qscore_tr,
    #            z_index, block_start)

if __name__ == '__main__':
    main()