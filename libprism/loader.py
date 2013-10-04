# Copyright (c) 2012-2013 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/downloads/sequencing/licenses/>.
#

import numpy as np

from libprism.panel import Haplotype, Member, GeneticMap, Block

###################################################################################################
## FUNCTIONS THAT ARE USED TO LOAD A FULL PANEL

def parse_metadata(metadata):
    """Returns a dictionary of metadata.

    Metadata needs to be in the format k1=v1:k2=v2:...:kn=vn.
    """
    tags = dict()
    fields = metadata.split(':')
    for field in fields:
        key, value = field.split('=')
        tags[key] = value

    return tags

def count_haplotypes(pos_file_path):
    pos_file = open(pos_file_path)
    line = pos_file.readline()
    fields = line.strip().split(',')
    ref_panel = fields[-1].split()
    pos_file.close()

    return len(ref_panel)

def load_members_from_haplotypes(haplotypes, positions):
    H = len(haplotypes)
    assert H % 2 == 0
    N = H / 2

    members = [Member(positions, haplotypes[2*x], haplotypes[2*x+1]) for x in xrange(N)]
    return members

def load_haplotypes_and_genetic_map_from_panel_gz(pos_file_path, positions):
    H = count_haplotypes(pos_file_path)

    haplotypes = [Haplotype(positions) for x in xrange(H)]
    genetic_map = GeneticMap(positions)

    pos_file = open(pos_file_path)
    for line in pos_file:
        fields = line.split(',')
        pos = int(fields[0])
        if pos not in positions: continue

        for h, x in zip(haplotypes, fields[-1].split()):
            h[pos] = x
        genetic_map.set_distance(pos, float(fields[3]))

    return haplotypes, genetic_map

def load_blocks_and_haplotypes(block_file_path, pos_file_path, positions):
    blocks = list()
    hap0 = Haplotype(positions)
    hap1 = Haplotype(positions)

    # load blocks and subject heterozygous positions from block file
    block_file = open(block_file_path)
    for line in block_file:
        fields = line.strip().split()

        block_positions = list()
        metadata = parse_metadata(fields[3])

        # create a separate positions that restricts range?
        h0, h1 = Haplotype(positions), Haplotype(positions)
        for f in fields[4:]:
            fi = f.split(':')
            pos = int(fi[0])
            if not positions.first_position <= pos <= positions.last_position:
                continue
            block_positions.append(pos)

            p0, p1 = int(fi[1]), int(fi[2])
            
            h0[pos], h1[pos] = p0, p1
            hap0[pos], hap1[pos] = p0, p1

        if len(block_positions) > 0:
            block_positions.sort()
            first_position, last_position = block_positions[0], block_positions[-1]
            block = Block(positions, name=metadata['ID'], start=first_position,
                          end=last_position, h0=h0, h1=h1, num_snps=len(block_positions))
            blocks.append(block)
    block_file.close()

    pos_file = open(pos_file_path)
    for line in pos_file:
        fields = line.split(',')
        pos = int(fields[0])
        if not positions.first_position <= pos <= positions.last_position:
            continue

        if fields[2] == '(1/1)':
            hap0[pos] = 1
            hap1[pos] = 1
    pos_file.close()

    return blocks, hap0, hap1

###################################################################################################
## FUNCTIONS THAT ARE USED TO COMPUTE THE CLOSEST HAPLOTYPES TO THE SUBJECT

def compute_block_scores(haplotypes, subject, positions):
    # print len(positions)
    # print len(haplotypes)
    # print len(subject.blocks)
    FP = np.zeros((len(positions), len(haplotypes)), dtype=int)
    for k, h in enumerate(haplotypes):
        FP[:,k] = h.get_vector().squeeze()
    FP = FP.transpose()
    S = subject.get_haplotype_matrix()

    sorted_starts = [positions.get_position_index(block.start) for block in subject.blocks]
    sorted_starts.sort()
    segments = [(sorted_starts[i], sorted_starts[i+1]-1) for i in xrange(len(sorted_starts)-1)]
    segments += [(sorted_starts[-1], positions.get_position_index(positions.last_position))]

    block_scores = [[0,0] for s in segments]

    for b, (start, end) in enumerate(segments):
        for h in (0,1):
            # block_scores[b][h] is an Hx1 matrix
            block_scores[b][h] = np.sum(np.bitwise_xor(FP[:,start:end+1],S[start:end+1,h]),axis=1)

    return block_scores

def compute_haplotype_distances(block_distances):
    num_haplotypes = len(block_distances[0][0])
    D = np.zeros((num_haplotypes,2))
    for b, block_distance in enumerate(block_distances):
        D_prev = D.copy()

        for h in (0,1):
            if b == 0:
                D[:,h] = block_distance[h]
            else:
                D[:,h] = np.min(D_prev,axis=1)
                D[:,h] += block_distance[h]

    return np.min(D,axis=1)

def choose_smallest_haplotypes(haplotypes, distances, K):
    sorted_indices = [x[0] for x in sorted(enumerate(distances), key=lambda x:x[1], reverse=True)]

    min_distances = list()
    min_indices = list()
    haplotype_counts = dict()

    while len(min_indices) < K and len(sorted_indices) > 0:
        i = sorted_indices.pop()
        if haplotypes[i] not in haplotype_counts:
            haplotype_counts[haplotypes[i]] = 1
            min_indices.append(i)
            min_distances.append(distances[i])
        else:
            haplotype_counts[haplotypes[i]] += 1

    return haplotype_counts, min_indices, min_distances