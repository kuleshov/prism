#!/usr/bin/env python
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

import argparse
import re
from libprism.loader import parse_metadata

###############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--true-phase', required=True)
parser.add_argument('--blocks', required=True)
parser.add_argument('--threshold', type=float, default=0.8)

args = parser.parse_args()

###############################################################################

def most_common(lst):
    return max(set(lst), key=lst.count)

true_phase = dict()
with open(args.true_phase) as phase:
    for line in phase:
        fields = line.split()
        chrom = fields[0]
        pos = int(fields[1])
        ph = int(fields[2])

        true_phase[(chrom, pos)] = ph

global_pattern = ""
prev_pattern = ""
num_sw = 0
num_tr = 0
with open(args.blocks) as blocks:
    for line in blocks:
        pattern = ""
        fields = line.split()
        chrom = fields[0]
        md = parse_metadata(fields[3])

        for f in fields[4:]:
            fi = f.split(':')
            pos = int(fi[0])
            p0 = int(fi[1])

            if (chrom, pos) in true_phase:
                pattern += str(int(p0 == true_phase[(chrom, pos)]))

        if len(pattern) == 0:
            continue

        if prev_pattern == "":
            prev_pattern = pattern
            continue

        if float(md['TR']) < args.threshold:
            prev_pattern = pattern
            continue

        num_tr += 1
        if most_common(prev_pattern) != most_common(pattern):
            num_sw += 1

        prev_pattern = pattern

print 'Switch accuracy:', (num_tr - num_sw) / float(num_tr)