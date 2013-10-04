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

args = parser.parse_args()

###############################################################################

def most_common(lst):
    return max(set(lst), key=lst.count)

def correct_pattern(pattern):
    ones_at_ends_pattern=re.compile('^10|01$')
    zeros_at_ends_pattern=re.compile('^01|10$')
    pattern = ones_at_ends_pattern.sub('00', pattern)
    pattern = zeros_at_ends_pattern.sub('11', pattern)

    fives_pattern0 = re.compile('01010')
    fives_pattern1 = re.compile('10101')
    pattern = fives_pattern0.sub('00000', pattern)
    pattern = fives_pattern1.sub('11111', pattern)

    single_one_pattern = re.compile('010')
    single_zero_pattern = re.compile('101')
    pattern = single_one_pattern.sub('000', pattern)
    pattern = single_zero_pattern.sub('111', pattern)

    return pattern

###############################################################################

true_phase = dict()
with open(args.true_phase) as phase:
	for line in phase:
		fields = line.split()
		chrom = fields[0]
		pos = int(fields[1])
		ph = int(fields[2])

		true_phase[(chrom, pos)] = ph

num_switches = 0
num_long_switches = 0
num_transitions = 0
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

		# print md['ID'], '\t', pattern

		if len(pattern) > 1:
			corrected_pattern = correct_pattern(pattern)
			num_block_switches = sum([1 for i in xrange(1,len(pattern)) if pattern[i] != pattern[i-1]])
			num_long_block_switches = sum([1 for i in xrange(1,len(corrected_pattern)) 
										   if corrected_pattern[i] != corrected_pattern[i-1]])
			num_block_transitions = len(pattern) - 1

			num_switches += num_block_switches
			num_long_switches += num_long_block_switches
			num_transitions += num_block_transitions

print 'Switch accuracy:', (num_transitions - num_switches) / float(num_transitions)
print 'Long switch accuracy:', (num_transitions - num_long_switches) / float(num_transitions)