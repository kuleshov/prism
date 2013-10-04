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

from libprism.loader import parse_metadata

##############################################################################

parser = argparse.ArgumentParser()

parser.add_argument('--threshold', type=float)
parser.add_argument('--phased-local-blocks')
parser.add_argument('--global-blocks')

args = parser.parse_args()

##############################################################################

def write_block(out, global_block, id_):
	if not global_block:
		return

	global_block.sort()
	
	first_pos = global_block[0][0]
	last_pos = global_block[-1][0]

	out.write('%s\t%d\t%d\tID=%d' % ('chr22', first_pos, last_pos, id_))
	for pos, p0, p1 in global_block:
		out.write('\t%d:%d:%d' % (pos, p0, p1))

	out.write('\n')

current_block = list()
num_global_blocks = 0
out = open(args.global_blocks, 'w')
with open(args.phased_local_blocks) as blocks:
	for line in blocks:
		fields = line.split()

		chrom = fields[0]
		metadata = parse_metadata(fields[3])
		id_ = int(metadata['ID'])
		tr_score = float(metadata['TR'])

		if tr_score < args.threshold and current_block:
			write_block(out, current_block, num_global_blocks)
			current_block = list()
			num_global_blocks += 1

		for f in fields[4:]:
			fi = f.split(':')
			
			pos = int(fi[0])
			p0 = int(fi[1])
			p1 = int(fi[2])

			current_block.append((pos, p0, p1))

if current_block:
	write_block(out, current_block, num_global_blocks)

out.close()