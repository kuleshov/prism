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

parser.add_argument('--phase')
parser.add_argument('--unphased-local-blocks')
parser.add_argument('--phased-local-blocks')

args = parser.parse_args()

##############################################################################

block_phase = dict()
block_score = dict()
with open(args.phase) as phase:
	for line in phase:
		if line.startswith('#'): continue

		fields = line.strip().split()

		id_ = int(fields[0])
		ph = int(fields[1])
		score = float(fields[2])

		block_phase[id_] = ph
		block_score[id_] = score

def metadata_to_str(metadata):
	return ':'.join(['%s=%s' % (str(k), str(v)) for k, v in metadata.iteritems()])

phased_blocks = open(args.phased_local_blocks, 'w')

with open(args.unphased_local_blocks) as blocks:
	for line in blocks:
		fields = line.split()

		metadata = parse_metadata(fields[3])
		id_ = int(metadata['ID'])
		if id_ not in block_phase:
			continue

		metadata['TR'] = '%4f' % block_score[id_]
		ph = block_phase[id_]

		phased_blocks.write('\t'.join(fields[0:3]))
		phased_blocks.write('\t' + metadata_to_str(metadata))

		for f in fields[4:]:
			fi = f.split(':')

			pos = int(fi[0])
			p0 = (int(fi[1]) + ph) % 2
			p1 = (int(fi[2]) + ph) % 2

			phased_blocks.write('\t%d:%d:%d' % (pos, p0, p1))

		phased_blocks.write('\n')


			
