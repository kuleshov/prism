#!/usr/bin/python

import os
import numpy as np
import argparse
import bisect

from libprism.loader import parse_metadata

def main():

	######################################################################

	parser = argparse.ArgumentParser()

	parser.add_argument('--unphased-blocks')
	parser.add_argument('--path')
	parser.add_argument('--vcf')
	parser.add_argument('--chr')
	parser.add_argument('--reference-panel')
	parser.add_argument('--reference-legend')
	parser.add_argument('--reference-sample')
	parser.add_argument('--phased-blocks')
	parser.add_argument('--population', default='ALL')
	parser.add_argument('--leave_out', help='IDs of samples to leave out')

	args = parser.parse_args()

	######################################################################
	## LOAD TRUE VCF BASES

	def load_letters_by_allele(vcf_filename, chrom):
		letter_by_allele = dict()
		one_two_positions = set()

		cmd = "cat %s | grep -P '^%s\\t' | cut -f 2,4,5,10 | cut -d ':' -f 1 |"\
			  "grep -v 1/1 | grep -v 1\|1 | grep -v 0/0 | grep -v 0\|0 | grep -v 1\|2 |"\
			  "grep -v -P '1/\.'" % (vcf_filename, chrom)
		vcf_file = os.popen(cmd)

		for j, line in enumerate(vcf_file):
			fields = line.strip().split()

			pos = int(fields[0])
			if fields[3] == "0/1" or fields[3] == "1/0"\
			   or fields[3] == "0|1" or fields[3] == "1|0":
				letter_by_allele[pos] = (fields[1], fields[2])
			elif fields[3] == "1/2":
				letters = fields[2].split(',')
				letter_by_allele[pos] = (letters[0], letters[1])
				one_two_positions.add(pos)
			else:
				print "ERROR: Could not process this line in the VCF:"
				print line
				exit(1)

		return letter_by_allele, one_two_positions

	letter_by_allele, one_two_positions = load_letters_by_allele(args.vcf, args.chr)

	######################################################################

	# collect positions that need to be phased
	unphased_positions = set()

	# and collect blocks
	blocks = list()
	block_lines = list()

	print "Loading unphased blocks..."
	block_file = open(args.unphased_blocks)
	for line in block_file:
		block = list()
		fields = line.split()
		for f in fields[4:]:
			fi = f.split(':')
			pos, p0, p1 = int(fi[0]), int(fi[1]), int(fi[2])
			unphased_positions.add(pos)
			block.append((pos,p0,p1))
		blocks.append(block)
		block_lines.append(line)
	block_file.close()
	print "Loaded %d blocks covering %d unphased positions." % (len(blocks), len(unphased_positions))

	######################################################################
	# collect reference indices at phased positions

	print "Loading reference indices at phased positions..."
	ref_index_file = open(args.path)
	reference_indices_for_phased = dict()
	phased_positions = list()
	for line in ref_index_file:
		if line.startswith('---'): continue
		if line.startswith('#'): continue
		fields = line.split()
		pos, i1, i2, p0, p1 = (int(x) for x in fields)
		phased_positions.append(pos)
		reference_indices_for_phased[pos] = i1, i2
	ref_index_file.close()
	print "Loaded indices at %d positions." % len(reference_indices_for_phased.keys())

	######################################################################
	# assign each unphased position a reference index:

	print "Assigning reference indices to unphased positions..."
	num_pos_lost_to_recombination = 0
	reference_indices_for_unphased = dict()
	for missing_pos in sorted(list(unphased_positions)):
		j = bisect.bisect_left(phased_positions, missing_pos)

		if j >= len(phased_positions)-1 or j==0:
			num_pos_lost_to_recombination += 1
			continue

		pos_prev = phased_positions[j]
		pos_next = phased_positions[j+1]

		# print reference_indices_for_phased[pos_prev], reference_indices_for_phased[pos_next]
		if reference_indices_for_phased[pos_prev] == reference_indices_for_phased[pos_next]:
			reference_indices_for_unphased[missing_pos] = reference_indices_for_phased[pos_prev]
		else:
			num_pos_lost_to_recombination += 1

	print "Assigned indices at %d positions." % len(reference_indices_for_unphased.keys())
	print "Lost %d positions to recombination." % num_pos_lost_to_recombination

	######################################################################
	# collect the phase at each unphased position

	print "Loading the pase for unphased positions..."

	# build set of indices to take from ref_panel
	indices = list()
	if args.leave_out:
		left_out_members = set(args.leave_out.strip().split(','))
	else:
		left_out_members = set()

	if args.chr in {'chr' + str(x) for x in xrange(1,25)}:
			sex_chromosome = False
	else:
	        sex_chromosome = True

	with open(args.reference_sample) as sample:
		headers = sample.readline()
		for i,line in enumerate(sample):
			sample_id, nationality, population, sex = line.strip().split()
			if args.population != 'ALL' and population != args.population: continue
			if sample_id in left_out_members: continue

			if not sex_chromosome or sex == 2:
				indices.append(2*i+1)
				indices.append(2*i+2)
			else:
				indices.append(2*i+1)

	cmd = "zcat %s | cut -d ' ' -f %s" % (args.reference_panel, ','.join([str(i) for i in indices]))
	ref_panel_hap = os.popen(cmd)


	legend = os.popen('zcat ' + args.reference_legend)
	headers = legend.readline()

	ref_values = dict()
	num_unphaseable = 0
	num_diff_from_ref = 0
	for ref_line, legend_line in zip(ref_panel_hap, legend):
		fields = legend_line.split()
		pos = int(fields[1])
		if pos in unphased_positions:
			# all positions in 1000g are 0/1, so if we see a 2, continue:
			if pos in one_two_positions:
				num_diff_from_ref += 1
				continue

			if pos not in letter_by_allele:
				continue

			# check if 0/1 corresponds to the same thing:
			if letter_by_allele[pos][0] != fields[2] or letter_by_allele[pos][1] != fields[3]:
				num_diff_from_ref += 1
				continue

			# if position is not phaseable, skip:
			if pos not in reference_indices_for_unphased:
				num_unphaseable += 1
				continue

			# if everything is okay,
			ref = ref_line.split()
			i0, i1 = reference_indices_for_unphased[pos]
			ref_values[pos] = int(ref[i0]), int(ref[i1])
	ref_panel_hap.close()
	legend.close()

	print "Assigned a phase to %d unphased positions." % len(ref_values.keys())
	print "Could not assign phase to %d positions because they had different alleles in the reference." % num_diff_from_ref
	print "Could not assign phase to %d positions because they did not have a reference index." % num_unphaseable

	######################################################################
	# assign a phase to each block:

	def metadata_to_str(metadata):
		return ':'.join(['%s=%s' % (str(k), str(v)) for k, v in metadata.iteritems()])

	def write_block(out, line, phase):
		fields = line.strip().split()
		metadata = parse_metadata(fields[3])
		metadata['TR'] = 0.992 # empirical average

		out.write('\t'.join(fields[0:3]))
		out.write('\t' + metadata_to_str(metadata))

		for f in fields[4:]:
			fi = f.split(':')

			pos = int(fi[0])
			p0 = (int(fi[1]) + phase) % 2
			p1 = (int(fi[2]) + phase) % 2

			out.write('\t%d:%d:%d' % (pos, p0, p1))

		out.write('\n')

	print "Writing phased output..."
	num_written = 0
	num_phased = 0
	phase_file = open(args.phased_blocks, 'w')
	for line, block in zip(block_lines, blocks):
		haplotype = np.array([ref_values[pos] for pos,p0,p1 in block if pos in ref_values])
		block_haplotype = np.array([(p0,p1) for pos,p0,p1 in block if pos in ref_values])

		if len(block_haplotype) == 0:
			# no snps were imputed in that block
			continue

		distance0 = np.sum(np.bitwise_xor(haplotype, block_haplotype))
		distance1 = np.sum(np.bitwise_xor(haplotype[:,(1,0)], block_haplotype))

		if 0.4 < (float(distance0) / float(distance0 + distance1)) < 0.6:
			# cannot say for sure where the block goes
			continue

		num_written += 1
		num_phased += len(block)
		if distance0 < distance1:
			write_block(phase_file, line, 0)
		else:
			write_block(phase_file, line, 1)

	print "%d/%d blocks were phased." % (num_written, len(blocks))
	print "%d/%d positions were phased." % (num_phased, len(unphased_positions))
	phase_file.close()

if __name__ == '__main__':
	main()