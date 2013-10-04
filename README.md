Prism
=====

Prism is an implementation of the statistical phasing algorithm used in

> Kuleshov, Xie, Chen, et al., Cost-effective whole-genome haplotyping using a combination of dilution and statistical methods.

Whole-genome haplotyping (also known as genome phasing) is the problem of determining the differences between the paternal and maternal copies of a chromosome. Current sequencing methods are presently unable to adequately resolve these differences; in practice, one uses either statistical methods, which can approximately infer the most likely genomic phase by looking at patterns in the genome of many individuals, or molecular approaches, which are often expensive and time consuming.

Prism uses a new phasing algorithm that augments traditional molecular haplotyping approaches with statistical methods. It determines haplotypes statistically where typical molecular phasing methods fail. By using Prism, one can substantially reduce the amount of molecular phasing that needs to be performed, and thus substantially reduce haplotyping costs. The paper by Kuleshov, Xie, Chen, et al. demonstrates how Prism can be combined with an inexpensive dilution haplotyping method to produce haplotypes that match the quality of existing methods at a fraction of the cost. The resulting technology, statistically-aided dilution haplotyping (SDH) represents the first commercially available haplotyping product and is offered by Illumina Inc.

This document contains instructions for running Prism, as well as a step-by-step tutorial that goes over running Prism on a 1 Mbp region of chromosome 22.

Copyright
---------
Prism has been written by Volodymyr Kuleshov. The algorithm and source code are the property of Illumina Inc.

The source code of Prism is released under the Illumina Open Source License, available at https://github.com/sequencing/licenses.

Installation
------------

A package containing Prism is available at http://www.stanford.edu/~kuleshov/prism.tar.gz. To install the package, extract the contents, and run the setup.py script:

wget http://www.stanford.edu/~kuleshov/prism.tar.gz
tar -zxvf prism.tar.gz
cd prism
python setup.py install

Prism is written in Python and C, and requires numpy (>=1.5.1) and scipy (>=0.8.0) as well as a recent version of gcc. It makes use of the SciPy Weave package (http://docs.scipy.org/doc/scipy/reference/tutorial/weave.html) to compile inline C code into Python. It has been tested on a stock Ubuntu 12.04 Linux machine.

Instructions
------------

Prism takes as input a set of haplotype blocks assembled through a form of dilution haplotyping (which we call "local" blocks) and outputs for each block its most likely phase relative to its predecessor, as well as a confidence score. Given a desired accuracy threshold, one can combine the local blocks into global statistically-assembled superblocks. When used with the SDH technology described in Kuleshov, Xie, Chen, et al, these global blocks attain N50 lengths of 500 Kbp or more and contain on average fewer than one long switching error per Mbp.

Prism phases locally-assembled blocks based on a phased reference panel. It was designed to be used with the Thousand Genomes project pre-phased reference panel that is distributed as part of the IMPUTE2 software package.

Statistical phasing is performed as a two-stage process, each associated with one of two commands: `prism` and `prism-interleaving`. The `prism` command computes the phase for a subset of the local blocks that do not overlap among each other. It produces also a _path_ file, which contains a "path" through the phased reference panel, i.e. which describes the pair of reference panel haplotypes that most closely matches the subject at each genomic position. The `prism-interleaving`, then determines the phase of the remaining blocks based on the _path_ file.

More precisely, the `prism` command works as follows:

prism
--blocks <local-blocks>
--positions <reference panel at position that need to be phased>
--start <genomic coordinate>
--end <genomic coordinate>
--phase <output file containing the phase of each block>
--path <output file describing the path through the reference panel>
--K <number of reference haplotypes to use in current window>
--N <effective population size parameter>

Each input block must be tagged with an ID. The phase output file has the following format: 

ID PHASE SCORE 

where PHASE equals 0 or 1. If it equals 1, the block identified with ID should be flipped. The path output file uses the format

POS R1 R1 A1 A2

where POS is a genomic position from the positions file, R1 and R2 are the indices of the reference panel that best describe the sample at that position, and A1, A2 are the allele of the reference panel haplotypes at that position.

The `prism-interleaving` command requires the following input

prism-interleaving
--unphased-blocks <set of interleaving blocks to be phased>
--path <output file describing the path through the reference panel>
--vcf <vcf file of the subject>
--chr <chromosome to consider>
--reference-panel <the reference panel for the chromosome; this requires a .hap.gz file in IMPUTE2 format>
--reference-legend <the reference panel for the chromosome; this requires a .legend.gz file in IMPUTE2 format>
--reference-sample <the reference panel for the chromosome; this requires a .sample file in IMPUTE2 format>
--phased-blocks <output file containing output blocks>

The phased blocks returned by `prism-interleaving` are in the same format as the input blocks, but some of them have been flipped to match the phase that was deemed the most likely. The algorithm used behind `prism-interleaving` is very simple; it simply assigns to each block the phase that is closest to the haplotypes defined in the "path" input file.

A set of locally-phased blocks is transformed into global blocks in two stages. Initially, the blocks must be split into two subsets, in a way that the first subset contains no overlapping blocks.

Demonstration
-------------

The following is a step-by-step tutorial on how to use Prism over a small ~1Mbp region of chromosome 22 of HapMap sample NA12878. We ask the reader to run these steps on a stock Ubuntu 12.04 Linux machine, and recommend using an m1.large Amazon EC2 instance, so that there are no issues with performance. Depending on the chosen parameter settings, Prism can use several gigabytes of RAM over the course of its execution.

### Setting up

On a stock Ubuntu 12.04 machine, extract the contents of the prism.tar.gz package:

tar -zxvf prism.tar.gz

### Phasing a set of local blocks

### Examining the accuracy

### Phasing the interleaving local blocks

### 

### Next steps