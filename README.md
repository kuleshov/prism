Prism
=====

Prism is an implementation of the statistical phasing algorithm used in

> Kuleshov, Xie, Chen, et al., Cost-effective whole-genome haplotyping using a combination of dilution and statistical methods.

Whole-genome haplotyping (also known as genome phasing) is the problem of determining the differences between the paternal and maternal copies of a chromosome. Current sequencing methods are presently unable to adequately resolve these differences; in practice, one uses either statistical methods, which can approximately infer the most likely genomic phase by looking at patterns in the genome of many individuals, or molecular approaches, which are often expensive and time consuming.

Prism uses a new phasing algorithm that augments traditional molecular haplotyping approaches with statistical methods. It determines haplotypes statistically where typical molecular phasing methods fail. By using Prism, one can substantially reduce the amount of molecular phasing that needs to be performed, and thus substantially reduce haplotyping costs. The paper by Kuleshov, Xie, Chen, et al. demonstrates how Prism can be combined with an inexpensive dilution haplotyping method to produce haplotypes that match the quality of existing methods at a fraction of the cost. The resulting technology, statistically-aided dilution haplotyping (SDH) represents the first commercially available haplotyping product and is offered by Illumina Inc.

This document contains instructions for running Prism, as well as a step-by-step tutorial that goes over running Prism on a 1 Mbp region of chromosome 22.

Copyright
---------
The Prism phaser was written by Volodymyr Kuleshov. The algorithm and source code are the property of Illumina Inc.

The source code of Prism is released under the Illumina Open Source License, available at https://github.com/sequencing/licenses.

Installation
------------

A package containing Prism is available at http://www.stanford.edu/~kuleshov/prism.tar.gz. To install the package, extract the contents, and run the setup.py script:

	wget http://www.stanford.edu/~kuleshov/prism.tar.gz
	tar -zxvf prism.tar.gz
	cd prism
	python setup.py install

Prism is written in Python and C, and requires numpy (>=1.5.1) and scipy (>=0.8.0) as well as a relatively recent version of gcc. It makes use of the SciPy Weave package to compile inline C code into Python. It has been tested on a stock Ubuntu 12.04 Linux machine.

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

where `PHASE` equals 0 or 1. If it equals 1, the block identified with ID should be flipped. The path output file uses the format

	POS R1 R1 A1 A2

where `POS` is a genomic position from the positions file, `R1` and `R2` are the indices of the reference panel that best describe the sample at that position, and `A1`, `A2` are the allele of the reference panel haplotypes at that position.

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

The following is a step-by-step tutorial on how to use Prism over a small ~1Mbp region of chromosome 22 of HapMap sample NA12878. We suggest that the reader run these steps on a stock Ubuntu 12.04 Linux machine, for example on an new m1.large Amazon EC2 instance. Depending on the chosen parameter settings, Prism can use several gigabytes of RAM over the course of its execution.

### Setting up

We recommend manually setting up the packages required for running Prism:

	sudo apt-get update
	sudo apt-get install gcc
	sudo apt-get install python-setuptools
	sudo apt-get install python-numpy
	sudo apt-get install python-scipy

On a stock Ubuntu 12.04 machine, extract the contents of the prism.tar.gz package and run the installation script:

	tar -zxvf prism.tar.gz
	cd prism
	sudo python setup.py install

This will make availble two new commands: `prism` and `prism-interleaving`.

Next, download the data package with various input files for NA12878 chromosome 22 (for simplicity, we recommend downloading directly to the prism folder):

	wget http://www.stanford.edu/~kuleshov/prism-data.tar.gz
	tar -zxvf prism-data.tar.gz

The `prism-data` package contains the following files:
* `na12878.chr22.local.blocks`: Set of non-interleaving blocks obtained from the local phasing stage of SDH.
* `na12878.chr22.local.interleaving.blocks`: Set of interleaving blocks from the local phasing stage of SDH. Together with the above file, these represent all the local blocks for chromosome 22. Note that there are many more interleaving blocks; however, they are all very short. They mostly represent single SNPs that fell within a larger local block, but could not be connected to it via a long fragment.
* `na12878.chr22.positions`: File representing the set of all positions to be considered for global phasing. Each position also contains its genetic distance from the very first position in the file, as well as the haplotypes of the reference panel.
* `na12878.chr22.vcf`: The VCF of the subject.
* `na12878.chr22.true-phase`: The true phase of the subject, determined by trio phasing.
* `impute2-panel/`: Data for chromosome 22 from the IMPUTE2 phasing panel.

### Phasing a set of local blocks

The first step in a typical Prism workflow would be to phase the non-interleaving blocks in `na12878.chr22.local.blocks`. These are long and accurate blocks obtained by local phasing. We may assess their accuracy using the `evaluate_blocks.py` helper script in the `scripts/` subfolder:

	python scripts/evaluate_blocks.py 
		--blocks prism-data/na12878.chr22.local.blocks 
		--true-phase prism-data/na12878.chr22.true-phase

The long switch accuracy for these blocks should be 99.95%. Unfortunately, these blocks are relatively short (N50 length of about 60 Kbp); let's now use the `prism` phaser to compute the most likely phase of each block relative to the other within a small genomic region: `chr22:18890344-19737410`.

	prism 
		--blocks prism-data/na12878.chr22.local.blocks 
		--positions prism-data/na12878.chr22.positions 
		--start 18890344 
		--end 19737410 
		--K 75 
		--path region.path 
		--phase region.phase

We manually specified a `--K` value of 75, indicating that `prism` should perform inference with the 75 haplotypes that are the closest to the subject within that particular region. Although larger values of `K` may yield better performance, the running time and memory usage of the algorithm increases superlinearly with `K`. A value of 75 should be suitable for a demonstration, but if the phaser takes too much time, we suggest moving to a value of 50 or 40.

Over the course of the execution, `prism` will automatically compile its C functions using SciPy Weave.

On a m1.large Amazon EC2 instance, the above command takes about 5 minutes to finish with `--K 75`.

After `prism` finishes successfully, it will have created two files: `region.phase` and `region.path`. For now, we may disregard the latter; the phase inferred for every block (identified by its ID) is found in the first file. To apply this phase to our local blocks, we use the `phase_local_blocks.py` helper script, which will flip the blocks the way the phaser specified:

	python scripts/phase_local_blocks.py 
		--phase region.phase 
		--unphased-local-blocks prism-data/na12878.chr22.local.blocks 
		--phased-local-blocks phased.local.blocks

We may assess the accuracy of the resulting statistical phasing using the `evaluate_statistical_accuracy.py` helper script:

	python scripts/evaluate_statistical_accuracy.py 
		--true-phase prism-data/na12878.chr22.true-phase 
		--blocks phased.local.blocks 
		--threshold 0

We should see a switch accuracy of 94.8%. This is slightly higher than we would normally expect for non-overlapping blocks. Typically, such blocks will be relatively far apart from each other, and will exhibit less linkage disequilibrium (LD) than, say, interleaving blocks. However, the region we chose seems to have somewhat higher LD than usual.

In addition to reporting the phase of each block, `prism` also provided us with confidence scores for each block. We may use these scores to disregard the suggested phase of some blocks. For instance, we may look only at the phase between blocks that have a confidence score above 0.8:

	python scripts/evaluate_statistical_accuracy.py 
		--true-phase prism-data/na12878.chr22.true-phase 
		--blocks phased.local.blocks 
		--threshold 0.8

This should give us a much higher confidence score: 97.4%. Typically, across the entire genome, the statistical phasing accuracy at this stage should be slightly above 96%.

### Phasing the interleaving local blocks

Next, we need to phase the remaining local blocks: ones in `na12878.chr22.local.interleaving.blocks`. For that we will need the `region.path` file generated at the previous stage, as well as the full Thousand Genomes reference panel for chromosome 22. The `region.path` specifies which pair of haplotypes from the reference panel best describe the subject at each genomic position that was considered by the `prism` command (i.e. the ones in the file `na12878.chr22.positions`). To determine the phase of each interleaving block, we will compare it to the haplotypes specified by these two files, and pick the phase that best matches them:

	prism-interleaving 
		--unphased-blocks prism-data/na12878.chr22.local.interleaving.blocks 
		--path region.path 
		--vcf prism-data/na12878.chr22.vcf 
		--chr chr22 
		--reference-panel prism-data/impute2-panel/ALL_1000G_phase1integrated_v3_chr22_impute.hap.gz 
		--reference-legend prism-data/impute2-panel/ALL_1000G_phase1integrated_v3_chr22_impute.legend.gz 
		--reference-sample prism-data/impute2-panel/ALL_1000G_phase1integrated_v3.sample 
		--phased-blocks phased.interleaving.blocks

This command will take about a minute to complete. It will produce a file in which the interleaving blocks whose phase could have been determined; they will be already flipped accordingly.

Let's now combine them with the blocks obtained at the previous stage:

	cat phased.local.blocks phased.interleaving.blocks | sort -n -k 3 > phased.blocks

We may now estimate again the statistical accuracy of our phasing:

	python scripts/evaluate_statistical_accuracy.py 
		--true-phase prism-data/na12878.chr22.true-phase 
		--blocks phased.blocks 
		--threshold 0.8

The statistical phasing accuracy should have now improved to 98%. This is typically what we would expect at this stage, using a threshold of 0.8. Interleaving blocks have a high statistical phasing accuracy because they are physically closer to neighboring blocks.

### Generating global blocks

We now have assigned a phase and a confidence score to every block obtained from the local stage. We may use this information to combine these haplotypes into global super-blocks. We will do so by going from left to right across the genome and adding the next local block to the current global block if its confidence score is above our minimum threshold. If it is not, we stop the current block, and start a new one. We recommend using a threshold of 0.8 for a good tradeoff between block size and quality.

	python scripts/get_global_blocks.py 
		--phased-local-blocks phased.blocks 
		--global-blocks global.blocks 
		--threshold 0.8

We can now evaluate the accuracy of these blocks:

	python scripts/evaluate_accuracy.py 
		--true-phase prism-data/na12878.chr22.true-phase 
		--blocks global.blocks

This should give a long switch accuracy of 99.8%, which is low by our standards, suggesting that an error slipped in somewhere in this region. A look at global.blocks shows we have two blocks of 250 and 600 Kbp:

cut -f 2,3 global.blocks | awk '{print $2 - $1 + 1}'

A more closer inspection would show that the second block has one long switch error in it, which accounts for the 0.2% loss in accuracy:

	1111111111111101111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
	00000000000000000000000000000000000000000000000000000000000000011111111111111111111111011111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111

Strangely, this error was not detected through the confidence scores. A closer look at the local blocks can explain why:

	49 	0000000000000000000000000000000
	50 	000000000
	51 	00000001111111111
	52 	1111111111
	53 	111
	54 	011111111
	55 	111111111111111111111111
	56 	111111111111111111111111111111111111
	57 	11
	58 	11111111111111111

It turns out that a switch error occurred in the middle of local block #51. Thus, the statistical algorithm actually made a sensible decision, as it aligned blocks 52-58 to continue the phase that was started in block 51.

Fortunately, local blocks have a high accuracy (recall that we measured it to be 99.95%) and errors as in the above case are quite rare. Most genomic windows, such as `chr22:24700000-25500000` or `chr22:20300000-21000000` exhibit no global errors. Across the chromsome 22, the long switch accuracy of the global blocks equals 99.9%.

### Next steps

The above tutorial showed how phase a 1Mbp window of chromosome 22. Extending this approach to the whole genome is conceptually simple, and requires partitioning the genome into many such overlapping windows. These windows can then all be phased in parallel on a compute cluster, and merged together to obtain an entire phased genome. For greater accuracy, one may use a larger panel size, although we did not see any improvements in accuracy beyond `-K 200`. Finally, it may be interesting to play with the `--N` parameter (effective population size in the Li & Stephens (2002) model), which may lead to a better accuracy on some samples, especially one that are very heterozygous.