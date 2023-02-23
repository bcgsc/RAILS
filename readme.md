[![Release](https://img.shields.io/github/release/bcgsc/RAILS.svg)](https://github.com/bcgsc/RAILS/releases)
[![Downloads](https://img.shields.io/github/downloads/bcgsc/RAILS/total?logo=github)](https://github.com/bcgsc/RAILS/releases/download/v1.5.1/rails_v1-5-1.tar.gz)
[![Issues](https://img.shields.io/github/issues/bcgsc/RAILS.svg)](https://github.com/bcgsc/RAILS/issues)
[![link](https://img.shields.io/badge/RAILScobbler-manuscript-brightgreen)](https://doi.org/10.21105/joss.00116)
Thank you for your [![Stars](https://img.shields.io/github/stars/bcgsc/RAILS.svg)](https://github.com/bcgsc/RAILS/stargazers)

![Logo](https://github.com/bcgsc/RAILS/blob/master/rails-logo.png)

# RAILS v1.5.1 and Cobbler v0.6.1
## Rene L. Warren, 2014-2023

### Contents
--------
1. [Name](#name)
2. [Description](#des)
3. [What's new](#new)
4. [Implementation and requirements](#imp)
5. [Community guidelines](#guide)
6. [Installation](#install)
7. [Dependencies](#dep)
8. [Test data](#test)
9. [Citing RAILS/Cobbler](#citing)
10. [Usage](#usage)
11. [Algorithm](#algo)
12. [Runs on human](#runs)
13. [License preamble](#license)
--------


### Name <a name=name></a>
-------------

<pre>
RAILS: Radial Assembly Improvement by Long Sequence Scaffolding

Cobbler: Gap-filling with long sequences
</pre>

### Description <a name=des></a>
-------------

RAILS and Cobbler are genomics application for scaffolding and automated finishing of genome assemblies with long DNA sequences.
They can be used to scaffold & finish high-quality draft genome assemblies with any long, preferably high-quality, sequences such as scaftigs/contigs from another genome draft. 

They both rely on accurate, long DNA sequences to patch gaps in existing genome assembly drafts.

Cobbler is a utility to automatically patch gaps (ambiguous regions in a draft assembly, represented by N's)
It does so by first aligning the long sequences to the assembly, tallying the alignments and replacing N's with the sequences from these long DNA sequences.

RAILS is an all-in-one scaffolder and gap-filler. Its process is similar to that of Cobbler. It scaffolds your genome draft with the help of long DNA sequences (contig sequences are ordered/oriented using alignment information). The newly created gaps are automatically filled with the DNA sequence of the provided long DNA sequence.

You can test the software by executing "runme.sh" in the test folder. A simulated SARS genome assembly is provided to test the software. 

### What's new in v1.5.1 <a name=new></a>

Remove requirement on samtools when running in "stream" mode


### What's new in v1.5.0

Ability to stream the .sam output of your favorite aligner directly into cobbler/RAILS (tested with minimap2/human data -- see runRAILSminimapSTREAM.sh)


### What's new in v1.4.2

Improved documentation, minor fixes, support for minimap2 (see runRAILSminimap.sh in the test folder)


### What's new in v1.4.1

1. Save in memory gap sequence from highest-matching read for both cobbler and RAILS
2. Track the number of reads support in cobbler (-l) and RAILS, and allow cutoff when scaffolding (-l and -a), with latter (RAILS)
3. Remove the hardcoded two-hit requirement for a read in RAILS. Instead, process two best hits for each read aligning different sequences
4. Implement grace (-g) option, which effectively simulate read trimming (valuable for Nanopore read mapping (suggested -g 250 to -g 500))
5. bug fixes (-list.tsv (cobbler) reported some instances of gap-fill regions not fixed in the assembly). cobbler gap-fill table now lists #supporting reads for each gap filled


### Implementation and requirements <a name=imp></a>
-------------

RAILS and Cobbler are implemented in PERL and run on any OS where PERL is installed.
Both tools require samtools (tested with v1.8) to read sequence alignment bamfiles. 
The runRAILS.sh pipeline requires bwa (see Dependencies below for tested version).
The runRAILSminimap.sh and runRAILSminimapSTREAM.sh pipelines require minimap2.
Please make sure these tools are in your PATH before running the above pipelines.


### Community guidelines <a name=guide></a>
-------------

I encourage the community to contribute to the development of this software, by providing suggestions for improving the code and/or directly contributing to the open source code for these tools. Users and developers may report software issues, bug fix requests, comments, etc, at <https://github.com/warrenlr/RAILS>


### Installation <a name=install></a>
-------------

Download the tar ball, gunzip and extract the files on your system using:

<pre>
gunzip rails_v1-5-1.tar.gz
tar -xvf rails_v1-5-1.tar
</pre>

Pleasure ensure that both cobbler.pl and RAILS are in your PATH.

Alternatively, individual tools are available for download/cloning within the github repository


### Dependencies <a name=dep></a>
-------------

Make sure you have installed bwa (Version: 0.7.15-r1140) or minimap2 (2.15-r905) and that they are in your PATH.
Make sure you have installed samtools (Version: 1.8) and that it is in your PATH.

Other versions of bwa, minimap2 & samtools may or may not be compatible and they have not been tested. Users may choose to use other versions than the ones specified here, at they see fit, but are expected to thoroughly test the behavior on their own.

Compatible tools may be used, but have not been tested fully (eg. sambamba)


### Test data <a name=test></a>
-------------

<pre>
Go to ./test
(cd test)

You may need to change both runme.sh and runmeHuman.sh to specify the path of samtools on your system

1. SARS:
execute runme.sh
(./runme.sh)

2. Human:
execute runmeHuman.sh (will take a while to run with bwa mem (~12h). With minimap2, this test will take ~1h.)
(./runmeHuman.sh)
</pre>


### Citing RAILS/Cobbler <a name=citing></a>
-------------

Thank you for your [![Stars](https://img.shields.io/github/stars/bcgsc/RAILS.svg)](https://github.com/bcgsc/RAILS/stargazers) and for using, developing and promoting this free software!

If you use RAILS or Cobbler for you research, please cite:

<pre>
Warren RL. 2016. RAILS and Cobbler: Scaffolding and automated finishing
of draft genomes using long DNA sequences. The Journal of Open Source
Software. doi: 10.21105/joss.00116
</pre>
[![link](https://img.shields.io/badge/RAILScobbler-manuscript-brightgreen)](https://doi.org/10.21105/joss.00116)


### Usage <a name=usage></a>
-------------

<pre>
./runRAILS.sh
Usage: runRAILS.sh <FASTA assembly .fa> <FASTA long sequences .fa> <anchoring sequence length eg. 250> <min sequence identity 0.95> <path to samtools>

this pipeline will:
1. reformat the assembly file $1
2. rename the long sequence file $2
3. Build a database index with bwa
4. Align the reformatted long sequences to your re-formatted baseline assembly
5. Run Cobbler to gap-fill regions of ambiguity
6. Reformat Cobbler's .fa file
7. Build a database index of it with bwa
8. Align the reformatted long sequences to your re-formatted cobbler assembly
9. Run RAILS to generate a newly scaffolded assembly draft

Usage: ./cobbler.pl [v0.6.1]
-f  Assembled Sequences to further scaffold (Multi-FASTA format NO LINE BREAKS, required)
-q  File of filenames containing long Sequences queried (Multi-FASTA format NO LINE BREAKS, required)
-s  File of filenames containing full path to BAM file(s) (use v0.2 for reading SAM files) or simply type: stream for streaming the .sam output of minimap2 or favorite aligner
-p  Full path to samtools (known to work/tested with v1.8, required if reading BAM files)
-d  Anchoring bases on contig edges (ie. minimum required alignment size on contigs, default -d 1000, optional)
-i  Minimum sequence identity fraction (0 to 1), default -i 0.9, optional
-l  Minimum number of long sequence support per gap, default -l 1, optional
-g  Grace length (bp), default -g 1, optional
-t  LIST of names/header, long sequences to avoid using for merging/gap-filling scaffolds (optional)
-b  Base name for your output files (optional)
-v  Runs in verbose mode (-v 1 = yes, default = no, optional)
IMPORTANT: the order of files in -q and -s MUST match!


Usage: ./RAILS [v1.5.1]
-f  Assembled Sequences to further scaffold (Multi-Fasta format, required)
-q  File of filenames containing long Sequences queried (Multi-Fasta format, required)
-s  File of filenames containing full path to BAM file(s) or simply type: stream for streaming the .sam output of minimap2 or favorite aligner
-p  Full path to samtools (known to work/tested with v1.8, required if reading BAM files)
-d  Anchoring bases on contig edges (ie. minimum required alignment size on contigs, default -d 1000, optional)
-i  Minimum sequence identity fraction (0 to 1), default -i 0.9, optional
-t  LIST of names/header, long sequences to avoid using for merging/gap-filling scaffolds (optional)
-l  Minimum number of links to compute scaffold (default -l 1, optional)
-a  Maximum link ratio between two best contig pairs *higher values lead to least accurate scaffolding* (default -a 0.99, optional)
-g  Grace length (bp), default -g 1, optional
-b  Base name for your output files (optional)
-v  Runs in verbose mode (-v 1 = yes, default = no, optional)
IMPORTANT: the order of files in -q and -s MUST match!


</pre>

### Algorithm <a name=algo></a>
-------------

The pipeline is detailed in the provided script runRAILS.sh. PLEASE ensure the draft assembly is FASTA-formatted with one sequence per line (NO LINE BREAKS)

Cobbler's process:

The assembly draft sequence supplied to Cobbler is first broken up at the ambiguous regions of the assembly (Ns) to create scaftigs.
In the runRAILS.sh, these scaftigs are renamed, tracking their scaffold of origin (renumbered incrementally) and their position within it (also numbered incrementally).
A bwa index is created and the long sequence file, also re-numbered, is aligned to the scaftigs.
Cobbler is supplied with the alignment file (-s sam file) and the long reads files (-q option), specifying the minimum length of anchoring bases (-d) aligning at the edge of scaftigs and the minimum sequence identity of the alignment (-i). When 1 or more long sequences align unambiguously to the 3'end of a scaftig and the 5'end of its neighbour, the gap is patched with the sequence of that long sequence. If no long sequences are suitable, or the -d and -i conditions are not met, the original Ns are placed back between those scaftigs.

RAILS process:

In RAILS, the process is similar as for Cobbler, except that the draft assembly is not broken up at Ns, since the goal is to merge distinct sequences into larger ones.  Long sequences are aligned to the draft assembly sequences, orienting and ordering sequences and simulateneously filling the gaps between them, using DNA bases from the long sequences.

Scaffolding in RAILS is done using the LINKS scaffolder code (Warren et al. 2015), the unpublished scaffolding engine in the widely-used SSAKE assembler (Warren et al. 2007), and foundation of the SSPACE-LongRead scaffolder (Boetzer and Pirovano, 2014).

The grace (-g) parameter may be used to set the MAXIMUM length of unaligned bases allowed at the end of each (long) sequencing read alignment to the draft genome assembly. For example, setting -g 250 tells cobbler/RAILS to consider a sequencing read with a soft-clip of up to 250 bp in 5' or 3' 

Output: For both Cobbler and RAILS, a summary of the gaps closed and their lengths is provided (.tsv) as a text file.
A fasta file (.fa) of the finished and/or scaffolded draft is generated for both along with a log file reporting basic success statistics.

<pre>
Boetzer M, Pirovano W. 2014. SSPACE-LongRead: scaffolding bacterial draft genomes using long read sequence information. BMC Bioinformatics.15:211. DOI: 10.1186/1471-2105-15-211

Warren RL, Yang C, Vandervalk BP, Behsaz B, Lagman A, Jones SJ, Birol I. 2015. LINKS: Scalable, alignment-free scaffolding of draft genomes with long reads. GigaScience 4:35. DOI: 10.1186/s13742-015-0076-3

Warren RL, Sutton GG, Jones SJM, Holt RA.  2007.  Assembling millions of short DNA sequences using SSAKE.  Bioinformatics. 23(4):500-501. DOI: 10.1093/bioinformatics/btl629
</pre>


### Runs on human <a name=runs></a>
-------------

On a human HG004 ABySS draft assembly, cobbler filled over 65% of the gaps using 1, 2.5, 5, 15 kb long DNA sequences simulated from the human genome reference. The Pearson correlation between the predicted gap sizes and the size of patched gaps is R=0.8150


**Table 1.** Patching gaps with Cobbler (v0.2) using simulated 1, 2.5, 5, 15kbp simulated long sequences from human genome reference GRCh38.

Metric | Value
---- | ----:
Total gaps | 148,091
Number of gaps patched | 95,523
Proportion of gaps patched | 65.1%
Average length (bp) | 343.39
Length st.dev +/- | 931.12
Total bases added | 32,801,755
Largest gap resolved (bp) | 13,662
Shortest gap resolved (bp) | 1 

RAILS (v1.1) was used to further contiguate the human baseline assembly draft and automatically close gaps within in:

**Table 2.** RAILS scaffolding and gap-filling summary on a human assembly baseline, using simulated 1, 2.5, 5, 15kbp simulated long sequences from human genome reference GRCh38.

Metric | Value
---- | ----:
Number of merges induced | 6,029
Average closed gap length (bp) | 1,136.71
Closed gap length st.dev +/- | 2,511.69
Total bases added | 6,853,222
Largest gap resolved (bp) | 14,471
Shortest gap resolved (bp) | 1

6,029 merges resulted from RAILS scaffolding of the baseline human assembly draft (1,695 >= 500bp)
The scaffold N50 length increased from 5.6 to 7.3 Mbp, a 30% increase in N50 length.


**Table 3.** Assembly statistics on human genome scaffolding and finishing post Cobbler and RAILS (reporting sequences 500 bp and larger).

Stage | n:500 | n:N50 | n:NG50 |        NG50 (bp) |     N50 (bp)|       max (bp) |      sum (bp)
--------- | ------: | -----: | -----: | ---------: | ---------: | ---------: |  -------:
Baseline | 65,905 |        145 |   164 |   5,144,025 |     5,597,244 |       26.41e6 |       2.794e9
Cobbler | 65,905 |        145 |   161 |   5,312,196 |     5,658,133 |       26.66e6 |       2.827e9
RAILS | 64,210 |        113 |   125 |   6,935,685 |     7,266,542 |       32.14e6 |       2.836e9


### License preamble <a name=license></a>
-------------

RAILS and Cobbler Copyright (c) 2014-2023 British Columbia Cancer Agency Branch.  All rights reserved.

RAILS and Cobbler are released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

