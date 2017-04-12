# RAILS v1.1 and Cobbler v0.2
## Rene L. Warren, 2014-2016
## email: rwarren at bcgsc.ca

### Name
------------

RAILS: Radial Assembly Improvement by Long Sequence Scaffolding
Cobbler: Gap-filling with long sequences


### Description
------------

RAILS and Cobbler are genomics application for scaffolding and automated finishing of genome assemblies with long DNA sequences.
They can be used to scaffold & finish high-quality draft genome assemblies with any long, preferably high-quality, sequences such as scaftigs/contigs from another genome draft. 

They both rely on accurate, long DNA sequences to patch gaps in existing genome assembly drafts.

Cobbler is a utility to automatically patch gaps (ambiguous regions in a draft assembly, represented by N's)
It does so by first aligning the long sequences to the assembly, tallying the alignments and replacing N's with the sequences from these long DNA sequences.

RAILS is an all-in-one scaffolder and gap-filler. Its process is similar to that of Cobbler. It scaffolds your genome draft with the help of long DNA sequences (contig sequences are ordered/oriented using alignment information). The newly created gaps are automatically filled with the DNA sequence of the provided long DNA sequence.

You can test the software by executing "runme.sh" in the test folder. A simulated SARS genome assembly is provided to test the software. 

### Implementation and requirements
------------

RAILS and Cobbler are implemented in PERL and run on any OS where PERL is installed.


### Community guidelines:
------------

I encourage the community to contribute to the development of this software, by providing suggestions for improving the code and/or directly contributing to the open source code for these tools. Users and developers may report software issues, bug fix requests, comments, etc, at <https://github.com/warrenlr/RAILS>


### Install
------------

Download the tar ball, gunzip and extract the files on your system using:
<pre>
gunzip rails_v1-1.tar.gz
tar -xvf rails_v1-1.tar
</pre>
Alternatively, individual tools are available within the github repository


### Dependencies
------------

Make sure you have installed bwa (Version: 0.7.15-r1140) and that is is in your path.


### Test data
------------

Go to ./test
(cd test)

1. SARS:
execute runme.sh
(./runme.sh)

2. Human:
execute runmeHuman.sh (will take a while to run)
(./runmeHuman.sh)

###Citing RAILS/Cobbler
------------

Thank you for using, developing and promoting this free software.
If you use RAILS or Cobbler for you research, please cite:

<pre>
Warren RL. 2016. RAILS and Cobbler: Scaffolding and automated finishing
of draft genomes using long DNA sequences. The Journal of Open Source
Software. doi: 10.21105/joss.00116
</pre>


### Usage
------------
<pre>
./runRAILS.sh
Usage: runRAILS.sh <FASTA assembly .fa> <FASTA long sequences .fa> <anchoring sequence length eg. 250> <min sequence identity 0.95>

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

Usage: ./cobbler.pl [v0.2]
-f  Assembled Sequences to further scaffold (Multi-Fasta format, required)
-q  Long Sequences queried (Multi-Fasta format, required)
-s  SAM file
-d  Anchoring bases on contig edges (ie. minimum required alignment size on contigs, default -d 1000, optional)
-i  Minimum sequence identity, default -i 0.9, optional
-t  LIST of names/header, long sequences to avoid using for merging/gap-filling scaffolds (optional)
-b  Base name for your output files (optional)
-v  Runs in verbose mode (-v 1 = yes, default = no, optional)

Usage: ./RAILS [v1.1]
-f  Assembled Sequences to further scaffold (Multi-Fasta format, required)
-q  Long Sequences queried (Multi-Fasta format, required)
-s  SAM file
-d  Anchoring bases on contig edges (ie. minimum required alignment size on contigs, default -d 1000, optional)
-i  Minimum sequence identity, default -i 0.9, optional
-t  LIST of names/header, long sequences to avoid using for merging/gap-filling scaffolds (optional)
-b  Base name for your output files (optional)
-v  Runs in verbose mode (-v 1 = yes, default = no, optional)
</pre>

### How it works
------------

The pipeline is detailed in the provided script runRAILS.sh

Cobbler's process:

The assembly draft sequence supplied to Cobbler is first broken up at the ambiguous regions of the assembly (Ns) to create scaftigs.
In the runRAILS.sh, these scaftigs are renamed, tracking their scaffold of origin (renumbered incrementally) and their position within it (also numbered incrementally).
A bwa index is created and the long sequence file, also re-numbered, is aligned to the scaftigs.
Cobbler is supplied with the alignment file (-s sam file) and the long reads files (-q option), specifying the minimum length of anchoring bases (-d) aligning at the edge of scaftigs and the minimum sequence identity of the alignment (-i). When 1 or more long sequences align unambiguously to the 3'end of a scaftig and the 5'end of its neighbour, the gap is patched with the sequence of that long sequence. If no long sequences are suitable, or the -d and -i conditions are not met, the original Ns are placed back between those scaftigs.

RAILS process:

In RAILS, the process is similar as for Cobbler, except that the draft assembly is not broken up at Ns, since the goal is to merge distinct sequences into larger ones.  Long sequences are aligned to the draft assembly sequences, orienting and ordering sequences and simulateneously filling the gaps between them, using DNA bases from the long sequences.

Scaffolding in RAILS is done using the LINKS scaffolder code (Warren et al. 2015), the unpublished scaffolding engine in the widely-used SSAKE assembler (Warren et al. 2007), and foundation of the SSPACE-LongRead scaffolder (Boetzer and Pirovano, 2014).

Output: For both Cobbler and RAILS, a summary of the gaps closed and their lengths is provided (.tsv) as a text file.
A fasta file (.fa) of the finished and/or scaffolded draft is generated for both along with a log file reporting basic success statistics.


Boetzer M, Pirovano W. 2014. SSPACE-LongRead: scaffolding bacterial draft genomes using long read sequence information. BMC Bioinformatics.15:211. DOI: 10.1186/1471-2105-15-211

Warren RL, Yang C, Vandervalk BP, Behsaz B, Lagman A, Jones SJ, Birol I. 2015. LINKS: Scalable, alignment-free scaffolding of draft genomes with long reads. GigaScience 4:35. DOI: 10.1186/s13742-015-0076-3

Warren RL, Sutton GG, Jones SJM, Holt RA.  2007.  Assembling millions of short DNA sequences using SSAKE.  Bioinformatics. 23(4):500-501. DOI: 10.1093/bioinformatics/btl629


### Runs on the human genome
------------

On a human draft assembly, cobbler patched over 65% of the gaps using 1, 2.5, 5, 15 kb long DNA sequences simulated from the human genome reference. The Pearson correlation between the predicted gap sizes and the size of patched gaps is R=0.8150


**Table 1.** Patching gaps in a genome assembly draft with Cobbler, using simulated 1, 2.5, 5 and 15 kbp simulated long sequences from human genome reference GRCh38.

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

RAILS was used to further contiguate the human baseline assembly draft and automatically close gaps within in:

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


### License Preamble
------------

RAILS and Cobbler Copyright (c) 2014-2016 British Columbia Cancer Agency Branch.  All rights reserved.

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

