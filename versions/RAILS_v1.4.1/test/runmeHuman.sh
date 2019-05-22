#!/bin/bash
#RLW 2016,2019
#-------------
##Example of automated finishing (cobbler) and scaffolding (RAILS) of HG004 genome-in-a-bottle human genome assembly
#-------------
echo Downloading the draft assembly...
wget http://www.bcgsc.ca/downloads/supplementary/RAILS/hsapiens-8.fa.gz
#Dowloading the simulated long sequences, derived from the human genome reference (1, 2.5, 5, 15kbp reads)
echo Dowloading simulated long sequences...
wget http://www.bcgsc.ca/downloads/supplementary/RAILS/reads.fa.gz
#-------------
echo Decompressing files...
gunzip hsapiens-8.fa.gz
gunzip reads.fa.gz
#-------------
echo Running Cobbler and RAILS, Make sure you have bwa or minimap2 in your path...
#./runRAILS.sh hsapiens-8.fa reads.fa 250 0.95 /gsc/btl/linuxbrew/bin/samtools
echo Running the minimap2 pipeline. To run the bwa pipeline uncomment line above in the script
./runRAILSminimap.sh hsapiens-8.fa reads.fa 250 0.95 0 1 nil /gsc/btl/linuxbrew/bin/samtools
echo Process completed.
