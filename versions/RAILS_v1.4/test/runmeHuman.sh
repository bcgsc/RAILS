##Example of automated finishing (cobbler) and scaffolding (RAILS) done on HG004 genome in a bottle human genome assembly
#Downloading the draft assembly
wget ftp://ftp.bcgsc.ca/supplementary/RAILS/hsapiens-8.fa.gz
#Dowloading the simulated long sequences, derived from the human genome reference (1, 2.5, 5, 15kbp reads)
wget ftp://ftp.bcgsc.ca/supplementary/RAILS/reads.fa.gz
#Decompressing files
gunzip hsapiens-8.fa.gz
gunzip reads.fa.gz
#Running Cobbler and RAILS, Make sure you have bwa in your path
./runRAILS.sh hsapiens-8.fa reads.fa 250 0.95
