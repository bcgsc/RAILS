#!/bin/bash
#RLW 2016,2019

if [ $# -ne 9 ]; then
        echo "Usage: $(basename $0) <FASTA assembly .fa> <FASTA long sequences .fa> <anchoring sequence length eg. 250> <min sequence identity 0.95> <max. softclip eg. 250bp> <min. number of read support eg. 2> <long read type eg.: ont, pacbio, nil> <path to samtools> <number of threads>"
        exit 1
fi

echo Resolving ambiguous bases -Ns- in $1 assembly using long sequences $2
#-------------------------
echo reformatting file $1
### WARNING: MAKE SURE YOUR INPUT FASTA IS ONE SEQUENCE PER LINE, WITH NO LINE BREAKS!
echo WARNING: MAKE SURE YOUR INPUT FASTA IS ONE SEQUENCE PER LINE WITH NO LINE BREAKS!
cat $1 | perl -ne 'if(/^\>/){$scafnum++;}else{my $len=length($_);my @scaftigs=split(/N+/i,$_);my $scaftignum=0;foreach my $scaftig(@scaftigs){ my $len=length($scaftig);$scaftignum++;  print ">wga$scafnum";print "."; print "$scaftignum,$len\n$scaftig\n";}}' > $1-formatted.fa
echo reformatting file $2
cat $2 | perl -ne 'if(/^\>/){$ct++;}else{my $len=length($_);print ">seq$ct,$len\n$_";}' > $2-formatted.fa
#--------------------------
echo Aligning long sequences $2-formatted.fa to your contigs..

if [ $7 == 'ont' ]; then
   echo Running minimap2 with preset map-ont
   minimap2 -x map-ont -I50g -N 10 -a -t $9 $1-formatted.fa $2-formatted.fa | samtools view -Sb - > $2_vs_$1_gapfilling.bam

elif [ $7 == 'pacbio' ]; then
   echo Running minimap2 with preset map-pb
   minimap2 -x map-pb -I50g -N 10 -a -t $9 $1-formatted.fa $2-formatted.fa | samtools view -Sb - > $2_vs_$1_gapfilling.bam

else
   echo Running minimap2 with no preset
   minimap2 -I50g -N 10 -a -t $9 $1-formatted.fa $2-formatted.fa | samtools view -Sb - > $2_vs_$1_gapfilling.bam

fi

#--------------------------
echo Gap-filling $1-formatted.fa using $2-formatted.fa
echo $2-formatted.fa > $2-formatted.fof
echo $2_vs_$1_gapfilling.bam > $2_vs_$1_gapfilling.fof
echo Running cobbler.pl -f $1 -s $2_vs_$1_gapfilling.fof -l $6 -g $5 -d $3 -i $4 -b $2_vs_$1_$3_$4_gapsFill -q $2-formatted.fof -p $8 ...
./cobbler.pl -f $1 -s $2_vs_$1_gapfilling.fof -l $6 -g $5 -d $3 -i $4 -b $2_vs_$1_$3_$4_gapsFill -q $2-formatted.fof -p $8
echo Process terminated.
#--------------------------
echo RAILS scaffolding $1.gapsFill.fa sequences and gap-filling using long seqs $2 -- anchoring sequence threshold $3 bp 
echo reformatting file $1.gapsFill.fa
cat $2_vs_$1_$3_$4_gapsFill.fa | perl -ne 'if(/^\>/){$ct++;}else{my $len=length($_);print ">wga$ct,$len\n$_";}' > $2_vs_$1_$3_$4_gapsFill-formatted.fa
#--------------------------
echo Aligning long sequences $2-formatted.fa to your contigs..

if [ $7 == 'ont' ]; then
   echo Running minimap2 with preset map-ont
   minimap2 -x map-ont -I50g -N 10 -a -t $9 $2_vs_$1_$3_$4_gapsFill-formatted.fa $2-formatted.fa | samtools view -Sb - > $2_vs_$1_scaffolding.bam

elif [ $7 == 'pacbio' ]; then
   echo Running minimap2 with preset map-pb
   minimap2 -x map-pb -I50g -N 10 -a -t $9 $2_vs_$1_$3_$4_gapsFill-formatted.fa $2-formatted.fa | samtools view -Sb - > $2_vs_$1_scaffolding.bam

else
   echo Running minimap2 with no preset
   minimap2 -I50g -N 10 -a -t $9 $2_vs_$1_$3_$4_gapsFill-formatted.fa $2-formatted.fa | samtools view -Sb - > $2_vs_$1_scaffolding.bam

fi

#--------------------------
echo Scaffolding $2_vs_$1_$3_$4_gapsFill-formatted.fa using $2-formatted.fa and filling new gaps with sequences in $2-formatted.fa
echo $2-formatted.fa > $2-formatted.fof
echo $2_vs_$1_scaffolding.bam > $2_vs_$1_scaffolding.fof
echo Running RAILS -f $2_vs_$1_$3_$4_gapsFill-formatted.fa -s $2_vs_$1_scaffolding.fof -l $6 -g $5 -d $3 -i $4 -b $2_vs_$1_$3_$4_rails -q $2-formatted.fof -p $8 ...
./RAILS -f $2_vs_$1_$3_$4_gapsFill-formatted.fa -s $2_vs_$1_scaffolding.fof -l $6 -g $5 -d $3 -i $4 -b $2_vs_$1_$3_$4_rails -q $2-formatted.fof -p $8
#--------------------------
echo RAILS process terminated.
