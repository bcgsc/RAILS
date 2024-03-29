#!/bin/bash
#RLW 2016
if [ $# -ne 6 ]; then
        echo "Usage: $(basename $0) <FASTA assembly .fa> <FASTA long sequences .fa> <anchoring sequence length eg. 250> <min sequence identity 0.95> <path to samtools> <number of threads>"
        exit 1
fi

echo Resolving ambiguous bases -Ns- in $1 assembly using long sequences $2
echo reformatting file $1
### WARNING: MAKE SURE YOUR INPUT FASTA IS ONE SEQUENCE PER LINE, WITH NO LINE BREAKS!
echo WARNING: MAKE SURE YOUR INPUT FASTA IS ONE SEQUENCE PER LINE WITH NO LINE BREAKS!
cat $1 | perl -ne 'if(/^\>/){$scafnum++;}else{my $len=length($_);my @scaftigs=split(/N+/i,$_);my $scaftignum=0;foreach my $scaftig(@scaftigs){ my $len=length($scaftig);$scaftignum++;  print ">wga$scafnum";print "."; print "$scaftignum,$len\n$scaftig\n";}}' > $1-formatted.fa
echo reformatting file $2
cat $2 | perl -ne 'if(/^\>/){$ct++;}else{my $len=length($_);print ">seq$ct,$len\n$_";}' > $2-formatted.fa
echo Building sequence database index out of your $1-formatted.fa assembly contigs..
bwa index $1-formatted.fa
echo Aligning long sequences $2-formatted.fa to your contigs..
### YOU MAY CONSIDER: SETTING THE MORE STRINGENT bwa mem -x intractg OPTION AND ADJUSTING -t to higher values for speed
bwa mem -a -t $6 $1-formatted.fa $2-formatted.fa | samtools view -Sb - > $2_vs_$1_gapfilling.bam
echo Scaffolding $1-formatted.fa using $2-formatted.fa and filling gaps with sequences in $2-formatted.fa
echo $2-formatted.fa > $2-formatted.fof
echo $2_vs_$1_gapfilling.bam > $2_vs_$1_gapfilling.fof
cobbler.pl -f $1 -s $2_vs_$1_gapfilling.fof -d $3 -i $4 -b $2_vs_$1_$3_$4_gapsFill -q $2-formatted.fof -p $5
echo Process terminated.
echo RAILS scaffolding $1.gapsFill.fa sequences using long seqs $2 -- anchoring sequence threshold $3 bp 
echo reformatting file $1.gapsFill.fa
cat $2_vs_$1_$3_$4_gapsFill.fa | perl -ne 'if(/^\>/){$ct++;}else{my $len=length($_);print ">wga$ct,$len\n$_";}' > $2_vs_$1_$3_$4_gapsFill-formatted.fa
echo Building sequence database index out of your $2_vs_$1_$3_$4_gapsFill-formatted.fa assembly contigs..
bwa index $2_vs_$1_$3_$4_gapsFill-formatted.fa
echo Aligning long sequences $2-formatted.fa to your contigs..
### YOU MAY CONSIDER: SETTING THE MORE STRINGENT bwa mem -x intractg OPTION AND ADJUSTING -t to higher values for speed
bwa mem -a -t $6 $2_vs_$1_$3_$4_gapsFill-formatted.fa $2-formatted.fa | samtools view -Sb - > $2_vs_$1_scaffolding.bam
echo Scaffolding $2_vs_$1_$3_$4_gapsFill-formatted.fa using $2-formatted.fa and filling new gaps with sequences in $2-formatted.fa
echo $2-formatted.fa > $2-formatted.fof
echo $2_vs_$1_scaffolding.bam > $2_vs_$1_scaffolding.fof
RAILS -f $2_vs_$1_$3_$4_gapsFill-formatted.fa -s $2_vs_$1_scaffolding.fof -d $3 -i $4 -b $2_vs_$1_$3_$4_rails -q $2-formatted.fof -p $5
echo RAILS process terminated.
