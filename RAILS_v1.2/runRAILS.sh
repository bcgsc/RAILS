#!/bin/bash
#RLW 2016
if [ $# -ne 4 ]; then
        echo "Usage: $(basename $0) <FASTA assembly .fa> <FASTA long sequences .fa> <anchoring sequence length eg. 250> <min sequence identity 0.95>"
        exit 1
fi
###Change line below to point to path of bwa executables
export PATH=/gsc/btl/linuxbrew/bin:$PATH
echo Resolving ambiguous bases -Ns- in $1 assembly using long sequences $2
echo reformatting file $1
cat $1 | perl -ne 'if(/^\>/){$scafnum++;}else{my $len=length($_);my @scaftigs=split(/N+/i,$_);my $scaftignum=0;foreach my $scaftig(@scaftigs){ my $len=length($scaftig);$scaftignum++;  print ">wga$scafnum";print "."; print "$scaftignum,$len\n$scaftig\n";}}' > $1-formatted.fa
echo reformatting file $2
cat $2 | perl -ne 'if(/^\>/){$ct++;}else{my $len=length($_);print ">seq$ct,$len\n$_";}' > $2-formatted.fa
cat $2 | perl -ne 'if(/^\>/){  }else{my $len=length($_);my @scaftigs=split(/N+/i,$_);foreach my $scaftig(@scaftigs){ $ct++;my $len=length($scaftig);  print ">seq$ct,$len\n$scaftig\n";}}' > $2-scaftigformatted.fa
echo Building sequence database index out of your $1-formatted.fa assembly contigs..
bwa index $1-formatted.fa
echo Aligning long sequences $2-scaftigformatted.fa to your contigs..
#bwa mem -x intractg -a -t24 $1-formatted.fa $2-formatted.fa  | /gsc/btl/linuxbrew/bin/samtools view -Sb - > $2_vs_$1_gapfilling.bam
bwa mem -a -t24 $1-formatted.fa $2-scaftigformatted.fa  | /gsc/btl/linuxbrew/bin/samtools view -Sb - > $2_vs_$1_gapfilling.bam
echo Scaffolding $1-formatted.fa using $2-scaftigformatted.fa and filling gaps with sequences in $2-formatted.fa
./cobbler.pl -f $1 -s $2_vs_$1_gapfilling.bam -d $3 -i $4 -b $2_vs_$1_$3_gapsFill -q $2-scaftigformatted.fa
echo Process terminated.
echo RAILS scaffolding $1.gapsFill.fa sequences using long seqs $2 -- anchoring sequence threshold $3 bp 
echo reformatting file $1.gapsFill.fa
cat $2_vs_$1_$3_gapsFill.fa | perl -ne 'if(/^\>/){$ct++;}else{my $len=length($_);print ">wga$ct,$len\n$_";}' > $2_vs_$1_$3_gapsFill-formatted.fa
echo Building sequence database index out of your $2_vs_$1_$3_gapsFill-formatted.fa assembly contigs..
bwa index $2_vs_$1_$3_gapsFill-formatted.fa
echo Aligning long sequences $2-formatted.fa to your contigs..
bwa mem -a -t24 $2_vs_$1_$3_gapsFill-formatted.fa $2-formatted.fa | /gsc/btl/linuxbrew/bin/samtools view -Sb - > $2_vs_$1_scaffolding.bam
echo Scaffolding $2_vs_$1_$3_gapsFill-formatted.fa using $2-formatted.fa and filling new gaps with sequences in $2-formatted.fa
./RAILS -f $2_vs_$1_$3_gapsFill-formatted.fa -s $2_vs_$1_scaffolding.bam -d $3 -i $4 -b $2_vs_$1_$3_rails -q $2-formatted.fa
echo RAILS process terminated.
