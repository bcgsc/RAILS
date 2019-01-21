#!/bin/bash
#RLW 2016
if [ $# -ne 7 ]; then
        echo "Usage: $(basename $0) <FASTA assembly .fa> <FASTA long sequences .fa> <anchoring sequence length eg. 250> <min sequence identity 0.95> <number of threads> <output dir> <path to samtools>"
        exit 1
fi

if [ -f $6 ]; then
	echo "ERROR File $3 exists"
	exit 1
fi

if [ ! -d $6 ]; then
	mkdir $6
fi

ASSEMBLY=`realpath $1`
LONGREADS=`realpath $2`
OUT=`realpath $6`
OUT=`echo $OUT/`

echo Resolving ambiguous bases -Ns- in $ASSEMBLY assembly using long sequences $LONGREADS
echo reformatting file $ASSEMBLY
perl -ne 'if(/^\>/){print $n . $_;}else{$_=~s/\n//;print $_;$n="\n"}END{print "\n";}' $ASSEMBLY > $OUT$(basename $ASSEMBLY)-wolb
cat $OUT$(basename $ASSEMBLY)-wolb | perl -ne 'if(/^\>/){$scafnum++;}else{my $len=length($_);my @scaftigs=split(/N+/i,$_);my $scaftignum=0;foreach my $scaftig(@scaftigs){ my $len=length($scaftig);$scaftignum++;  print ">wga$scafnum";print "."; print "$scaftignum,$len\n$scaftig\n";}}' > $OUT$(basename $ASSEMBLY)-formatted.fa
echo reformatting file $LONGREADS
cat $LONGREADS | perl -ne 'if(/^\>/){$ct++;}else{my $len=length($_);print ">seq$ct,$len\n$_";}' > $OUT$(basename $LONGREADS)-formatted.fa
echo Building sequence database index out of your $OUT$(basename $ASSEMBLY)-formatted.fa assembly contigs..
bwa index $OUT$(basename $ASSEMBLY)-formatted.fa
echo Aligning long sequences $OUT$(basename $LONGREADS)-formatted.fa to your contigs..
### YOU MAY CONSIDER: SETTING THE MORE STRINGENT bwa mem -x intractg OPTION
bwa mem -a -t$5 $OUT$(basename $ASSEMBLY)-formatted.fa $OUT$(basename $LONGREADS)-formatted.fa | samtools view -Sb - > $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_gapfilling.bam
echo Scaffolding $OUT$(basename $ASSEMBLY)-formatted.fa using $OUT$(basename $LONGREADS)-formatted.fa and filling gaps with sequences in $OUT$(basename $LONGREADS)-formatted.fa
echo $OUT$(basename $LONGREADS)-formatted.fa > $OUT$(basename $LONGREADS)-formatted.fof
echo $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_gapfilling.bam > $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_gapfilling.fof
cobbler.pl -f $OUT$(basename $ASSEMBLY)-wolb -s $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_gapfilling.fof -d $3 -i $4 -b $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_$3_$4_gapsFill -q $OUT$(basename $LONGREADS)-formatted.fof -p $7
echo Process terminated.
echo RAILS scaffolding $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_$3_$4_gapsFill.fa sequences using long seqs $OUT$(basename $LONGREADS)-formatted.fa -- anchoring sequence threshold $3 bp 
echo reformatting file $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_$3_$4_gapsFill.fa
cat $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_$3_$4_gapsFill.fa | perl -ne 'if(/^\>/){$ct++;}else{my $len=length($_);print ">wga$ct,$len\n$_";}' > $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_$3_$4_gapsFill-formatted.fa
echo Building sequence database index out of your $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_$3_$4_gapsFill-formatted.fa assembly contigs..
bwa index $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_$3_$4_gapsFill-formatted.fa
echo Aligning long sequences $OUT$(basename $LONGREADS)-formatted.fa to your contigs..
### YOU MAY CONSIDER: SETTING THE MORE STRINGENT bwa mem -x intractg OPTION
bwa mem -a -t$5 $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_$3_$4_gapsFill-formatted.fa $OUT$(basename $LONGREADS)-formatted.fa | $7 view -Sb - > $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_scaffolding.bam
echo Scaffolding $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_$3_$4_gapsFill-formatted.fa using $OUT$(basename $LONGREADS)-formatted.fa and filling new gaps with sequences in $OUT$(basename $LONGREADS)-formatted.fa
echo $OUT$(basename $LONGREADS)-formatted.fa > $OUT$(basename $LONGREADS)-formatted.fof
echo $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_scaffolding.bam > $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_scaffolding.fof
RAILS -f $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_$3_$4_gapsFill-formatted.fa -s $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_scaffolding.fof -d $3 -i $4 -b $OUT$(basename $LONGREADS)_vs_$(basename $ASSEMBLY)_$3_$4_rails -q $OUT$(basename $LONGREADS)-formatted.fof -p $7
echo RAILS process terminated.
