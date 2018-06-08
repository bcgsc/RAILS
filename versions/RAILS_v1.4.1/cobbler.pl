#!/usr/bin/env perl

#AUTHOR
#   Rene Warren
#   rwarren at bcgsc.ca


#NAME
#RAILS: Radial Assembly Improvement by Long Sequence Scaffolding 
#Scaffolding and gap-closure using alignment of long sequences

#SYNOPSIS

#DOCUMENTATION
#   readme.md distributed with this software
#   We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca
#   If you use RAILS, the RAILS code or ideas, please cite our work
#  

#LICENSE
#   LINKS, RAILS and Cobbler Copyright (c) 2014-2018 Canada's Michael Smith Genome Science Centre.  All rights reserved.

use strict;
use Getopt::Std;
use Net::SMTP;
use vars qw($opt_f $opt_s $opt_d $opt_i $opt_v $opt_b $opt_t $opt_q $opt_l $opt_g $opt_p);
getopts('f:s:d:v:b:t:i:q:g:l:p:');
my ($base_name,$anchor,$seqid,$verbose,$minreads,$grace)=("",1000,0.9,0,1,1);

my $version = "[v0.5.1]";
my $dev = "rwarren\@bcgsc.ca";
my $SAMPATH = "";

#-------------------------------------------------

if(! $opt_f || ! $opt_s || ! $opt_q || ! $opt_p){
   print "Usage: $0 $version\n";
   print "-f  Assembled Sequences to further scaffold (Multi-FASTA format NO LINE BREAKS, required)\n"; 
   print "-q  File of filenames containing long Sequences queried (Multi-FASTA format NO LINE BREAKS, required)\n";
   print "-s  File of filenames containing full path to BAM file(s) (use v0.2 for reading SAM files)\n";
   print "-p  Full path to samtools (known to work/tested with v1.8, required)\n";
   print "-d  Anchoring bases on contig edges (ie. minimum required alignment size on contigs, default -d $anchor, optional)\n";
   print "-i  Minimum sequence identity fraction (0 to 1), default -i $seqid, optional\n";
   print "-l  Minimum number of long sequence support per gap, default -l $minreads, optional\n";
   print "-g  Grace length (bp), default -g $grace, optional\n";
   print "-t  LIST of names/header, long sequences to avoid using for merging/gap-filling scaffolds (optional)\n"; 
   print "-b  Base name for your output files (optional)\n";
   print "-v  Runs in verbose mode (-v 1 = yes, default = no, optional)\n"; 
   die   "IMPORTANT: the order of files in -q and -s MUST match!\n";
}

my $file = $opt_f;
my $fof = $opt_s;
my $queryfof = $opt_q;
$anchor = $opt_d if($opt_d);
$seqid = $opt_i if($opt_i);
$verbose = $opt_v if($opt_v);
my $listfile = $opt_t if($opt_t);
$base_name = $opt_b if($opt_b);
$grace = $opt_g if($opt_g);
$minreads = $opt_l if($opt_l);
$SAMPATH = $opt_p if($opt_p);

my $assemblyruninfo="";


if(! -e $file){
   die "Invalid file: $file -- fatal\n";
}
if(! -e $fof){
   die "Invalid file: $fof -- fatal\n";
}
if(! -e $SAMPATH){
   die "Invalid : $SAMPATH -- fatal\n";
}

### Naming output files
if ($base_name eq ""){

   $base_name = $file . ".scaff_s-" . $fof . "_q-" . $queryfof . "_d" . $anchor . "_i" . $seqid . "_l" . $minreads . "_g" . $grace . "_t" . $listfile;
   my $pid_num = getpgrp(0);
   $base_name .= "_pid" . $pid_num;
}

my $log = $base_name . ".log";
my $newassemblyfile = $base_name . ".fa";
my $tsvfile = $base_name . "-list.tsv";

open (LOG, ">$log") || die "Can't write to $log -- fatal\n";


#-------------------------------------------------

my $init_message = "\nRunning: $0 $version\n-f $file\n-q $queryfof\n-s $fof\n";

$init_message .= "-d $anchor\n-i $seqid\n-l $minreads\n-g $grace\n-t $listfile\n";

print $init_message;
print LOG $init_message;
$assemblyruninfo=$init_message . "\n";

#-------------------------------------------------

my $date = `date`;
chomp($date);

my $reading_reads_message = "\n=>Reading bam: $date\n";
print $reading_reads_message;
print LOG $reading_reads_message;
$assemblyruninfo.=$reading_reads_message;
my $tigpair;
my $initpos=0;
my $totalpairs=0;

### READ Query read FOF
my @qryfilearray;
open(QRYFOF,$queryfof) || die "Can't open $queryfof for reading -- fatal.\n";
while(<QRYFOF>){
   chomp;
   push @qryfilearray, $_;
}
close QRYFOF;

my $ctline=0;
open(FOF,$fof) || die "Can't read $fof -- fatal.\n";
while(<FOF>){
   chomp;
   my $bamfile = $_;
   my $rh = &readSeqMemory($qryfilearray[$ctline]);
   print "Parsing alignment file $bamfile...\n";
   $tigpair = &readBam($tigpair,$bamfile,$anchor,$seqid,$listfile,$initpos,$rh,$grace);
   print "done.\n";
   $ctline++;
}
close FOF;

my $date = `date`;
chomp($date);
my $patchmsg = "done.\nFixing ambiguous bases (Ns): $date\n";
print $patchmsg;
print LOG $patchmsg;
$assemblyruninfo.=$patchmsg;
my ($gsl,$totalgap) = &patchGaps($file,$tigpair,$newassemblyfile,$tsvfile,$minreads);

my $date = `date`;
chomp($date);
my ($avg,$sum,$max,$min) = &average($gsl);
my $sd = &stdev($gsl);
my $final_message = "done: $date\n\n--------------- $0 Summary ---------------\nNumber of gaps patched : %i out of %i (%.2f %%)\nAverage length (bp) : %.2f\nLength st.dev +/- : %.2f\nTotal bases added : %i\nLargest gap resolved (bp) : %i\nShortest gap resolved (bp) : %i\n---------------------------------------------\n";
my @arrsg=@$gsl;
my $numgaps = $#arrsg+1;
my $percentclosed = $numgaps / $totalgap *100;
printf $final_message, ($numgaps,$totalgap,$percentclosed,$avg,$sd,$sum,$max,$min);
printf LOG $final_message, ($numgaps,$totalgap,$percentclosed,$avg,$sd,$sum,$max,$min);

$assemblyruninfo .= "done: $date\n\n--------------- $0 Summary ---------------\nNumber of gaps patched : $numgaps out of $totalgap ($percentclosed %) \nAverage length (bp) : $avg\nLength st.dev +/- : $sd\nTotal bases added : $sum\nLargest gap resolved (bp) : $max\nShortest gap resolved (bp) : $min\n---------------------------------------------\n";

exit;

###for dev. test purposes
eval{
   my $wdir = `pwd`;
   chomp($wdir);
   my $smtp = Net::SMTP->new('mailhost');
   $smtp->mail("RAILS\@bcgsc.ca");
   $smtp->to($dev);
   $smtp->data();
   $smtp->datasend("Subject: Your $0 run\n");
   $smtp->datasend("At: $wdir\n");
   $smtp->datasend($assemblyruninfo);
   $smtp->dataend();
   $smtp->quit;
};

exit;


#-----------------
sub readSeqMemory{

      my $file = shift;

      my $fh;
      my $prev="NA";
      my $seq="";
      open(FA,$file) || die "Cannot open $file for reading -- fatal.\n";
      while(<FA>){
	      chomp;
	      if (/\>(\S+)/){
		      my $head=$1;
		      if($prev ne $head && $prev ne "NA"){
			      $fh->{$prev} = $seq;
		      }
		      $prev = $head;
		      $seq='';
	      }elsif(/^(\S+)$/){
		      $seq .= uc($1);
	      }
      }
      $fh->{$prev} = $seq;

      close FA;

      return $fh;
}

#----------------
sub patchGaps{
   my ($file,$tigpair,$newfile,$gaplist,$minreads) = @_;

   my $tignames;
   my $head ="";
   my $ctseq=0;
   open(IN,$file) || die "Error reading $file -- fatal.\n";
   open(OUT,">$newfile") || die "Error writing $newfile -- fatal.\n";;
   open(TSV,">$gaplist") || die "Error reading $gaplist -- fatal.\n";
   my $filledct=0;
   my $totalgap=0;
   my @gapspatched;
   print "\nSequences processed:\n";
   print TSV "scaffold\tscaftig\tgapLength\tgapFilledLength\treadSupportCount\n";   
   while(<IN>){
      chomp;	   
      if(/^\>(\S+)/){
         print OUT "$_\n";
         $head=$1;
         $ctseq++;
         $tignames->{$ctseq}=$head;
         print "\r$ctseq";
         $|++;
      }else{
         my @scaftigs = split(/N+/i,$_); 
	 my @gaps = split(/[ABCDEFGHIJKLMOPQRSTUVWXYZ]+/i,$_);###Anything but Ns
         my $numgap = $#gaps;
	 $totalgap += $numgap if($numgap>0) ;
         #print "@gaps $numgap\n";#XXXXX

         my $scaftignum=0;
         my $gappos=1;######ASSUMES SCAFFOLDS NEVER START WITH Ns
	 foreach my $scaftig(@scaftigs){ 
            print OUT "$scaftig";
            my $len = length($scaftig);
            $scaftignum++;
	    my $num = $ctseq . "." . $scaftignum;
            my $nextscaftig = $scaftignum+1;
            #print "$num\n";
	    if(defined $tigpair->{$num}){
               my $list = $tigpair->{$num};
               my $next = $ctseq . "." . $nextscaftig;

	       if(defined $tigpair->{$num}{$next}{'seq'} && length($tigpair->{$num}{$next}{'seq'})>0 && $tigpair->{$num}{$next}{'distr'}>=$minreads){###MIN READ LOGIC
                  print OUT "$tigpair->{$num}{$next}{'seq'}"; ### gap-filling
                  $filledct++;
	          my $gaplen = length($gaps[$gappos]);
	          my $filledlen = length($tigpair->{$num}{$next}{'seq'});
                  print TSV "$ctseq\t$scaftignum\t$gaplen\t$filledlen\t$tigpair->{$num}{$next}{'distr'}\n";
                  push @gapspatched, $filledlen if($filledlen > 0);
	          #print "Scaftig $num -- $next\n $tigpair->{$num}{$next}{'distr'}\n$tigpair->{$num}{$next}{'seq'}\n$tigpair->{$num}{$next}{'configuration'}\n$tigpair->{$num}{$next}{'origin'}\nGAP:$gappos :: $gaps[$gappos]\n" if(defined $tigpair->{$num}{$next});
               }else{### Does not pass filters, put back the Ns
                  if($gaps[$gappos] ne ""){
                     print OUT "$gaps[$gappos]";
                     my $gaplen = length($gaps[$gappos]);
                     print TSV "$ctseq\t$scaftignum\t$gaplen\t\t$tigpair->{$num}{$next}{'distr'}\n";### will still indicate read support
                  }
               }
            }else{### Does not pass filters, put back the Ns
               print OUT "$gaps[$gappos]";
	       my $gaplen = length($gaps[$gappos]);
	       print TSV "$ctseq\t$scaftignum\t$gaplen\t\n";
	    }
	    $gappos++;
         }	    
         print OUT "\n";
      }
   }
   close IN;
   close OUT;
   close TSV;
   print "\ndone.\n";
   my $endmessage = "Filled $filledct out of $totalgap gaps (gaps are defined by any stretch of Ns in your assembly)\nGap-filled assembly: $newfile\nList of gap lengths: $gaplist\n";
   print LOG "$endmessage";
   print "$endmessage";
   $assemblyruninfo .= $endmessage;

   return \@gapspatched,$totalgap;
}

#---------------
sub readBam{

   my ($tigpair,$bamfile,$anchor,$seqid,$listfile,$initpos,$rh,$grace) = @_;

   my $mem;
   if(-f $listfile){ 
      open(IN,$listfile) || die "Can't read $listfile -- fatal.\n";
      while(<IN>){
         chomp;
         $mem->{$_}=1;
      }
      close IN;
   }
   my $bt;
   my $track_all;
#HS9_159:6:1308:13492:64472      272     scaffold43,6983,f43Z6983        6439    0       536M    *       0       0       *       *       NM:i:0  AS:i:536
#HS9_159:6:1308:13492:64472      0       scaffold30,32025,f30Z32025      25411   0       536M    *       0       0       GCTTATAAAAGAAGGTGCAATTGATCCTTGCCTTACGCCTACAAAGGAGGGTAGGTGCGATTGGTCCTTACATTCTTACGCCGCTTAGGAAGCTAGGCGAGATAGGATGGGTTCTAGAGCACCTAACTAGCTTTACACGCCGAATCCAGACCTGCCGGCTACCATCCGGATTCATACTAGATAACATAAAGGAGAGAACAACTGTTCAAAGAACAACTCGGAGAACATTTGTATCCGGTGGTTGGGGCATTGCGTGCTATACCAACTACCTCAGGTGCGCGAGGTCTCATTCCTTTTCCAAGCCCAATAAAGAAAAAATATCATTAGTGATGGTGAATCCCGTTTATATAAGTAAGTTGCATTCTTATCTAAGTAAGTGGGCTTTCCTAAGTCACTTATTGGGTGGGGGGCCCCTGTCGAGTGAGCCATCCTTCCTCACCCTCTCTTTTGTTGGGCGAGCCATCTTTCCTTTTATACGATTCGATCCAGTAGATAAGGAAGACCGACCGAGAACAACCAATGGCCTTCCCTGGGGG        *       NM:i:0  AS:i:536        XS:i:536
#HS9_159:6:1308:13492:64472      272     scaffold22,90777,f22Z90777      90233   0       536M    *       0       0       *       *       NM:i:0  AS:i:536
   my $t;
   my $ct=0;

   my %options = ();

   print join(
    "\t",
    'qname',
    'qstart',
    'qend',
    'qalen',
        'qlen',
    'rname',
    'rstart',
    'rend',
    'ralen',
        'rlen',
        'edit_dist',
   ) . "\n" if $options{header};


   my %rlength = ();

   my $ERRLOG = $bamfile.".bampreprocessor.err.log".$$.time();
   my $cmd = "$SAMPATH view $bamfile 2>$ERRLOG|";
   open(IN,$cmd) || die "Error reading $bamfile -- fatal.\n";
   while(<IN>){

      chomp;
      $ct++;

      my @a=split(/\t/);
      my @b=split(/\,/,$a[2]);
      my @c=split(/\,/,$a[0]);

      if ($options{rlen} && /^\@SQ\s+SN:(\S+)\s+LN:(\S+)/) {
         $rlength{$1} = $2;
      }
      next unless @a >= 10;
      my $line = $_;
      my $qname = $a[0];
      my $rname = $a[2];
      my $rstart = $a[3];
      my $cigar = $a[5];
      my $qseq = $a[9];
      # Query
      my $qstart = 1;
      $_ = $cigar;
      s/^(\d+)[SH]/$qstart += $1/eg;
      my $qalen = 0;
      $_ = $cigar;
      s/(\d+)[M=XI]/$qalen += $1/eg;
      my $qend = $qstart + $qalen - 1;
      $_ = $cigar;
      my $end_clip_len = 0;
      s/(\d+)[SH]$/$end_clip_len += $1/eg;
      my $qlen = $c[1];
      #if ($qalen > 0) {
      #   $qlen = ($qstart-1) + $qalen + $end_clip_len;
      #} elsif ($qseq ne "*") {
      #   $qlen = length($a[9]);
      #}

      # Reference
      my $ralen = 0;
      $_ = $cigar;
      s/(\d+)[M=XDN]/$ralen += $1/eg;
      my $rend = $rstart + $ralen - 1;
      my $rlen = $b[1];
      #if ($options{rlen} && exists($rlength{$rname})) {
      #   $rlen = $rlength{$rname};
      #}

      # Calculate edit distance including clipping
      my $edit_dist = '';
      if ($line =~ /NM:i:(\d+)/) {
         $edit_dist = $1;# + $qstart - 1 + $end_clip_len;
      }

#      if ($rname eq '*') {
#        # case: query sequence is unmapped
#        print join("\t", $qname, $qstart, $qend, $qalen, $qlen) . "\n";
#      } else {
#        print join("\t", $qname, $qstart, $qend, $qalen, $qlen, $rname, $rstart, $rend, $ralen, $rlen);
#        print "\t$edit_dist" if length($edit_dist) > 0;
#        print "\n";
#      }

      my $read = $a[0] . "-" . $ct;
      my $si=0;
      $si = ($qalen - $edit_dist) / $qalen if($qalen);


      if($si >= $seqid && $qalen >= $anchor && (( $rstart <= $grace &&  ($qlen-$qend)<= $grace) || ($qstart<=$grace && ($rlen-$rend)<=$grace )    )){     ### this indicates anchoring bases, within $anchor of edges


         print "$si >= $seqid && $qalen >= $anchor && (( $rstart <= $grace &&  ($qlen-$qend)<= $grace) || ($qstart<=$grace && ($rlen-$rend)<=$grace\n" if($verbose);
         my $dir;
         my $start;
         my $end;
         ###Coordinates on the scaffolds
         if($rstart <= $grace &&  ($qlen-$qend)<= $grace){
            $start = $rend;
            $end = $rstart;
         }else{
            $start = $rstart;
            $end = $rend;
         }
         my $orient="";
         if($a[1]==272 || $a[1]==16 || $a[1]==2064){ ### matches on negative strand
            $orient="r";
            my $tmpstart = $qlen - $qend;
            my $tmpend = $qlen - $qstart;
            $qstart = $tmpstart;
            $qend = $tmpend;
         }else{
            $orient="f";
         }
         ###tracks from a read perspective
	 my ($numtig,$scaftignum,$sz)=($2,$1,$3) if($a[2]=~/\D+((\d+)\.\d+),(\d+)/);### scaffoldNUMBER,LENGTH eg. wga1,1301
         $t->{$a[0]}{$scaftignum}{'orient'}= $dir . $orient ;
         $t->{$a[0]}{$scaftignum}{'real'}=$read;###my $read = $a[0] . "-" . $ct;
         $t->{$a[0]}{$scaftignum}{'length'}=$qlen;
         $track_all->{$read}{'tig'}=$numtig;
         $track_all->{$read}{'scaftig'}=$scaftignum;
         $track_all->{$read}{'start'}=$start;
         $track_all->{$read}{'end'}=$end;
         $track_all->{$read}{'multiple'}=1;
         $track_all->{$read}{'sam'}=$line;
         $track_all->{$read}{'orient'}=$orient;
         $track_all->{$read}{'qalen'}=$qalen;
         $track_all->{$read}{'qstart'}=$qstart;
         $track_all->{$read}{'qend'}=$qend;
         $track_all->{$read}{'si'}=$si; ### added 11APR2018  the read with most matching bases is chosen for gapfill (patch seq)
#         print "$line\n\n";
      }
   }
   close IN;###End SAM parse
   my ($occ,$same)=(0,0);###TRACK REDUNDANCY

   foreach my $rd(keys %$t){
      my $scafflist=$t->{$rd};
      my $num = keys(%$scafflist);
      my $prevscaff = "NA";
      foreach my $scaff(sort {$a<=>$b} keys %$scafflist){
         if($prevscaff ne "NA"){
      #if($num==2){###maps on two different scaftigs only
         #print "$num!\n";
            my @arr;
            my $totalreadlength=0;
            my $current = $scafflist->{$scaff}{'real'};
            my $prev = $scafflist->{$prevscaff}{'real'};
            $totalreadlength = $scafflist->{$scaff}{'length'};

            my ($p_s,$p_t)=($1,$2) if($track_all->{$prev}{'scaftig'}=~/(\d+)\.(\d+)/);
            my ($c_s,$c_t)=($1,$2) if($track_all->{$current}{'scaftig'}=~/(\d+)\.(\d+)/);
            my $prev_match = $p_s . "." . ($p_t + 1);
            my $curr_match = $c_s . "." . ($c_t + 1);

	 #print "$track_all->{$current}{'tig'} == $track_all->{$prev}{'tig'} $track_all->{$current}{'scaftig'} $track_all->{$prev}{'scaftig'}\n";
         #print "$track_all->{$current}{'scaftig'} ... $track_all->{$prev}{'scaftig'}\n";
            if($track_all->{$current}{'tig'} == $track_all->{$prev}{'tig'} && (($track_all->{$current}{'scaftig'} eq $prev_match)  || ($track_all->{$prev}{'scaftig'} eq $curr_match ))){###ADDED OCT2016, make sure on same scaffold and consecutive

               my ($one,$two)=($track_all->{$current}{'scaftig'},$track_all->{$prev}{'scaftig'});
               if($p_t < $c_t){
                 ($one,$two)=($track_all->{$prev}{'scaftig'},$track_all->{$current}{'scaftig'})
               }
    	       #print "$one,$two\n";
               ### this will track the best anchoring long reads for the merge/gapfill
               my $m1 = $track_all->{$current}{'qalen'} * $track_all->{$current}{'si'};
               my $m2 = $track_all->{$prev}{'qalen'} * $track_all->{$prev}{'si'};
               my $matchbases = $m1 + $m2;

               if(! defined $mem->{$rd} && ($track_all->{$current}{'qstart'} > $track_all->{$prev}{'qend'} || $track_all->{$prev}{'qstart'} > $track_all->{$current}{'qend'}) && $rh->{$rd} ne ""){### WILL TRACK BEST ANCHORING BASES
                $tigpair->{$one}{$two}{'distr'}++;
                if($matchbases > $bt->{$one}{$two}{'bestmatch'}){
                  $bt->{$one}{$two}{'bestmatch'} = $matchbases;
                  my $pos=0;
                  $pos = $track_all->{$prev}{'qend'} if($track_all->{$current}{'qstart'} > $track_all->{$prev}{'qend'});
                  $pos = $track_all->{$current}{'qend'} if($track_all->{$prev}{'qstart'} > $track_all->{$current}{'qend'});            
                  my $gapseqlen = $totalreadlength - ($track_all->{$prev}{'qalen'} + $track_all->{$current}{'qalen'});
                  print ">$rd   $pos @ $gapseqlen \n$rh->{$rd}\n\n" if($verbose);
                  my $patch = substr($rh->{$rd},$pos,$gapseqlen-1);

                  print "GAP:$patch\n" if($verbose);

                  ###JUST SOME TEST CODE
                  if(defined $tigpair->{$track_all->{$prev}{'scaftig'}}{$track_all->{$current}{'scaftig'}}{'seq'}){
                     $occ++;
                     $same if($patch eq $tigpair->{$track_all->{$prev}{'scaftig'}}{$track_all->{$current}{'scaftig'}}{'seq'});

                     print "$prev ($track_all->{$prev}{'scaftig'})...$current ($track_all->{$current}{'scaftig'})\n$tigpair->{$track_all->{$prev}{'tig'}}{$track_all->{$current}{'tig'}}{'seq'}\nNEW GAP:\n$patch\n" if($verbose);

                     #if($patch ne $tigpair->{$track_all->{$prev}{'scaftig'}}{$track_all->{$current}{'scaftig'}}{'seq'}){print "NOT SAME\n\n";}else{print "SAME\n\n";}
                  }

                  if($track_all->{$prev}{'orient'} eq $track_all->{$current}{'orient'}){### ff or rr
                     $patch = &reverseComplement($patch) if($track_all->{$prev}{'orient'} eq "r");
	             $tigpair->{$one}{$two}{'seq'}=lc($patch);
                     $tigpair->{$one}{$two}{'origin'}=$rd;
                  }
                  print ">>>> $track_all->{$prev}{'scaftig'} $track_all->{$current}{'scaftig'} $patch\n" if($verbose);

                  print "$track_all->{$prev}{'sam'}\n$track_all->{$current}{'sam'}\n   x====x ($totalreadlength)  ($track_all->{$prev}{'qstart'}-$track_all->{$prev}{'qend'}:$track_all->{$prev}{'qalen'} $track_all->{$prev}{'orient'}) AND ($track_all->{$current}{'qstart'}-$track_all->{$current}{'qend'}:$track_all->{$current}{'qalen'} $track_all->{$current}{'orient'})\n===x    x==== sc$track_all->{$prev}{'tig'} ($track_all->{$prev}{'start'}-$track_all->{$prev}{'end'}:$track_all->{$prev}{'qalen'} $track_all->{$prev}{'orient'}) AND sc$track_all->{$current}{'tig'} ($track_all->{$current}{'start'}-$track_all->{$current}{'end'}:$track_all->{$current}{'qalen'} $track_all->{$current}{'orient'}) \n\n" if($verbose);
                }###save for bestmatch only
               }
	    }
	 }#IF PREV NE NA
	 $prevscaff = $scaff;
      }#foreach scatigs, ordered
   }

   print "\nRedundant same contig combo linking:$occ\nSame gap sequence fill:$same\n\n";

   return $tigpair;
}

#-----------------------
sub reverseComplement{
   $_ = shift;
   $_ = uc();
   tr/ATGCYRKMBDHV/TACGRYMKVHDB/;
   return (reverse());
}

#----------------
sub average{
        my $data = shift;
        if (not @$data) {
                die("Empty arrayn -- maybe the scaffold merging step did not necessitate gap filling.");
        }
        my $total = 0;
        my $max = 0;
        my $min = 1000000;
        foreach (@$data) {
                $total += $_;
                $max = $_ if($_ > $max);
                $min = $_ if($_ < $min);
        }
        my $average = $total / @$data;
        return $average,$total,$max,$min;
}

#----------------
sub stdev{
        my $data = shift;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

## We hope this code is useful to you -- Please send comments & suggestions to rwarren at bcgsc.ca
