#!/usr/bin/perl

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
#   LINKS, RAILS and Cobbler Copyright (c) 2014-2017 Canada's Michael Smith Genome Science Centre.  All rights reserved.

use strict;
use Getopt::Std;
use Net::SMTP;
use vars qw($opt_f $opt_s $opt_d $opt_i $opt_v $opt_b $opt_t $opt_q);
getopts('f:s:d:v:b:t:i:q:');
my ($base_name,$frag_dist,$seqid,$verbose)=("",250,0.9,0);

my $version = "[v0.3]";
my $dev = "rwarren\@bcgsc.ca";
my $SAMPATH = "/gsc/btl/linuxbrew/bin/samtools";

#-------------------------------------------------

if(! $opt_f || ! $opt_s || ! $opt_q){
   print "Usage: $0 $version\n";
   print "-f  Assembled Sequences to further scaffold (Multi-FASTA format NO LINE BREAKS, required)\n"; 
   print "-q  Long Sequences queried (Multi-FASTA format NO LINE BREAKS, required)\n";
   print "-s  BAM file (use v0.2 for reading SAM files)\n";
   print "-d  Anchoring bases on contig edges (ie. minimum required alignment size on contigs, default -d $frag_dist, optional)\n";
   print "-i  Minimum sequence identity, default -i $seqid, optional\n";
   print "-t  LIST of names/header, long sequences to avoid using for merging/gap-filling scaffolds (optional)\n"; 
   print "-b  Base name for your output files (optional)\n";
   die   "-v  Runs in verbose mode (-v 1 = yes, default = no, optional)\n"; 
}

my $file = $opt_f;
my $longfile = $opt_s;
my $queryfile = $opt_q;
$frag_dist = $opt_d if($opt_d);
$seqid = $opt_i if($opt_i);
$verbose = $opt_v if($opt_v);
my $listfile = $opt_t if($opt_t);
$base_name = $opt_b if($opt_b);

my $assemblyruninfo="";


if(! -e $file){
   die "Invalid file: $file -- fatal\n";
}
if(! -e $longfile){
   die "Invalid file: $longfile -- fatal\n";
}


### Naming output files
if ($base_name eq ""){

   $base_name = $file . ".scaff_s-" . $longfile . "_q-" . $queryfile . "_d" . $frag_dist . "_i" . $seqid . "_t" . $listfile;

   my $pid_num = getpgrp(0);
   $base_name .= "_pid" . $pid_num;
}

my $log = $base_name . ".log";
my $newassemblyfile = $base_name . ".fa";
my $tsvfile = $base_name . "-list.tsv";

open (LOG, ">$log") || die "Can't write to $log -- fatal\n";


#-------------------------------------------------

my $init_message = "\nRunning: $0 $version\n-f $file\n-s $longfile\n";

$init_message .= "-d $frag_dist\n-i $seqid\n-t $listfile\n";

print $init_message;
print LOG $init_message;
$assemblyruninfo=$init_message . "\n";

#-------------------------------------------------

my $date = `date`;
chomp($date);

my $reading_reads_message = "\n=>Reading sam: $date\n";
print $reading_reads_message;
print LOG $reading_reads_message;
$assemblyruninfo.=$reading_reads_message;
my $tigpair;
my $initpos=0;
my $totalpairs=0;

my $delta = $frag_dist;

my $rh = &readSeqMemory($queryfile);
$tigpair = &readSam($longfile,$frag_dist,$seqid,$listfile,$delta,$initpos,$rh);
my $date = `date`;
chomp($date);
my $patchmsg = "done.\nFixing ambiguous bases (Ns): $date\n";
print $patchmsg;
print LOG $patchmsg;
$assemblyruninfo.=$patchmsg;
my ($gsl,$totalgap) = &patchGaps($file,$tigpair,$newassemblyfile,$tsvfile);

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

#exit;

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
   my ($file,$tigpair,$newfile,$gaplist) = @_;

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
   print TSV "scaffold\tscaftig\tgapLength\tgapFilledLength\n";   
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

	       if(defined $tigpair->{$num}{$next}{'seq'} && length($tigpair->{$num}{$next}{'seq'})>0){
                  print OUT "$tigpair->{$num}{$next}{'seq'}";
                  $filledct++;
               }else{
                  if($gaps[$gappos] ne ""){
                     print OUT "$gaps[$gappos]";
		  }
	       }
	       my $gaplen = length($gaps[$gappos]);
	       my $filledlen = length($tigpair->{$num}{$next}{'seq'});
               print TSV "$ctseq\t$scaftignum\t$gaplen\t$filledlen\n";
               push @gapspatched, $filledlen if($filledlen > 0);
	       #print "Scaftig $num -- $next\n $tigpair->{$num}{$next}{'distr'}\n$tigpair->{$num}{$next}{'seq'}\n$tigpair->{$num}{$next}{'configuration'}\n$tigpair->{$num}{$next}{'origin'}\nGAP:$gappos :: $gaps[$gappos]\n" if(defined $tigpair->{$num}{$next});
            }else{
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
sub readSam{

   my ($samfile,$frag_dist,$seqid,$listfile,$delta,$initpos,$rh) = @_;

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
   my $cutoff = $delta;
   my ($track_all,$tigpair);
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
   my $min=1;

   my $ERRLOG = $samfile.".bampreprocessor.err.log".$$.time();
   my $cmd = "$SAMPATH view $samfile 2>$ERRLOG|";
   open(IN,$cmd) || die "Error reading $samfile -- fatal.\n";
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


      if($si >= $seqid && $qalen >= $cutoff && (( $rstart <= $min &&  ($qlen-$qend)<= $min) || ($qstart<=$min && ($rlen-$rend)<=$min )    )){     ### this indicates anchoring bases, within $cutoff of edges


#         print "$si >= $seqid && $qalen >= $cutoff && (( $rstart <= $min &&  ($qlen-$qend)<= $min) || ($qstart<=$min && ($rlen-$rend)<=$min\n";
         my $dir;
         my $start;
         my $end;
         ###Coordinates on the scaffolds
         if($rstart <= $min &&  ($qlen-$qend)<= $min){
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
               if(! defined $mem->{$rd} && ! defined $bt->{$one}{$two} && ($track_all->{$current}{'qstart'} > $track_all->{$prev}{'qend'} || $track_all->{$prev}{'qstart'} > $track_all->{$current}{'qend'}) && $rh->{$rd} ne ""){### ADDED NO NULL SEQUENCES 30JAN2016
                  $bt->{$one}{$two}=1;### only one kollector merge per site 
                  my $pos=0;
                  $pos = $track_all->{$prev}{'qend'} if($track_all->{$current}{'qstart'} > $track_all->{$prev}{'qend'});
                  $pos = $track_all->{$current}{'qend'} if($track_all->{$prev}{'qstart'} > $track_all->{$current}{'qend'});            
                  my $gapseqlen = $totalreadlength - ($track_all->{$prev}{'qalen'} + $track_all->{$current}{'qalen'});
                  print ">$rd   $pos @ $gapseqlen \n$rh->{$rd}\n\n" if($verbose);
                  my $patch = substr($rh->{$rd},$pos,$gapseqlen-1);

               #$patch = &reverseComplement($patch) if($track_all->{$prev}{'orient'} eq "-" && $track_all->{$current}{'orient'} eq "-"  );
           

                  print "GAP:$patch\n" if($verbose);

                  ###JUST SOME TEST CODE
                  if(defined $tigpair->{$track_all->{$prev}{'scaftig'}}{$track_all->{$current}{'scaftig'}}{'seq'}){
                     $occ++;
                     print "$prev ($track_all->{$prev}{'scaftig'})...$current ($track_all->{$current}{'scaftig'})\n$tigpair->{$track_all->{$prev}{'tig'}}{$track_all->{$current}{'tig'}}{'seq'}\nNEW GAP:\n$patch\n";
                     if($patch ne $tigpair->{$track_all->{$prev}{'scaftig'}}{$track_all->{$current}{'scaftig'}}{'seq'}){print "NOT SAME\n\n";}else{print "SAME\n\n";$same++;}
                  }

                  if($track_all->{$prev}{'orient'} eq $track_all->{$current}{'orient'}){### ff or rr
                     $tigpair->{$one}{$two}{'distr'}++;
                     $patch = &reverseComplement($patch) if($track_all->{$prev}{'orient'} eq "r");
	             $tigpair->{$one}{$two}{'seq'}=lc($patch);
                     #$tigpair->{$one}{$two}{'configuration'}=$track_all->{$prev}{'orient'} . $track_all->{$current}{'orient'};
                     $tigpair->{$one}{$two}{'origin'}=$rd;
                  }
                  print ">>>> $track_all->{$prev}{'scaftig'} $track_all->{$current}{'scaftig'} $patch\n" if($verbose);

                  print "$track_all->{$prev}{'sam'}\n$track_all->{$current}{'sam'}\n   x====x ($totalreadlength)  ($track_all->{$prev}{'qstart'}-$track_all->{$prev}{'qend'}:$track_all->{$prev}{'qalen'} $track_all->{$prev}{'orient'}) AND ($track_all->{$current}{'qstart'}-$track_all->{$current}{'qend'}:$track_all->{$current}{'qalen'} $track_all->{$current}{'orient'})\n===x    x==== sc$track_all->{$prev}{'tig'} ($track_all->{$prev}{'start'}-$track_all->{$prev}{'end'}:$track_all->{$prev}{'qalen'} $track_all->{$prev}{'orient'}) AND sc$track_all->{$current}{'tig'} ($track_all->{$current}{'start'}-$track_all->{$current}{'end'}:$track_all->{$current}{'qalen'} $track_all->{$current}{'orient'}) \n\n" if($verbose);
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
