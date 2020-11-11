#!/usr/bin/env perl

# Simple script that takes an input fasta sequence
# and generates paired end reads

use strict;
use warnings;

my $usage = "Usage: $0 <infile.fa> <read length> <number of pairs> <inner mate distance> <seed>\n";
my $fasta = shift or die $usage;
my $len = shift or die $usage;
my $num = shift or die $usage;
my $inner_mate = shift or die $usage;
my $seed = shift or die $usage;

srand($seed);

my $seq = '';
open(IN,'<',$fasta) || die "Could not open $fasta: $!\n";
while(<IN>){
   chomp;
   next if /^>/;
   $seq .= $_;
}
close(IN);

my $seq_len = length($seq);
my $limit = $seq_len - $len -$len - $inner_mate + 1;

if ($len > $seq_len){
   die "Your read length is longer than the input sequence\n";
}

# on Illumina 1.5+ B is the worst quality
# and h is the best
my $fake_qual = 'J' x $len;

my $name = 'l' . $len . '_' . 'n' . $num . '_' . 'd' . $inner_mate . '_' . $seed;
my $first_out = $name . '_1.fq.gz';
my $second_out = $name . '_2.fq.gz';

open(my $read1, '|-', "gzip >$first_out") || die "Could not write output to $first_out: $!\n";
open(my $read2, '|-', "gzip >$second_out") || die "Could not write output to  $second_out: $!\n";

for (1 .. $num){

   my $first_start = int(rand($limit));
   if ($first_start > $limit){
      while( $first_start > $limit ){
         $first_start = int(rand($limit));
      }
   }

   my $first_read = substr($seq,$first_start,$len);
   my $first_pos = $first_start + 1;
   print $read1 "\@$_:$first_pos\n$first_read\n+\n$fake_qual\n";

   my $second_start = $first_start + $inner_mate;
   my $second_read = substr($seq,$second_start,$len);
   $second_read = reverse($second_read);
   $second_read =~ tr/ACGT/TGCA/;
   my $second_pos = $second_start + 1;
   print $read2 "\@$_:$first_pos\n$second_read\n+\n$fake_qual\n";
}

close($read1);
close($read2);

exit(0);

