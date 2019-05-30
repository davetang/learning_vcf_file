#!/usr/bin/env perl

# Add additional samples to an existing VCF file; this script only adds genotypes!

use strict;
use warnings;

my $usage = "Usage: $0 <infile.vcf> <number of sample> <seed>\n";
my $infile = shift or die $usage;
my $num = shift or die $usage;
my $seed = shift or die $usage;

srand($seed);

if ($infile =~ /\.gz$/){
   open(IN, '-|', "gunzip -c $infile") || die "Could not open $infile: $!\n";
} else {
   open(IN, '<', $infile) || die "Could not open $infile: $!\n";
}

my @sample = 1 .. $num;

while(<IN>){
   chomp;
   if (/^#/){
      if (/CHROM/){
         print join("\t", "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", @sample), "\n";
      } else {
         print "$_\n";
      }
      next;
   }
   #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
   my @line_split = split(/\t/);
   my @gt = ();
   foreach my $sample (@sample){
      push(@gt, int(rand(2)) == 1 ? "1/1" : "0/0");
   }
   print join("\t", @line_split[0..8], @gt), "\n";
}
close(IN);

exit(0);

