#!/usr/bin/env perl

# Convert a VCF file into a (UCSC Genome Browser) BED file; not a binary PED file

use strict;
use warnings;

my $usage = "Usage: $0 <infile.vcf>\n";
my $infile = shift or die $usage;

if ($infile =~ /\.gz$/){
   open(IN, '-|', "gunzip -c $infile") || die "Could not open $infile: $!\n";
} else {
   open(IN, '<', $infile) || die "Could not open $infile: $!\n";
}

while(<IN>){
   chomp;
   next if (/^#/);
   #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
   my ($chr, $pos, $id, $ref, $alt, $qual, $filter, @rest) = split(/\t/);
   my $start = $pos - 1;
   my $end = $pos;
   print join("\t", $chr, $start, $end, $id, $qual, '+'), "\n";
}
close(IN);

exit(0);

