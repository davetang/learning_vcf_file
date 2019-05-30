#!/usr/bin/env perl

use strict;
use warnings;

open(IN, "<", "raw/sample.vcf") || die "Could not open raw/sample.vcf: $!\n";
while(<IN>){
   chomp;
   next if /^#/;
   my @line_split = split(/\t/);
   my $tally = 0;
   #print "$line_split[9]\n";
   foreach my $genotype (@line_split[10..18]){
      if ($genotype eq $line_split[9]){
         ++$tally;
      }
   }
   print "$_\n" if $tally == 0;
}

close(IN);


