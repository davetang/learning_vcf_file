#!/usr/bin/env perl

# Filter based on the fields in the INFO column

use strict;
use warnings;

my $usage = "Usage: $0 <infile.vcf> <field> <threshold>\n";
my $infile = shift or die $usage;
my $field = shift or die $usage;
my $threshold = shift or die $usage;

open(IN,'<',$infile) || die "Could not open $infile: $!\n";
while(<IN>){
   chomp;
   if (/^#/){
      print "$_\n";
      next;
   }
   
   my $skip = 0;
   my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$lib) = split(/\t/);
   my @info = split(/;/, $info);
   foreach my $i (@info){
      if ($i =~ /=/){
         my ($key, $value) = split(/=/, $i);
         next unless $key eq $field;
         if ($value < $threshold){
            $skip = 1;
         }
      }
   }
   if ($skip == 1){
      next;
   } else {
      print "$_\n";
   }
}
close(IN);

exit(0);
