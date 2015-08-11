#!/usr/bin/env perl

# Filter a VCF file based on the CHROM, QUAL, and FILTER columns
# or based on the fields in the INFO column

use strict;
use warnings;

my $usage = "Usage: $0 <infile.vcf> <column | field> <threshold>\n";
my $infile = shift or die $usage;
my $field = shift or die $usage;
my $threshold = shift or die $usage;

if ($infile =~ /\.gz$/){
   open(IN,'-|',"gunzip -c $infile") || die "Could not open $infile: $!\n";
} else {
   open(IN,'<',$infile) || die "Could not open $infile: $!\n";
}

while(<IN>){
   chomp;
   if (/^#/){
      print "$_\n";
      next;
   }

   # by default skip
   my $skip = 1;
   my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$lib) = split(/\t/);

   # if one of the columns is used as a "field" then
   # simply check the value in the corresponding column
   #CHROM  POS     ID      REF     ALT     QUAL    FILTER
   if ($field eq 'CHROM'){
      if ($threshold eq $chrom){
         print "$_\n";
         next;
      } else {
         next;
      } 
   } elsif ($field eq 'QUAL'){
      if ($threshold > $qual){
         print "$_\n";
         next;
      } else {
         next;
      } 
   } elsif ($field eq 'FILTER'){
      if ($threshold eq $filter){
         print "$_\n";
         next;
      } else {
         next;
      } 
   }

   # check the fields in the INFO column iteratively
   # and if condition is fulfilled, don't skip
   my @info = split(/;/, $info);
   foreach my $i (@info){
      if ($i =~ /=/){
         my ($key, $value) = split(/=/, $i);
         next unless $key eq $field;
         if ($value > $threshold){
            $skip = 0;
         }
      }
   }

   # the default is to skip, in case the user
   # inputs a field that cannot be found
   if ($skip == 1){
      next;
   } else {
      print "$_\n";
   }
}
close(IN);

exit(0);
