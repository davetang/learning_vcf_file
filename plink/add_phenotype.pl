#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "Usage: <infile.ped> <infile.panel>\n";
my $ped = shift or die $usage;
my $panel = shift or die $usage;

my %phenotype = ();
open(IN, '<', $panel) || die "Could not open $panel: $!\n";
while(<IN>){
   chomp;
   next if /^sample/;
   my ($sample, $pop, $super_pop, $gender, $phenotype) = split(/\t/);
   $phenotype{$sample} = $phenotype;
}
close(IN);

open(IN, '<', $ped) || die "Could not open $ped: $!\n";
while(<IN>){
   chomp;
   my @line = split(/\t/);
   my $indi = $line[1];
   if (exists $phenotype{$indi}){
      $line[5] = $phenotype{$indi};
   } else {
      die "Missing $indi in $panel\n";
   }
   print join("\t", @line), "\n";
}
close(IN);


exit(0);

__END__

