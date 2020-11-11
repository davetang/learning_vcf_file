#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "Usage: $0 <bp> <seed>\n";
my $bp = shift or die $usage;
my $seed = shift or die $usage;

# set seed for reproducibility
srand($seed);

my $random_seq = random_seq($bp);
if ($bp > 60){
   $random_seq = add_newline($random_seq, 60);
   print ">$bp | $seed\n$random_seq";
} else {
   print ">$bp | $seed\n$random_seq\n";
}

exit(0);

sub random_seq {
   my ($bp) = @_;
   my @nuc = qw/ A C G T /;
   my $seq = '';
   for (1 .. $bp){
      my $rand_ind = int(rand(scalar(@nuc)));
      $seq .= $nuc[$rand_ind];
   }
   return($seq);
}

sub add_newline {
   my ($line, $pos) = @_;
   my $mod = '';
   for (my $i = 0; $i < length($line); $i += $pos){
      $mod .= substr($line, $i, $pos) . "\n";
   }
   return($mod);
}

exit(0);

