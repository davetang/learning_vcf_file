#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "Usage: $0 <bp> <seed>\n";
my $num = shift or die $usage;
my $seed = shift or die $usage;

#set seed for reproducibility
srand($seed);

my $random_seq = random_seq($num);
print ">$num | $seed\n$random_seq\n";

exit(0);

sub random_seq {
   my ($num) = @_;
   my @nuc = qw/ A C G T /;
   my $seq = '';
   for (1 .. $num){
      my $rand_ind = int(rand(scalar(@nuc)));
      $seq .= $nuc[$rand_ind];
   }
   return($seq);
}

exit(0);
