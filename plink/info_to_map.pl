#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "Usage: <infile.info>\n";
my $infile = shift or die $usage;

my $chr = '';

if ($infile =~ /^(\d+)\_*.*\.info$/){
   $chr = $1;
} else {
   die "Could not work out chromosome from $infile\n";
}

open(IN, '<', $infile) || die "Could not open $infile: $!\n";
while(<IN>){
   chomp;
   my ($id, $bp) = split(/\t/);
   print join("\t", $chr, $id, 0, $bp), "\n";
}
close(IN);

exit(0);

__END__

