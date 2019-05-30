#!/usr/bin/env perl

# Takes a fasta file (expects only 1 sequence!, i.e. not a multi-fasta file) as input and outputs random SNVs

use strict;
use warnings;

my $usage = "Usage: $0 <infile.fa> <num> <seed>\n";
my $infile = shift or die $usage;
my $num = shift or die $usage;
my $seed = shift or die $usage;

my $seq = read_fasta($infile);
my $header = get_fasta_header($infile);
my $chr = get_chr($header);
random_snv($seq, $num);

exit(0);

sub get_fasta_header {
   my ($infile) = @_;
   open(IN, '<', $infile) || die "Could not open $infile: $!\n";
   chomp(my $header = <IN>);
   close(IN);
   return($header);
}

sub get_chr {
   my ($header) = @_;
   $header =~ s/>(\w+)\s.*/$1/;
   return($header);
}

sub read_fasta {
   my ($infile) = @_;
   my $seq = '';
   open(IN, '<', $infile) || die "Could not open $infile: $!\n";
   while(<IN>){
      chomp;
      next if /^>/;
      $seq .= $_;
   }
   close(IN);
   return($seq);
}

sub random_snv {
   my ($seq, $num) = @_;
   my $seq_len = length($seq);

   # set seed for reproducibility
   srand($seed);
   my @rand_pos = map { int(rand($seq_len)) } ( 1 .. $num );
   @rand_pos = sort({$a <=> $b} @rand_pos);

   # substr is zero based
   foreach my $pos (@rand_pos){

      my $ref = substr($seq, $pos-1, 1);
      my $alt = mutate_base($ref);
      print join("\t", $chr, $pos-1 ,$pos, $ref, $alt), "\n";

   }

}

sub mutate_base{
   my ($nuc) = @_;
   my %sub = ('A' => ['C','G','T'],
              'C' => ['A','G','T'],
              'G' => ['A','C','T'],
              'T' => ['A','C','G']);
   return $sub{$nuc}[rand(3)];
}

