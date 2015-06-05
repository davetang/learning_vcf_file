#!/usr/bin/env perl

# Takes a fasta file as input and randomly mutates bases

use strict;
use warnings;

my $usage = "Usage: $0 <infile.fa> <mutation percent> <seed>\n";
my $infile = shift or die $usage;
my $mutation = shift or die $usage;
my $seed = shift or die $usage;

#set seed for reproducibility
srand($seed);

my $seq = read_fasta($infile);
my $header = get_fasta_header($infile);
my $mutated = mutate_seq(\$seq, \$mutation);

print "$header\n$mutated\n";

exit(0);

sub get_fasta_header {
   my ($infile) = @_;
   open(IN,'<',$infile) || die "Could not open $infile: $!\n";
   chomp(my $header = <IN>);
   close(IN);
   return($header);
}

sub read_fasta {
   my ($infile) = @_;
   my $seq = '';
   open(IN,'<',$infile) || die "Could not open $infile: $!\n";
   while(<IN>){
      chomp;
      next if /^>/;
      $seq .= $_;
   }
   close(IN);
   return($seq);
}

sub mutate_seq {
   my ($seq,$pc) = @_;
   my $outfile = $infile;
   $outfile =~ s/\.fa$/_mutated.txt/;
   open(OUT,'>',$outfile) || die "Could not open $outfile for writing: $!\n";
   my $s = $$seq;
   my @nuc = qw/ A C G T /;
   my %mutated = ();
   for (1 .. $$pc * length($s)){
      my $rand_ind = int(rand(length($s)));
      my $base = substr($s, $rand_ind, 1);
      my $rand_nuc = $base;
      while($rand_nuc eq $base){
         $rand_nuc = $nuc[int(rand(scalar(@nuc)))];
      }
      substr($s, $rand_ind, 1, $rand_nuc);
      $mutated{$rand_ind}->{'FROM'} = $base;
      $mutated{$rand_ind}->{'TO'} = $rand_nuc;
   }
   my @pos = sort {$a <=> $b} keys %mutated;
   foreach my $pos (@pos){
      print OUT join("\t", $pos+1, $mutated{$pos}->{'FROM'}, $mutated{$pos}->{'TO'}), "\n";
   }
   return($s);
}
