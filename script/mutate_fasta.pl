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

#original sequence
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
   my $s = $$seq;
   my @nuc = qw/ A C G T /;
   my @mutation = qw/insert delete point/;
   for (1 .. $$pc * length($s)){
      #pick a random base
      my $rand_ind = int(rand(length($s)));
      #what type of mutation should we introduce
      my $type = int(rand(scalar(@mutation)));
      $type = $mutation[$type];

      if ($type eq 'point'){
         my $base = substr($s, $rand_ind, 1);
         my $rand_nuc = $base;
         while($rand_nuc eq $base){
            $rand_nuc = $nuc[int(rand(scalar(@nuc)))];
         }
         substr($s, $rand_ind, 1, $rand_nuc);
      }

      elsif ($type eq 'delete'){
         my $a = substr($s, 0, $rand_ind-1);
         my $b = substr($s, $rand_ind, length($s));
         $s = $a . $b;
      }

      elsif ($type eq 'insert'){
         my $rand_nuc = $nuc[int(rand(scalar(@nuc)))];
         my $a = substr($s, 0, $rand_ind);
         my $b = substr($s, $rand_ind, length($s));
         $s = $a . $rand_nuc . $b;
      }

   }
   return($s);
}
