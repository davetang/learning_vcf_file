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

   my @track = 1 .. length($s);
   my $point = 0;
   my $insert = 0;
   my $delete = 0;

   for (1 .. $$pc * length($s)){
      #pick a random base to mutate
      my $rand_ind = int(rand(length($s)));
      #what type of mutation should we introduce
      my $type = int(rand(scalar(@mutation)));
      $type = $mutation[$type];

      if ($type eq 'point'){
         ++$point;
         my $base = substr($s, $rand_ind, 1);
         my $rand_nuc = $base;
         while($rand_nuc eq $base){
            $rand_nuc = $nuc[int(rand(scalar(@nuc)))];
         }
         substr($s, $rand_ind, 1, $rand_nuc);
         $track[$rand_ind] = "$type: $base -> $rand_nuc";
#         print "$type on ", $rand_ind+1, " $base -> $rand_nuc\n";
      }

      elsif ($type eq 'delete'){
         ++$delete;
         my $a = substr($s, 0, $rand_ind-1);
         my $b = substr($s, $rand_ind, length($s));
         $s = $a . $b;
         splice @track, $rand_ind, 1;
      }

      elsif ($type eq 'insert'){
         ++$insert;
         my $rand_nuc = $nuc[int(rand(scalar(@nuc)))];
         my $a = substr($s, 0, $rand_ind);
         my $b = substr($s, $rand_ind, length($s));
         $s = $a . $rand_nuc . $b;
         splice @track, $rand_ind, 0, "$type: $rand_nuc";
#         print "$type on ", $rand_ind+1, " $rand_nuc\n";
      }

   }
   for (my $i=0; $i<scalar(@track); ++$i){
      warn $i+1,"\t$track[$i]\n" if $track[$i] !~ /^\d/;
   }
   warn "Point: $point\nDelete: $delete\nInsert: $insert\n";
   return($s);
}
