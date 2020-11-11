#!/usr/bin/env perl

# Takes a fasta file as input and randomly mutates bases

use strict;
use warnings;

my $usage = "Usage: $0 <infile.fa> <mutation percent> <seed>\n";
my $infile = shift or die $usage;
my $percent = shift or die $usage;
my $seed = shift or die $usage;

# set seed for reproducibility
srand($seed);

my $seq = read_fasta($infile);
my $header = get_fasta_header($infile);
my $mutated = mutate_seq($seq, $percent);
$mutated = add_newline($mutated, 60);

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

   my ($seq, $pc) = @_;
   my @nuc = qw/ A C G T /;
   my @mutation = qw/insert delete point/;

   my $total_num = int(length($seq) * $pc);
   my $point_num = int($total_num * .8);
   my $del_num = int($total_num * .1);
   my $ins_num = int($total_num * .1);
   my %mut_log = ();

   %mut_log = pick_locus(\%mut_log, $point_num, "point");
   %mut_log = pick_locus(\%mut_log, $ins_num, "ins");
   %mut_log = pick_locus(\%mut_log, $del_num, "del");

   foreach my $pos (sort {$b <=> $a} keys %mut_log){
      my $type = $mut_log{$pos};
      my $ref = substr($seq, $pos, 1);

      if ($type eq 'point'){
         my $alt = $ref;
         while($alt eq $ref){
            $alt = $nuc[int(rand(scalar(@nuc)))];
         }
         $mut_log{$pos} = "$ref\t$alt";
         substr($seq, $pos, 1, $alt);
      } elsif ($type eq 'del'){
         my $prev_ref = substr($seq, $pos-1, 1);
         $mut_log{$pos} = "${prev_ref}$ref\t$prev_ref";
         substr($seq, $pos, 1, '');
      } elsif ($type eq 'ins'){
         my $rand_nuc = $nuc[int(rand(scalar(@nuc)))];
         my $alt = $ref . $rand_nuc;
         $mut_log{$pos} = "$ref\t$alt";
         substr($seq, $pos, 1, $alt);
      } else {
         die "Unknown mutation\n";
      }
   }
   my $outfile = $infile;
   $outfile =~ s/\.fa$/_mutation.log/;
   open(my $fh, '>', $outfile) || die "Could not open $outfile for writing: $!\n";
   foreach my $pos (sort {$a <=> $b} keys %mut_log){
      # use 1 bp coordinates to be compatible with VCF files
      my $pos_one_base = $pos + 1;
      my $mut = $mut_log{$pos};
      my ($ref, $alt) = split(/\t/, $mut);
      # if deletion, report coorindate of previous base
      if (length($ref) > length($alt)){
         print $fh "$pos\t$mut\n";
      } else {
         print $fh "$pos_one_base\t$mut\n";
      }

   }
   close($fh);
   return($seq);
}

sub pick_locus {
   my ($mut_log, $num, $type) = @_;
   for (1 .. $num){
      # pick a random base to mutate
      # if we already picked that base, repeat until we find a new one
      my $rand_ind = int(rand(length($seq)));
      if (exists $mut_log->{$rand_ind}){
         while (exists $mut_log->{$rand_ind}){
            $rand_ind = int(rand(length($seq)));
         }
         $mut_log->{$rand_ind} = $type;
      } else {
         $mut_log->{$rand_ind} = $type;
      }
   }
   return(%{$mut_log});
}

sub add_newline {
   my ($line, $pos) = @_;
   my $mod = '';
   for (my $i = 0; $i < length($line); $i += $pos){
      $mod .= substr($line, $i, $pos) . "\n";
   }
   return($mod);
}

__END__

