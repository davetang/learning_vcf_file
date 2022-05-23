#!/usr/bin/env perl

# Takes a fasta file as input and randomly mutates bases

use strict;
use warnings;

my $usage = "Usage: $0 <infile.fa> <mutation percent> <seed> [max_indel_len (default: 10)]\n";
my $infile = shift or die $usage;
my $percent = shift or die $usage;
my $seed = shift or die $usage;

my $indel_max = 10;
my @nuc = qw/ A C G T /;
if (@ARGV > 0){
   $indel_max = shift;
}

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
   my $num = int(length($seq) * $pc);
   my %pos = pick_locus($seq, $num);
   my $outfile = $infile;
   $outfile =~ s/\.fa$/_mutation.log/;
   open(my $fh, '>', $outfile) || die "Could not open $outfile for writing: $!\n";

   foreach my $pos (sort {$b <=> $a} keys %pos){
      my $ref = substr($seq, $pos, 1);
      my $type = rand_mut(0.8, 0.1, 0.1);
      my $alt = $ref;

      if ($type eq 'SNV'){
         while($alt eq $ref){
            $alt = $nuc[int(rand(scalar(@nuc)))];
         }
         substr($seq, $pos, 1, $alt);
      } elsif ($type eq 'DEL'){
         my $len = rand_num($indel_max);
         $ref = substr($seq, $pos-1, $len + 1);
         $alt = substr($seq, $pos-1, 1);
         substr($seq, $pos, $len, '');
         $pos -= 1;
      } elsif ($type eq 'INS'){
         my $rand_nuc = rand_seq($indel_max);
         $alt = $ref . $rand_nuc;
         substr($seq, $pos, 1, $alt);
      } else {
         die "Unknown mutation\n";
      }
      # 1-based coordinates as per VCF files
      print $fh join("\t", $pos+1, $ref, $alt, $type), "\n";

   }
   close($fh);
   return($seq);
}

sub pick_locus {
   my ($seq, $num) = @_;
   my %pos = ();
   for (1 .. $num){
      my $r = int(rand(length($seq)));
      if (exists $pos{$r}){
         while (exists $pos{$r}){
            $r = int(rand(length($seq)));
         }
         $pos{$r} = 1;
      } else {
         $pos{$r} = 1;
      }
   }
   return(%pos);
}

sub rand_num {
   my ($max) = @_;
   int(rand($max))+1;
}

sub rand_seq {
   my ($max) = @_;
   my $len = rand_num($max);
   my $seq = '';
   for(1..$len){
      $seq .= $nuc[int(rand(scalar(@nuc)))];
   }
   return($seq);
}

sub rand_mut {
   my ($snv, $ins, $del) = @_;
   die "$snv + $ins + $del must equal 1\n" if $snv+$ins+$del != 1;
   my $rand = rand(1);
   if ($rand < $snv){
      return('SNV')
   } elsif ($rand < $snv+$ins){
      return('INS')
   } else {
      return('DEL')
   }
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

