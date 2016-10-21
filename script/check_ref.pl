#!/usr/bin/env perl
#
# Check the reference allele of a VCF file to a fasta file
#
# The script assumes that the reference genome name
# is used to name the fasta file and only contains
# alpha-numeric characters
#
# To obtain hg19 (and assuming using Linux):
# wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
# wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
# chmod 755 twoBitToFa 
# twoBitToFa hg19.2bit hg19.fa
# gzip hg19.fa
# rm hg19.2bit
#

use strict;
use warnings;

my $usage = "Usage: $0 <genome.fa> <infile.vcf>\n";
my $genome = shift or die $usage;
my $infile = shift or die $usage;

warn("Storing $genome into memory\n");
my %genome = store_genome($genome);
warn("Finished storing $genome into memory\n");

if ($infile =~ /\.vcf$/){
   open(IN,'<',$infile) || die "Could not open $infile: $!\n";
} elsif ($infile =~ /\.vcf\.gz$/){
   open(IN,'-|',"gunzip -c $infile") || die "Could not open $infile: $!\n";
} else {
   die "Could not recognise the format of $infile; is it a VCF or bgzipped VCF file?\n";
}

# tallies
my $mm = 0;
my $total = 0;

while(<IN>){

   chomp;
   next if /^#/;
   # $pos is one base
   my ($chr, $pos, $id, $ref, $alt, @rest) = split(/\t/);

   if (exists $genome{$chr}){
      ++$total;

      # obtain sequence from genome fasta
      my $l = length($ref);
      my $base = substr($genome{$chr}, $pos-1, $l);
      $base = uc($base);

      # compare sequence from fasta to reported sequence
      if ($base ne $ref){
         ++$mm;
         warn "Mismatch at position $pos: $infile has $ref but $genome has $base\n";
      }

   } else {

      die "Missing $chr in $genome\n";

   }

}
close(IN);

warn("There were $mm mismatches in $total variants\n");

exit(0);

sub store_genome {

   my ($infile) = @_;
   # warn "Storing $infile\n";
   my %genome = ();

   if ($infile =~ /\.fa$/){
      open(IN,'<',$infile) || die "Could not open $infile: $!\n";
   } elsif ($infile =~ /\.fa\.gz$/){
      open(IN,'-|',"gunzip -c $infile") || die "Could not open $infile: $!\n";
   } else {
      die "Could not recognise the file format of $infile\n";
   }

   my $chr = '';
   while(<IN>){
      chomp;
      if (/^>(.*)$/){
         $chr = $1;
      } else {
         $genome{$chr} .= $_;
      }
   }
   close(IN);

   # warn "Done storing $genome\n";
   return(%genome);

}

