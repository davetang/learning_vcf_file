#!/usr/bin/env perl
#
# Create a minimalistic VCF file by Dave Tang
#
# Provide script with a five column tab-delimited file containing
# "chrom, start, end, reference base, alternate base"
# and a reference genome fasta file
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
# Validate the VCF file using http://vcftools.sourceforge.net/perl_module.html#vcf-validator
#

use strict;
use warnings;

my $usage = "Usage: $0 <genome.fa> <infile.tsv> > infile.vcf\n";
my $genome = shift or die $usage;
my $infile = shift or die $usage;

my $reference = '';
if ($genome =~ /(\w+)\./){
   $reference = $1;
} else {
   die "Could not create reference name\n";
}

my $sample = '';
if ($infile =~ /(\w+)\./){
   $sample = $1;
} else {
   die "Could not create sample name\n";
}

my %genome = store_genome($genome);

my %chromosome = ();
my @position = ();
my %position = ();

open(IN,'<',$infile) || die "Could not open $infile: $!\n";
VARIANT: while(<IN>){
   chomp;
   my ($chr, $start, $end, $ref, $alt) = split(/\t/);
   next if /^#/;
   next if /^$/;

   if (exists $genome{$chr}){
      my $id = "$chr:$start-$end";
      my $base = substr($genome{$chr}, $end-1, 1);
      $base = uc($base);

      if ($base ne $ref){
         warn "Mismatch at $id: $infile has $ref but $genome has $base\n";
         next VARIANT;
      }

      push(@position, $id);
      $position{$id}->{'REF'} = $ref;
      $position{$id}->{'ALT'} = $alt;

      unless (exists $chromosome{$chr}){
         $chromosome{$chr} = length($genome{$chr});
      }
   } else {
      die "Missing $chr in $genome\n";
   }

}
close(IN);

print <<EOF;
##fileformat=VCFv4.0
##reference=$reference
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOF

foreach my $chr (sort {$a cmp $b} keys %chromosome){
   ##contig=<ID=chr12,length=133851895,assembly=hg19>
   print "##contig=<ID=$chr,length=$chromosome{$chr},assembly=$reference>\n";
}

print join("\t", '#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',$sample),"\n";

foreach my $id (@position){
   #chr1    861368  .       CG      C       .       .       IC=19   GT      1/1
   my $ref = $position{$id}->{'REF'};
   my $alt = $position{$id}->{'ALT'};
   my ($chrom, $start, $end) = ('','','');
   if ($id =~ /^(\w+):\d+-(\d+)$/){
      $chrom = $1;
      $end = $2;
   }
   print join("\t", $chrom, $end, '.', $ref, $alt, '.', '.', '.', 'GT', '1/1'), "\n";
}

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

