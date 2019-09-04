#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

# libraries are initialised during compile
BEGIN{
   my $libdir = $0;
   $libdir =~ s/\w+\.pl$//;
   push(@INC, $libdir);
}

# download Vcf module from VCFtools
# https://raw.githubusercontent.com/vcftools/vcftools/master/src/perl/Vcf.pm
use Vcf;

my $usage = "Usage: $0 <infile.vcf>\n";
my $infile = shift or die $usage;

my $vcf = Vcf->new(file => $infile);

$vcf->parse_header();
my @sample = $vcf->get_samples();

# I want to store variants by their chromosome and position on the chromosome; this will be the index
# the value will be a reference to the data structure used by the Vcf package
my %variant = ();

VARIANT: while (my $x = $vcf->next_data_hash()){

   # print Dumper \$x;
   # $VAR1 = \{
   #             'INFO' => {
   #                         'AC1' => '2',
   #                         'INDEL' => undef,
   #                         'MQ' => '60',
   #                         'FQ' => '-139.526',
   #                         'IDV' => '57',
   #                         'IMF' => '1',
   #                         'SGB' => '-0.693136',
   #                         'DP' => '57',
   #                         'MQ0F' => '0',
   #                         'DP4' => '0,0,35,0',
   #                         'AF1' => '1',
   #                         'VDB' => '1.20228e-08'
   #                       },
   #             'POS' => '58',
   #             'FILTER' => [
   #                           '.'
   #                         ],
   #             'CHROM' => '1000000',
   #             'ID' => '.',
   #             'gtypes' => {
   #                           'aln.bam' => {
   #                                          'GT' => '1/1',
   #                                          'PL' => '118,105,0'
   #                                        }
   #                         },
   #             'FORMAT' => [
   #                           'GT',
   #                           'PL'
   #                         ],
   #             'REF' => 'AT',
   #             'QUAL' => '77.4563',
   #             'ALT' => [
   #                        'A'
   #                      ]
   #           };

   my $chr = $$x{CHROM};
   my $pos = $$x{POS};
   my $ref = $$x{REF};
   my $alt = $$x{ALT};

   my $key = $chr . ':' . $pos;
   $variant{$key} = \$x;

   # if (exists $x->{"INFO"}

   # we expect normalised VCF files, i.e. only one variant per line
   if (scalar @$alt > 1){
      die "More than 1 alternate allele on line $.";
   }

   my $line = join("\t", $chr, $pos, $ref, @$alt);
   # ignore this section; I was trying to come up with different ways of storing the alleles
   # for my $gt (keys %{$$x{gtypes}}){
   # use sample array to keep same order each time
   for my $gt (@sample){

      # converts encoding into actual alleles
      my ($al1, $sep, $al2) = $vcf->parse_alleles($x, $gt);
      my $allele = "$al1$sep$al2";

      my $encoding = '';
      # "0" indicates two reference alleles
      if ($al1 eq $ref && $al2 eq $ref){
         $encoding = 0;
      # "2" indicates two alternate alleles
      } elsif ($al1 eq $alt->[0] && $al2 eq $alt->[0]){
         $encoding = 2;
      # "1" indicates heterozygous
      } elsif ($al1 eq $alt->[0] || $al2 eq $alt->[0]){
         $encoding = 1;
      # "3" indicates missing genotype
      } elsif ($al1 eq '.' && $al2 eq '.'){
         $encoding = 3;
      } else {
         die "Unexpected allele combination on line $.\n";
      }

      $line .= "\t$encoding";

   }

   # print "$line\n";

}

# print header
print join("\t", "CHROM", "POS", "REF", "ALT", @sample), "\n";

foreach my $index (keys %variant){
   # print Dumper $variant{$index};
   # foreach my $s (@sample){
   #    print join("\t", $index, $s, ${$variant{$index}}->{'gtypes'}->{$s}->{'GT'}), "\n";
   # }
   print_column($variant{$index});
}

sub print_column {
   my ($x) = @_;
   # print ${$x}->{"REF"}, "\n";

   my $chr = ${$x}->{"CHROM"};
   my $pos = ${$x}->{"POS"};
   my $ref = ${$x}->{"REF"};
   my $alt = ${$x}->{"ALT"};

   my $line = join("\t", $chr, $pos, $ref, @$alt);

   # we expect normalised VCF files, i.e. only one variant per line
   if (scalar @$alt > 1){
      die "More than 1 alternate allele on line $.";
   }
   foreach my $s (@sample){
      $line .= "\t${$x}->{'gtypes'}->{$s}->{'GT'}";
   }
   print "$line\n";

}

exit(0);

