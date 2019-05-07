#!/usr/bin/env perl

use strict;
use warnings;
# use Data::Dumper;

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
   my $line = join("\t", $chr, $pos, $ref, @$alt);

   for my $gt (keys %{$$x{gtypes}}){

      # converts encoding into actual alleles
      my ($al1, $sep, $al2) = $vcf->parse_alleles($x, $gt);

      my $allele = "\t$al1$sep$al2";
      $line .= $allele;

   }

   print "$line\n";

}

exit(0);

