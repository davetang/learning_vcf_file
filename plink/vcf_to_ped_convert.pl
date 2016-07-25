#!/usr/bin/env perl

=head1 LICENSE

	Copyright (c) 1999-2011 The European Bioinformatics Institute and
	Genome Research Limited.  All rights reserved.

	This software is distributed under a modified Apache license.
	For license details, please see

	http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

	Please email comments or questions to the 1000 genomes project information 
	list at <info@1000genomes.org>

=head1	AUTHOR

  	Ian Streeter (info@1000genomes.org)

=cut


use Getopt::Long;
use Net::FTP;
use Env qw( @PATH );
use List::Util qw (first max);

use strict;
use warnings;

my @populations;
my $sample_panel;
my $ftp_host;
my $vcf;
my $region;
my $tabix;
my $output_ped;
my $output_info;
my $output_dir;
my $help;
my $max_maf = 1;
my $min_maf = 0;
my $base_format = 'number';

GetOptions('population=s' => \@populations,
            'vcf=s' => \$vcf,
            'sample_panel_file=s' => \$sample_panel,
            'region=s' => \$region,
            'tabix=s' => \$tabix,
            'output_ped=s' => \$output_ped,
            'output_info=s' => \$output_info,
            'output_dir=s' => \$output_dir,
            'help!' => \$help,
            'max_maf=s' => \$max_maf,
            'min_maf=s' => \$min_maf,
            'base_format=s' => \$base_format,
            );

if ($help) {
    exec('perldoc', $0);
}

die("required arguments: vcf, sample_panel_file, region, population") if (! $vcf || ! $sample_panel || ! $region || ! @populations);
die("$output_dir is not a directory") if ($output_dir && ! -d $output_dir);

die("base_format must be 'number' or 'letter'") if ($base_format !~ /num/i && $base_format !~ /let/i);
my %base_codes = $base_format =~ /num/i ? ('A' => 1,   'C' => 2,   'G' => 3,   'T' => 4)
                                        : ('A' => 'A', 'C' => 'C', 'G' => 'G', 'T' => 'T');

my $is_compressed = $vcf =~ /\.b?gz(ip)?$/;
if ($is_compressed) {
  $tabix ||= first {-x $_} map {"$_/tabix"} @PATH;
  die("cannot find executable $tabix") if (! -x $tabix);
}
die("remote vcf file must be compressed by bgzip") if (!$is_compressed && $vcf =~ /ftp:\/\//);

my ($region_chromosome, $region_start, $region_end);
if ($region =~ /^(\w+):(\d+)-(\d+)$/) {
  ($region_chromosome, $region_start, $region_end) = ($1, $2, $3);
}
else {
  $region_chromosome = $region;
}

if (! $output_ped) {
    $output_ped = "$region.ped";
    $output_ped =~ s{:}{_};
}
if (! $output_info) {
    $output_info = "$region.info";
    $output_info =~ s{:}{_};
}
if ($output_dir) {
    $output_ped = $output_dir . '/' . $output_ped;
    $output_ped =~ s{//}{/}g;

    $output_info = $output_dir . '/' . $output_info;
    $output_info =~ s{//}{/}g;
}

my %individuals;
my @markers;
my %genotypes;

get_individuals();
get_markers_genotypes();

print_info();
print_ped();

print "Created ".$output_info." and ".$output_ped."\n";



sub get_markers_genotypes {

    my $vcf_opener = $is_compressed ? "$tabix -h $vcf $region |" : "<$vcf";
    open my $VCF, $vcf_opener
        or die("cannot open vcf $!");

    my %column_indices;
    my $found_chromosome = 0;

    LINE:
    while (my $line = <$VCF>) {
        next LINE if ($line =~ /^\#\#/);
        chomp $line;
        my @columns  = split(/\t/, $line);

        if ($line =~ /^\#/) {
            foreach my $i (0..$#columns) {
                $column_indices{$columns[$i]} = $i;
            }
            next LINE;
        }

        my ($chromosome, $position, $name, $ref_allele, $alt_alleles) = @columns;
        if ($chromosome ne $region_chromosome) {
          last LINE if ($found_chromosome);
          next LINE;
        }
        $found_chromosome = 1;

        next LINE if (defined $region_start && $position < $region_start);
        last LINE if (defined $region_end && $position > $region_end);

        my @allele_codes = map {$base_codes{$_} || 0} $ref_allele, (split(/,/, $alt_alleles));
        #my @allele_codes = map {($number_format ? $base_codes{$_} : $_) || 0} $ref_allele, (split(/,/, $alt_alleles));
        next LINE if ((scalar grep {$_} @allele_codes) < 2);

        my %marker_genotypes;
        my %alleles_present;
        my $total_alleles = 0;
        foreach my $population (keys %individuals) {
            INDIVIDUAL:
            foreach my $individual (@{$individuals{$population}}) {
                next INDIVIDUAL if (! $column_indices{$individual});
                my $genotype_string = $columns[ $column_indices{$individual} ];
                if ($genotype_string =~ /(\d+)(?:\/|\|)(\d+)/) {
                  my @genotype_codes = ($allele_codes[$1], $allele_codes[$2]);

                  foreach my $allele_code (grep {$_} @genotype_codes) {
                    $alleles_present{$allele_code} ++;
                    $total_alleles ++;
                  }
                  $marker_genotypes{$population}{$individual} = \@genotype_codes;
                }
                else {
                  $marker_genotypes{$population}{$individual} = [0,0];
                }
            }
        }

        next LINE if ((scalar keys %alleles_present) < 2);

        my $major_allele_frequency = (max values %alleles_present) / $total_alleles;
        next LINE if ((1-$major_allele_frequency) < $min_maf);
        next LINE if ((1-$major_allele_frequency) > $max_maf);

        foreach my $population (keys %marker_genotypes) {
            foreach my $individual (keys %{$marker_genotypes{$population}}) {
                push(@{$genotypes{$population}{$individual}}, $marker_genotypes{$population}{$individual});
            }
        }

        if ($name eq '.') {
            $name = "$chromosome:$position";
        }
        push(@markers, [$name,$position]);

    }

    close $VCF;

    if ($is_compressed) {
      my $exit_status = $? >>8;
      die("tabix exited with status $exit_status") if $exit_status;
    }

    return;
}

sub print_ped {

    open my $FILE, '>', $output_ped
        or die "cannot open $output_ped $!";
    foreach my $population (keys %genotypes) {
        my $pedigree_counter = 1;
        foreach my $individual (keys %{$genotypes{$population}}) {
            my $pedigree = $population . '_' . $pedigree_counter;
            print $FILE join("\t", $pedigree, $individual, 0, 0, 0, 0,);
            foreach my $genotype_codes (@{$genotypes{$population}->{$individual}}) {
                print $FILE "\t", $genotype_codes->[0], ' ', $genotype_codes->[1];
            }
            print $FILE "\n";
            $pedigree_counter ++;
        }
    }
    close $FILE;
    return;
}

sub print_info {

    open my $FILE, '>', $output_info
        or die "cannot open $output_info $!";
    foreach my $marker (@markers) {
        print $FILE join("\t", @$marker), "\n";
    }
    close $FILE;
    return;
}




sub get_individuals {

    my @sample_panel_lines;

    if ($sample_panel =~ /ftp:\/\/([\w.]+)(\/\S+)/) {
        my $ftp_host = $1;
        my $path = $2;

        my $ftp = Net::FTP->new($ftp_host);
        $ftp->login or die('Cannot login ' , $ftp->message);

        my $sample_panel_content;
        open my $PANEL, '>', \$sample_panel_content;
        $ftp->get($path, $PANEL) or die ('could not $sample_panel ' , $ftp->message);
        $ftp->quit;
        close $PANEL;

        @sample_panel_lines = split(/\n/, $sample_panel_content);
    }
    else {
        open my $FILE, '<', $sample_panel
            or die("cannot open $sample_panel $!");
        @sample_panel_lines = <$FILE>;
        close $FILE;
    }

    my %allowed_pops_hash;
    foreach my $pop (@populations) {
        $allowed_pops_hash{$pop} = 1;
        $individuals{$pop} = [];
    }

    foreach my $line (@sample_panel_lines) {
        my ($individual, $population) = split(/\s+/, $line);
        if ($allowed_pops_hash{$population}) {
            push(@{$individuals{$population}}, $individual);
        }
    }
}



##################################################################################################################################

=pod

=head1 NAME 

	vcf_to_ped_converter.pl

=head1 SYNOPSIS

        This script prepares the input files for Haploview	

=head1	VERSION

	1.0

=head1	REQUIRED ARGUMENTS

	-vcf		    Path to a locally or remotely accessible vcf file.
                            The vcf file must be compressed by bgzip and indexed by tabix if it is a remote file.
                            The vcf format is a tab format for presenting variation sites and 
			    genotypes data and is described at http://vcftools.sourceforge.net/specs.html.
                            This tool takes both vcf4.0 and vcf4.1 format files.
	-sample_panel_file  Path to a locally or remotely accessible sample panel file, listing all individuals (first column)
                            and their population (second column)
	-region		    Chromosomal region in the format of chr:start-end (e.g. 1:1000000-100500) or chr (e.g. 1)
	-population         A population name, which must appear in the second column of the sample panel file.
                            Can be specified more than once for multiple populations.

=head1	OPTIONAL ARGUMENTS

	-tabix		    Path to the tabix executable; default is to search the path for 'tabix'
                            tabix is not required if the vcf file is uncompressed and locally accessible
	-output_ped	    Name of the output ped file (linkage pedigree file);
                            default is region.ped (e.g. 1_100000-100500.ped)
        -output_info        Name of the output info file (marker information file);
                            default is region.info (e.g. 1_1000000-100500.info)
        -output_dir         Name of a directory to place the output_ped and output_info files
        -min_maf            Only include variations with a minor allele_frequency greater than or equal to this value
        -max_maf            Only include variations with a minor allele_frequency less than or equal to this value
        -base_format        Either 'letter' or 'number'. Genotypes in the ped file can be coded either ACGT or 1-4
                            where 1=A, 2=C, 3=G, T=4.  Default is 'number'.
	-help		    Print out help menu
			
=head1	OUTPUT FILES

        The file formats of the linkage pedigree and marker information files are described at
        http://http://www.broadinstitute.org/science/programs/medical-and-population-genetics/haploview/input-file-formats-0
        Both files are needed as input for Haploview.

=head1 EXAMPLE

perl ~/ReseqTrack/scripts/variation_data/vcf_to_ped_converter.pl -vcf ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr13.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.vcf.gz -sample_panel_file ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.sample_panel -region 13:32889611-32973805 -population GBR -population FIN -min_maf 0.1
