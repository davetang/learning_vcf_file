#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('h:v:l:', \%opts);

if ($opts{'h'} ||
    !exists $opts{'v'} ||
    !exists $opts{'l'}
){
   usage();
}

my $vcf = $opts{'v'};
my $log = $opts{'l'};

my %mut = ();
my $fh;
open($fh, '<', $log) || die "Could not open $log: $!\n";
while(<$fh>){
   chomp;
   my ($pos, $ref, $alt) = split(/\t/);
   $mut{$pos}->{'REF'} = $ref;
   $mut{$pos}->{'ALT'} = $alt;
}
close($fh);

my $fh2;
if ($vcf =~ /\.gz$/){
   open($fh2, '-|', "gunzip -c $vcf") || die "Could not open $vcf: $!\n";
} else {
   open($fh2, '<', $vcf) || die "Could not open $vcf: $!\n";
}

while(<$fh2>){
   chomp;
   next if /^#/;
   my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @sample) = split(/\t/);
   if (exists $mut{$pos}){
      if ($mut{$pos}->{'REF'} eq $ref && $mut{$pos}->{'ALT'} eq $alt){
         print join("\t", "Correct", $pos, $ref, $alt), "\n";
      } else {
         print join("\t", "Miscalled", $pos, $ref, $alt, $mut{$pos}->{REF}, $mut{$pos}->{ALT}), "\n";
      }
      delete $mut{$pos};
   } else {
      print join("\t", "FP", $pos, $ref, $alt), "\n";
   }

}
close($fh2);

foreach my $pos (keys %mut){
   my $ref = $mut{$pos}->{'REF'};
   my $alt = $mut{$pos}->{'ALT'};
   print join("\t", "Missed", $pos, $ref, $alt), "\n";
}

sub usage {
print STDERR <<EOF;
Usage: $0 -v file -l file

Where:   -v         VCF file
         -l         mutation log file
         -h         this helpful usage message

EOF
exit();
}

