#!/bin/bash

# manually download GATK
# I have a copy in my src folder
cp ~/src/GenomeAnalysisTK-3.5.tar.bz2 .
tar -xjf GenomeAnalysisTK-3.5.tar.bz2

# Picard
git clone https://github.com/broadinstitute/picard.git
cd picard
git clone https://github.com/samtools/htsjdk.git
cd htsjdk
ant htsjdk-jar
cd ..
ant -lib lib/ant package-commands
cd ..

# vt
git clone https://github.com/atks/vt.git
cd vt
make
make test
cd ..

# SnpSift
wget http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip
unzip snpEff_latest_core.zip

# create sequence dictionary
java -jar picard/dist/picard.jar CreateSequenceDictionary R=test_31.fa O=test_31.dict

# faidx
samtools faidx test_31.fa

# add read groups to the BAM file
java -jar picard/dist/picard.jar AddOrReplaceReadGroups \
INPUT=aln.bam \
OUTPUT=aln_rg.bam \
RGLB=test \
RGPL=illumina \
RGPU=test \
RGSM=test

# index
samtools index aln_rg.bam

# call variants
java -Xmx4G -jar GenomeAnalysisTK.jar -R test_31.fa -T HaplotypeCaller -I aln_rg.bam -o aln_rg.vcf

# convert BCF to VCF
bcftools convert -O v -o aln_consensus.vcf aln_consensus.bcf

# decompose and normalise
vt/vt decompose -s aln_consensus.vcf | vt normalize -r test_31.fa - > aln_consensus.vt.vcf
vt/vt decompose -s aln_rg.vcf | vt normalize -r test_31.fa - > aln_rg.vt.vcf

