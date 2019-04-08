#!/bin/bash

# Download GATK
# wget https://github.com/broadinstitute/gatk/releases/download/4.1.1.0/gatk-4.1.1.0.zip
# unzip gatk-4.1.1.0.zip

ref_size=1000000
ref_name=test.fa
mut_name=test_mut.fa
mut_log=test_mut.log
mut_perc=0.01
read_num=1000000
read_len=100
inner_dist=400
seed=31
read1=l${read_len}_n${read_num}_d${inner_dist}_${seed}_1.fq
read2=l${read_len}_n${read_num}_d${inner_dist}_${seed}_2.fq
bam=aln.bam
bcf=aln.bcf

# run parameterised R Markdown file to generate reference sequence and read pairs
Rscript -e "rmarkdown::render('variant.Rmd', params = list(ref_size = \"$ref_size\", ref_name = \"$ref_name\", mut_name = \"$mut_name\", mut_log = \"$mut_log\", mut_perc = \"$mut_perc\", read_num = \"$read_num\", read_len = \"$read_len\", inner_dist = \"$inner_dist\", seed = \"$seed\"))"

# use conda environment that has necessary programs installed
source activate learning_vcf

# index reference
bwa index $ref_name

# align and create sorted BAM file
bwa mem $ref_name $read1 $read2 | samtools sort -o $bam -

# add read groups to the BAM file
java -Xmx4G -jar gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar AddOrReplaceReadGroups \
--INPUT aln.bam \
--OUTPUT aln.bam.tmp \
--RGLB test_lib \
--RGPL illumina \
--RGPU test_plat \
--RGSM test_sample

# move BAM file with read groups
mv $bam.tmp $bam

# index BAM file
samtools index $bam

# mpileup
samtools mpileup -g -f $ref_name $bam > $bcf
# using bcftools does not produce the same results
# bcftools mpileup -O b -f $ref_name $bam > $bcf

# consensus
bcftools call -v -c -o bcftools.vcf -O v $bcf

# freebayes
freebayes -f $ref_name $bam > freebayes.vcf

# GATK
java -Xmx4G -jar gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar CreateSequenceDictionary -R $ref_name
java -Xmx4G -jar gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar HaplotypeCaller --input $bam --output gatk.vcf --reference $ref_name

# decompose and normalise SNVs
vt decompose -s bcftools.vcf | vt normalize -r $ref_name - > bcftools.vt.vcf
vt decompose -s freebayes.vcf | vt normalize -r $ref_name - > freebayes.vt.vcf
vt decompose -s gatk.vcf | vt normalize -r $ref_name - > gatk.vt.vcf

# move final results into result folder
mv bcftools.vt.vcf ../result
mv freebayes.vt.vcf ../result
mv gatk.vt.vcf ../result

echo Done

