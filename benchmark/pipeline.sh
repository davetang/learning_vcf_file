#!/usr/bin/env bash

set -euo pipefail

source ./conf.sh

rm -f ${sample}.fa* ${sample}.dict ${sample}.bam* ${fastq1} ${fastq2} ${sample}_mutated.fa ${sample}_mutation.log
rm -f ${sample}.bt.vcf.gz* ${sample}.fb.vcf.gz* ${sample}.hc.vcf.gz*

../script/generate_random_seq.pl ${ref_size} ${seed} > ${ref}
bwa index ${ref}

../script/mutate_fasta.pl ${ref} ${mut_pc} ${seed} > ${ref_mut}
../script/random_paired_end.pl ${ref_mut} ${read_len} ${read_no} ${inner_dist} ${seed}

bwa mem \
   -t ${threads} \
   -R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina\tPU:${sample}\tLB:${sample}" \
   ${ref} \
   ${fastq1} \
   ${fastq2} |
   samtools sort -@ ${threads} -O BAM |\
   tee ${sample}.bam |\
   samtools index - ${sample}.bam.bai

# see https://samtools.github.io/bcftools/howtos/variant-calling.html
bcftools mpileup \
   --threads ${threads} \
   -Ov \
   -f ${ref} \
   ${sample}.bam |\
   bcftools call \
   --threads ${threads} \
   -mv \
   -Ov |\
   vt decompose -s - |\
   vt normalize -r ${ref} - \
	> ${sample}.bt.vcf

bgzip -f ${sample}.bt.vcf
tabix -p vcf ${sample}.bt.vcf.gz

freebayes -f ${ref} ${sample}.bam |\
   vt decompose -s - |\
   vt normalize -r ${ref} - \
   > ${sample}.fb.vcf

bgzip -f ${sample}.fb.vcf
tabix -p vcf ${sample}.fb.vcf.gz

${gatk_path} \
   CreateSequenceDictionary \
   -R ${ref} \
   -O ${sample}.dict

samtools faidx ${ref}

${gatk_path} \
   HaplotypeCaller \
   -R ${ref} \
   -I ${sample}.bam \
   -O ${sample}.hc.vcf.gz

>&2 echo Done!

