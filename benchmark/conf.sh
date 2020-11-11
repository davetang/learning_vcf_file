#!/usr/bin/env bash

set -euo pipefail

# number of threads to use
threads=8

# set seed
seed=31

# size of reference sequence
ref_size=1000000

# sample name
sample=test_31

# name of reference sequence
ref=${sample}.fa

# name of mutated reference
ref_mut=${sample}_mutated.fa

# percentage of mutations
mut_pc=0.01

# read specficiations
read_no=1000000
read_len=100
inner_dist=400

fastq1=l${read_len}_n${read_no}_d${inner_dist}_${seed}_1.fq.gz
fastq2=l${read_len}_n${read_no}_d${inner_dist}_${seed}_2.fq.gz

gatk_path=../bin/gatk-4.1.9.0/gatk

