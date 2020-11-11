#!/usr/bin/env bash

set -euo pipefail

source ./conf.sh

rm -f ${sample}.fa* ${sample}.bam* ${fastq1} ${fastq2} ${sample}_mutated.fa

>&2 echo Done!

