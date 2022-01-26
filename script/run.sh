#!/bin/bash

java -Xmx1g -jar snpEff/SnpSift.jar concordance -v aln_consensus.vt.renamed.vcf aln_rg.vt.vcf > concordance_by_variant.txt

