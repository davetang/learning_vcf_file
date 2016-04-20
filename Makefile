all:
	./install.sh
	bpipe run pipeline.groovy

clean_tmp:
	rm -rf commandlog.txt .bpipe *.fq *.fa* *.sai

clean:
	rm -rf commandlog.txt .bpipe *.fq *.fa* *.sai *.bam *.bai *mutated.txt *.vcf *.log aln.* bcftools-1.3.tar.bz2 GenomeAnalysisTK-3.5.tar.bz2 htslib-1.3.tar.bz2 samtools-1.3.tar.bz2 GenomeAnalysisTK* snpEff snpEff_latest_core.zip bpipe* bwa htslib-1.3 samtools-1.3 bcftools-1.3 resources picard vt test_31.* aln_rg.vcf.idx vcflib

