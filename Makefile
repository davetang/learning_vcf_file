all:
	bpipe run pipeline.groovy

clean_tmp:
	rm -rf commandlog.txt .bpipe *.fq *.fa* *.sai

clean:
	rm -rf commandlog.txt .bpipe *.fq *.fa* *.sai *.bam *.bai *.bcf *mutated.txt *.vcf *.log
