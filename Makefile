all:
	./install.sh
	bpipe run pipeline.groovy

clean_tmp:
	rm -rf commandlog.txt .bpipe *.fq *.fa* *.sai

clean:
	rm -rf commandlog.txt .bpipe *.fq *.fa* *.sai *.bam *.bai *mutated.txt *.vcf *.log aln.*
