// Variables
//
// size of reference sequence
REF_SIZE=1000000
// seed
SEED=31
// basename
BASE="test_${SEED}"
// name of reference sequence
REF="${BASE}.fa"
// name of mutated reference
REF_MUT="${BASE}_mutated.fa"
// ~1% mutation
MUT_PC=0.01
// 1 million reads
READ_NO=1000000
// 100 bp reads
READ_LEN=100
// inner mate distance
INNER_DIST=400
// mpileup depth
DEPTH=300
// number of threads to use
THREAD=8

// Programs
TOOLS="../tools"
BWA="$TOOLS/bin/bwa"
BGZIP="$TOOLS/bin/bgzip"
SAMTOOLS="$TOOLS/bin/samtools"
BCFTOOLS="$TOOLS/bin/bcftools"
FREEBAYES="$TOOLS/bin/freebayes"
GATK="$TOOLS/gatk"

// produce() will create the file, which can be referred to as $output
// and this output file will serve as the input file in the next step of the pipeline,
// which can be referred to as $input
random_ref = {
   produce("$REF"){
      // Usage: generate_random_seq.pl <bp> <seed>
      exec "../../script/generate_random_seq.pl $REF_SIZE $SEED > $output"
   }
}

// $input will be $output from the previous step
// as defined in Bpipe.run {}
index_ref = {
   produce("*.amb", "*.ann", "*.bwt", "*.pac", "*.sa"){
      exec "$BWA index $input"
   }
   forward input
}

create_dict = {
   produce("${BASE}.dict"){
      exec "$GATK CreateSequenceDictionary --REFERENCE $REF --OUTPUT $output"
   }
   forward input
}

// The from statement reshapes the inputs to be the most recent output file(s) matching the given pattern for the following block.
// This is useful when a task needs an input that was produced earlier in the pipeline than the previous stage, or
// other similar cases where your inputs don't match the defaults that Bpipe assumes.
mutate_ref = {
   // def my_output = input;
   from("fa"){
      produce("$REF_MUT"){
         exec "../../script/mutate_fasta.pl $input $MUT_PC $SEED > $output"
      }
   }
}

random_read = {
   produce("*_1.fq.gz", "*_2.fq.gz"){
      // Usage: random_paired_end.pl <infile.fa> <read length> <number of pairs> <inner mate distance> <seed>
      exec "../../script/random_paired_end.pl $input1 $READ_LEN $READ_NO $INNER_DIST $SEED"
   }
}

bwa_align = {
   produce("aln.sam"){
      exec "$BWA mem -t $THREAD -R \"@RG\\tID:test\\tSM:test\\tPL:test\" $REF $input1 $input2 > $output"
   }
}

sam_to_bam = {
   produce("aln.bam"){
      exec "$SAMTOOLS view -@ $THREAD -bS $input | $SAMTOOLS sort -@ $THREAD - -o $output"
   }
}

index_bam = {
   transform("bam") to ("bam.bai") {
      exec "$SAMTOOLS index $input.bam"
   }
   forward input
}

// Call variants as per https://samtools.github.io/bcftools/howtos/variant-calling.html
call_variant = {
   produce("aln.bt.vcf"){
      exec "$BCFTOOLS mpileup -d $DEPTH -O v -f $REF $input | $BCFTOOLS call -mv -O v -o $output"
   }
}

haplotype_caller = {
   produce("aln.hc.vcf"){
      exec "$GATK HaplotypeCaller --reference $REF --input $input --output $output"
   }
}

freebayes = {
   produce("aln.fb.vcf"){
      exec "$FREEBAYES -f $REF $input > $output"
   }
}

bgzip = {
   produce("*.vcf.gz"){
      exec "$BGZIP $input"
   }
}

Bpipe.run {
   random_ref + [ create_dict,index_ref ] + mutate_ref + random_read + bwa_align + sam_to_bam + index_bam + [ call_variant+bgzip,haplotype_caller+bgzip,freebayes+bgzip ]
}

