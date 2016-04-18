// Variables
REF_SIZE=1000000
REF="test_31.fa"
REF_MUT="test_mutated.fa"
REF_MUT_LOG="test_mutated.log"
MUT_PC=0.01
//READ_NO=300000
READ_NO=1000000
READ_LEN=100
INNER_DIST=400

// Programs
BWA="bwa/bwa"
SAMTOOLS="samtools-1.3/samtools"
BCFTOOLS="bcftools-1.3/bcftools"

// produce() will create the file, which can be referred to as $output
// and this output file will serve as the input file in the next step of the pipeline,
// which can be referred to as $input
random_ref = {
   produce("$REF"){
      // Usage: generate_random_seq.pl <bp> <seed>
      exec "script/generate_random_seq.pl $REF_SIZE 31 > $output"
   }
}

// $input will be $output from the previous step
// as defined in Bpipe.run {}
index_ref = {
   forward input
   produce("*.amb", "*.ann", "*.bwt", "*.pac", "*.sa"){
      exec "$BWA index $input"
   }
}

mutate_ref = {
   // def my_output = input;
   produce("$REF_MUT"){
      exec "script/mutate_fasta.pl $input $MUT_PC 31 > $output 2> $REF_MUT_LOG"
   }
}

random_read = {
   produce("*_1.fq", "*_2.fq"){
      // Usage: random_paired_end.pl <infile.fa> <read length> <number of pairs> <inner mate distance> <seed>
      exec "script/random_paired_end.pl $input $READ_LEN $READ_NO $INNER_DIST 31"
   }
}

bwa_align = {
   produce("aln.sam"){
      exec "$BWA mem $REF $input1 $input2 > $output"
   }
}

sam_to_bam = {
   produce("aln.bam"){
      exec "$SAMTOOLS view -bS $input | $SAMTOOLS sort - -o $output"
   }
}

index_bam = {
   transform("bam") to ("bam.bai") {
      exec "$SAMTOOLS index $input.bam"
   }
   forward input
}

// samtools mpileup
// Collects summary information in the input BAMs,
// computes the likelihood of data given each possible,
// genotype and stores the likelihoods in the BCF format.
mpileup = {
   transform("bcf"){
      // -g          generate BCF output (genotype likelihoods)
      // -f FILE     faidx indexed reference sequence file [null]
      exec "$SAMTOOLS mpileup -g -f $REF $input > $output"
   }
}

consensus = {
   produce("aln_consensus.bcf"){
      exec "$BCFTOOLS call -c -o $output -O b $input"
   }
}

Bpipe.run { random_ref + index_ref + mutate_ref + random_read + bwa_align + sam_to_bam + index_bam + mpileup + consensus }
