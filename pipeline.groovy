SEED=31
REF_SIZE=1000000
REF="test_" + "$SEED" + ".fa"
REF_MUT="test_mutated.fa"
REF_MUT_LOG="test_mutated.log"
MUT_PC=0.01
//READ_NO=300000
READ_NO=1000000
READ_LEN=100
INNER_DIST=400

// produce() will create the file, which can be referred to as $output
// and this output file will serve as the input file, which can be
// referred to as $input, in the next step of the pipeline
random_ref = {
   produce("$REF"){
      // Usage: generate_random_seq.pl <bp> <seed>
      exec "script/generate_random_seq.pl $REF_SIZE $SEED > $output"
   }
}

index_ref = {
   exec "bwa index $input"
}

mutate_ref = {
   // Usage: ./mutate_fasta.pl <infile.fa> <mutation percent> <seed>
   exec "script/mutate_fasta.pl $REF $MUT_PC $SEED > $REF_MUT 2> $REF_MUT_LOG"
}

random_read = {
   // Usage: random_paired_end.pl <infile.fa> <read length> <number of pairs> <inner mate distance> <seed>
   exec "script/random_paired_end.pl $REF_MUT $READ_LEN $READ_NO $INNER_DIST $SEED"
}

bwa_align = {
   from(glob("*.fq")) {
      transform("sai","sai","bam"){
         multi "bwa aln $REF $input1.fq > $output1",
               "bwa aln $REF $input2.fq > $output2"

         exec """
            bwa sampe $REF $output1 $output2 $input1.fq $input2.fq | 
            samtools view -bS - | 
            samtools sort - $output.bam.prefix
         """
      }
   }
}

bwa_index = {
   from(glob("*.bam")){
      exec "samtools index $input.bam"
   }
}

// samtools mpileup
// Collects summary information in the input BAMs,
// computes the likelihood of data given each possible,
// genotype and stores the likelihoods in the BCF format.
mpileup = {
   from(glob("*.bam")){
      transform("bcf"){
         // -g          generate BCF output (genotype likelihoods)
         // -f FILE     faidx indexed reference sequence file [null]
         exec "samtools mpileup -g -f $REF $input.bam > $output"
      }
   }
}

consensus = {
   from(glob("*.bcf")){
      exec "bcftools call -c -o $output.bcf -O b $input.bcf"
   }
}

Bpipe.run { random_ref + index_ref + mutate_ref + random_read + bwa_align + bwa_index + mpileup + consensus }
