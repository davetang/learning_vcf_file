## README

A dumb benchmark of the variant calling tools: SAMtools mpileup/BCFtools, FreeBayes, and HaplotypeCaller. Why create a dumb benchmark? Sometimes it is useful to create simplified examples to observe how tools work.

Use Conda to create the working environment using `environment.yml`.

```bash
conda env create
```

GATK needs to be prepared separately; use `../bin/broad_tools.sh` to download and unzip the necessary files.

The script `pipeline.sh` will run the entire pipeline with 8 threads; make sure that the Conda environment has been activated and that `../bin/gatk-4.1.9.0/gatk` exists. Modify `conf.sh` accordingly if you wish to use less threads or have `gatk` downloaded elsewhere.

```bash
conda activate variant_tools
./pipeline.sh
```

Three separate bgzipped VCF files and their indexes (along with other files) should be produced at the end of the script.

Use `benchmark.pl` to compare called variants with the artificially introduced variants, i.e. real variants.

```bash
# variants called by BCFtools
./benchmark.pl -v test_31.bt.vcf.gz -l test_31_mutation.log | cut -f1 | sort | uniq -c
9379 Correct
 521 FP
  10 Miscalled
 602 Missed

./benchmark.pl -v test_31.fb.vcf.gz -l test_31_mutation.log | cut -f1 | sort | uniq -c
8789 Correct
 521 FP
 303 Miscalled
 899 Missed

./benchmark.pl -v test_31.hc.vcf.gz -l test_31_mutation.log | cut -f1 | sort | uniq -c
9439 Correct
 525 FP
  11 Miscalled
 541 Missed
```

FP are variants that don't exist in the real set of variants. Miscalled are variants at the correct position but incorrectly called variants. Missed are variants that were not detected, i.e. false negatives.

Most variants reported as "Missed" are repeated nucleotides.

```bash
./benchmark.pl -v test_31.hc.vcf.gz -l test_31_mutation.log | grep Missed | cut -f3-4 | sort | uniq -c | sort -k1rn
  79 G  GG
  73 C  CC
  68 GG G
  62 A  AA
  62 CC C
  57 TT T
  54 AA A
  46 T  TT
   4 GT G
   3 A  C
   3 AC A
   3 G  A
   2 A  AG
   2 AT A
   2 C  A
   2 T  C
   2 T  G
   2 TA T
   1 A  G
   1 A  T
   1 AG A
   1 C  T
   1 CA C
   1 CG C
   1 G  C
   1 G  GA
   1 G  T
   1 GC G
   1 T  A
   1 T  TA
   1 T  TC
   1 TC T
   1 TG T
```

In reality, these variants were not missed but just the reported position is different. Consider the example below.

```bash
Correct    7351    A     T
Correct    7371    T     A
FP         7453    CG    C
Missed     7454    GG    G
Correct    7504    A     C
Correct    7565    A     G
```

The variant was correctly identified as a G deletion, however, the first G was reported as deleted and not the second real variant. This should be counted as a correctly called variant since the G deletion was called.

