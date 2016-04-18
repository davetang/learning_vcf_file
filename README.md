VCF files
---------

Natural selection relies on three conditions:

1. There must be genetic variation among species
2. The genetic variation must be heritable
3. The genetic variation results in differing fitness

The _de facto_ file format for representing genetic variation is the Variant Call Format (VCF). A good starting point for learning about the VCF is this [poster](http://vcftools.sourceforge.net/VCF-poster.pdf). The binary equivalent of a VCF file is a BCF file, akin to the SAM and BAM format. BCFtools is used to view and manipulate VCF/BCF files. I have included an example BCF file (aln_consensus.bcf) in this repository to demonstrate the various utilities of BCFtools. If you are interested in how this file was generated refer to [Creating a test file](#creating-a-test-file).

# Usage

Typing ```bcftools``` without any parameters will output the usage and the subcommands.

~~~~{.bash}
bcftools 

Program: bcftools (Tools for variant calling and manipulating VCFs and BCFs)
Version: 1.3 (using htslib 1.3)

Usage:   bcftools [--version|--version-only] [--help] <command> <argument>

Commands:

 -- Indexing
    index        index VCF/BCF files

 -- VCF/BCF manipulation
    annotate     annotate and edit VCF/BCF files
    concat       concatenate VCF/BCF files from the same set of samples
    convert      convert VCF/BCF files to different formats and back
    isec         intersections of VCF/BCF files
    merge        merge VCF/BCF files files from non-overlapping sample sets
    norm         left-align and normalize indels
    plugin       user-defined plugins
    query        transform VCF/BCF into user-defined formats
    reheader     modify VCF/BCF header, change sample names
    view         VCF/BCF conversion, view, subset and filter VCF/BCF files

 -- VCF/BCF analysis
    call         SNP/indel calling
    consensus    create consensus sequence by applying VCF variants
    cnv          HMM CNV calling
    filter       filter VCF/BCF files using fixed thresholds
    gtcheck      check sample concordance, detect sample swaps and contamination
    roh          identify runs of autozygosity (HMM)
    stats        produce VCF/BCF stats

 Most commands accept VCF, bgzipped VCF, and BCF with the file type detected
 automatically even when streaming from a pipe. Indexed VCF and BCF will work
 in all situations. Un-indexed VCF and BCF and streams will work in most but
 not all situations.

~~~~

# Viewing a BCF file

Use the appropriately named subcommand ```view```.

~~~~{.bash}
bcftools view aln_consensus.bcf | grep -v "^#" | head
1000000 58      .       AT      A       77.4563 .       INDEL;IDV=57;IMF=1;DP=57;VDB=1.20228e-08;SGB=-0.693136;MQ0F=0;AF1=1;AC1=2;DP4=0,0,35,0;MQ=60;FQ=-139.526        GT:PL   1/1:118,105,0
1000000 68      .       CTTTT   CTTT    70.4562 .       INDEL;IDV=68;IMF=1;DP=68;VDB=7.54492e-06;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,51,0;MQ=60;FQ=-188.527        GT:PL   1/1:111,154,0
1000000 225     .       CTT     CT      169.457 .       INDEL;IDV=78;IMF=0.928571;DP=84;VDB=0.0449154;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,79,0;MQ=60;FQ=-272.528   GT:PL   1/1:210,238,0
1000000 336     .       A       G       221.999 .       DP=112;VDB=0.756462;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,102,0;MQ=60;FQ=-281.989    GT:PL   1/1:255,255,0
1000000 378     .       T       C       221.999 .       DP=101;VDB=0.704379;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,99,0;MQ=60;FQ=-281.989     GT:PL   1/1:255,255,0
1000000 451     .       AGG     AGGG    214.458 .       INDEL;IDV=127;IMF=0.969466;DP=131;VDB=0.0478427;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,87,42;MQ=60;FQ=-289.528 GT:PL   1/1:255,255,0
1000000 915     .       G       GC      214.458 .       INDEL;IDV=179;IMF=0.913265;DP=196;VDB=0.929034;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,89,101;MQ=60;FQ=-289.528 GT:PL   1/1:255,255,0
1000000 1009    .       G       C       221.999 .       DP=203;VDB=0.259231;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,94,101;MQ=60;FQ=-281.989    GT:PL   1/1:255,255,0
1000000 1062    .       ATT     AT      214.458 .       INDEL;IDV=187;IMF=0.958974;DP=195;VDB=0.244824;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,97,92;MQ=60;FQ=-289.528  GT:PL   1/1:255,255,0
1000000 1207    .       T       G       221.999 .       DP=177;VDB=0.628515;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,94,79;MQ=60;FQ=-281.989     GT:PL   1/1:255,255,0
~~~~

# BCF to VCF

Use the ```convert``` subcommand.

~~~~{.bash}
bcftools convert -O v -o aln_consensus.vcf aln_consensus.bcf

# we can use cat to view the file
cat aln_consensus.vcf | grep -v "^#" | head
1000000 58      .       AT      A       77.4563 .       INDEL;IDV=57;IMF=1;DP=57;VDB=1.20228e-08;SGB=-0.693136;MQ0F=0;AF1=1;AC1=2;DP4=0,0,35,0;MQ=60;FQ=-139.526        GT:PL   1/1:118,105,0
1000000 68      .       CTTTT   CTTT    70.4562 .       INDEL;IDV=68;IMF=1;DP=68;VDB=7.54492e-06;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,51,0;MQ=60;FQ=-188.527        GT:PL   1/1:111,154,0
1000000 225     .       CTT     CT      169.457 .       INDEL;IDV=78;IMF=0.928571;DP=84;VDB=0.0449154;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,79,0;MQ=60;FQ=-272.528   GT:PL   1/1:210,238,0
1000000 336     .       A       G       221.999 .       DP=112;VDB=0.756462;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,102,0;MQ=60;FQ=-281.989    GT:PL   1/1:255,255,0
1000000 378     .       T       C       221.999 .       DP=101;VDB=0.704379;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,99,0;MQ=60;FQ=-281.989     GT:PL   1/1:255,255,0
1000000 451     .       AGG     AGGG    214.458 .       INDEL;IDV=127;IMF=0.969466;DP=131;VDB=0.0478427;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,87,42;MQ=60;FQ=-289.528 GT:PL   1/1:255,255,0
1000000 915     .       G       GC      214.458 .       INDEL;IDV=179;IMF=0.913265;DP=196;VDB=0.929034;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,89,101;MQ=60;FQ=-289.528 GT:PL   1/1:255,255,0
1000000 1009    .       G       C       221.999 .       DP=203;VDB=0.259231;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,94,101;MQ=60;FQ=-281.989    GT:PL   1/1:255,255,0
1000000 1062    .       ATT     AT      214.458 .       INDEL;IDV=187;IMF=0.958974;DP=195;VDB=0.244824;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,97,92;MQ=60;FQ=-289.528  GT:PL   1/1:255,255,0
1000000 1207    .       T       G       221.999 .       DP=177;VDB=0.628515;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,94,79;MQ=60;FQ=-281.989     GT:PL   1/1:255,255,0
~~~~

# Filtering for different types of mutations

The ```view``` subcommand lets you select specific types of variants.

## SNPs

~~~~{.bash}
bcftools view -v snps aln_consensus.bcf | grep -v "^#" | head
1000000 336     .       A       G       221.999 .       DP=112;VDB=0.756462;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,102,0;MQ=60;FQ=-281.989    GT:PL   1/1:255,255,0
1000000 378     .       T       C       221.999 .       DP=101;VDB=0.704379;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,99,0;MQ=60;FQ=-281.989     GT:PL   1/1:255,255,0
1000000 1009    .       G       C       221.999 .       DP=203;VDB=0.259231;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,94,101;MQ=60;FQ=-281.989    GT:PL   1/1:255,255,0
1000000 1207    .       T       G       221.999 .       DP=177;VDB=0.628515;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,94,79;MQ=60;FQ=-281.989     GT:PL   1/1:255,255,0
1000000 1281    .       C       A       221.999 .       DP=154;VDB=0.286069;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,66,80;MQ=60;FQ=-281.989     GT:PL   1/1:255,255,0
1000000 1405    .       A       T       221.999 .       DP=203;VDB=0.0898873;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,104,89;MQ=60;FQ=-281.989   GT:PL   1/1:255,255,0
1000000 1669    .       G       C       221.999 .       DP=191;VDB=0.656207;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,108,73;MQ=60;FQ=-281.989    GT:PL   1/1:255,255,0
1000000 1775    .       C       A       221.999 .       DP=225;VDB=0.413906;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,101,115;MQ=60;FQ=-281.989   GT:PL   1/1:255,255,0
1000000 2036    .       T       A       221.999 .       DP=193;VDB=0.227246;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,83,98;MQ=60;FQ=-281.989     GT:PL   1/1:255,255,0
1000000 2180    .       G       C       221.999 .       DP=211;VDB=0.123382;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,97,105;MQ=60;FQ=-281.989    GT:PL   1/1:255,255,0
~~~~

## INDELs

~~~~{.bash}
bcftools view -v indels aln_consensus.bcf | grep -v "^#" | head
1000000 58      .       AT      A       77.4563 .       INDEL;IDV=57;IMF=1;DP=57;VDB=1.20228e-08;SGB=-0.693136;MQ0F=0;AF1=1;AC1=2;DP4=0,0,35,0;MQ=60;FQ=-139.526        GT:PL   1/1:118,105,0
1000000 68      .       CTTTT   CTTT    70.4562 .       INDEL;IDV=68;IMF=1;DP=68;VDB=7.54492e-06;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,51,0;MQ=60;FQ=-188.527        GT:PL   1/1:111,154,0
1000000 225     .       CTT     CT      169.457 .       INDEL;IDV=78;IMF=0.928571;DP=84;VDB=0.0449154;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,79,0;MQ=60;FQ=-272.528   GT:PL   1/1:210,238,0
1000000 451     .       AGG     AGGG    214.458 .       INDEL;IDV=127;IMF=0.969466;DP=131;VDB=0.0478427;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,87,42;MQ=60;FQ=-289.528 GT:PL   1/1:255,255,0
1000000 915     .       G       GC      214.458 .       INDEL;IDV=179;IMF=0.913265;DP=196;VDB=0.929034;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,89,101;MQ=60;FQ=-289.528 GT:PL   1/1:255,255,0
1000000 1062    .       ATT     AT      214.458 .       INDEL;IDV=187;IMF=0.958974;DP=195;VDB=0.244824;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,97,92;MQ=60;FQ=-289.528  GT:PL   1/1:255,255,0
1000000 1278    .       TA      TAA     214.458 .       INDEL;IDV=144;IMF=0.929032;DP=155;VDB=0.252598;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,65,80;MQ=60;FQ=-289.528  GT:PL   1/1:255,255,0
1000000 1328    .       AT      A       129.457 .       INDEL;IDV=177;IMF=0.988827;DP=179;VDB=1.83715e-25;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,36,18;MQ=60;FQ=-197.527       GT:PL   1/1:170,163,0
1000000 1380    .       TA      TAA     214.458 .       INDEL;IDV=180;IMF=0.957447;DP=188;VDB=5.28227e-08;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,80,68;MQ=60;FQ=-289.528       GT:PL   1/1:255,255,0
1000000 1449    .       GT      G       214.458 .       INDEL;IDV=210;IMF=0.972222;DP=216;VDB=0.783773;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,101,109;MQ=60;FQ=-289.528        GT:PL   1/1:255,255,0
~~~~

# Extracting INFO field/s

The VCF has various information fields; use the ```query``` subcommand to extract specific field/s.

~~~~{.bash}
bcftools query -f 'DP=%DP\tAF1=%AF1\tAC1=%AC1\tMQ=%MQ\n' aln_consensus.bcf | head
DP=57   AF1=1   AC1=2   MQ=60
DP=68   AF1=1   AC1=2   MQ=60
DP=84   AF1=1   AC1=2   MQ=60
DP=112  AF1=1   AC1=2   MQ=60
DP=101  AF1=1   AC1=2   MQ=60
DP=131  AF1=1   AC1=2   MQ=60
DP=196  AF1=1   AC1=2   MQ=60
DP=203  AF1=1   AC1=2   MQ=60
DP=195  AF1=1   AC1=2   MQ=60
DP=177  AF1=1   AC1=2   MQ=60
~~~~

Combining with the ```view``` subcommand:

~~~~{.bash}
bcftools view -v snps aln_consensus.bcf | bcftools query -f 'DP=%DP\tAF1=%AF1\tAC1=%AC1\tMQ=%MQ\n' - | head
DP=112  AF1=1   AC1=2   MQ=60
DP=101  AF1=1   AC1=2   MQ=60
DP=203  AF1=1   AC1=2   MQ=60
DP=177  AF1=1   AC1=2   MQ=60
DP=154  AF1=1   AC1=2   MQ=60
DP=203  AF1=1   AC1=2   MQ=60
DP=191  AF1=1   AC1=2   MQ=60
DP=225  AF1=1   AC1=2   MQ=60
DP=193  AF1=1   AC1=2   MQ=60
DP=211  AF1=1   AC1=2   MQ=60
~~~~

# Creating a test file

The ```aln_consensus.bcf``` file was created from a simple pipeline. Firstly a random reference sequence was generated; genetic variants are created by modifying the reference sequence, i.e. introducing mutations, into a mutated copy and sequence reads were derived from the mutated reference sequence. Lastly, the reads were mapped back to the original non-mutated reference sequence. The ```pipeline.groovy``` file contains the pipeline, which is written in [Groovy](http://www.groovy-lang.org/) and processed by Bpipe. I have a [blog post](http://davetang.org/muse/2015/06/04/paired-end-alignment-using-bpipe/) that provides more information.

To create ```aln_consensus.bcf```, simply clone this repository and type ```make```.

~~~~{.bash}
git clone https://github.com/davetang/learning_vcf_file.git
make
~~~~

This will download and install all the necessary programs from online and run the pipeline.

## Adjusting parameters

All the variables are defined in ```pipeline.groovy```, which can be adjusted.

~~~~{.java}
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
~~~~

## Consensus caller

~~~~{.bash}
bcftools call -c -o aln_consensus.bcf -O b aln.bcf
~~~~

# Using GATK for calling variants

We'll use another variant caller to call variants and compare them to the variants called by BCFtools. The [HaplotypeCaller](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. Firstly, [download](https://www.broadinstitute.org/gatk/download/) and extract GATK; you'll need to register an account and to agree to the terms and conditions.

~~~~{.bash}
tar -xjf GenomeAnalysisTK-3.5.tar.bz2 
~~~~

Then we need to setup Picard to prepare our reference fasta file:

~~~~{.bash}
git clone https://github.com/broadinstitute/picard.git
cd picard
git clone https://github.com/samtools/htsjdk.git
cd htsjdk
# install ant on Debian/Ubuntu
# sudo apt-get install ant
ant htsjdk-jar
cd ..
ant -lib lib/ant package-commands
cd ..
~~~~

Some necessary steps before running HaplotypeCaller:

~~~~{.bash}
java -jar picard/dist/picard.jar CreateSequenceDictionary R=test_31.fa O=test_31.dict
samtools faidx test_31.fa

# add read groups to the BAM file
java -jar picard/dist/picard.jar AddOrReplaceReadGroups \
INPUT=aln.bam \
OUTPUT=aln_rg.bam \
RGLB=test \
RGPL=illumina \
RGPU=test \
RGSM=test

# check out the header
# to see the read groups we added
samtools view -H aln_rg.bam 
@HD     VN:1.5  SO:coordinate
@SQ     SN:1000000      LN:1000000
@RG     ID:1    LB:test PL:illumina     SM:test PU:test
@PG     ID:bwa  PN:bwa  VN:0.7.13-r1126 CL:bwa/bwa mem test_31.fa l100_n1000000_d400_31_1.fq l100_n1000000_d400_31_2.fq

# index
samtools index aln_rg.bam
~~~~

Now to call variants:

~~~~{.bash}
java -Xmx4G -jar GenomeAnalysisTK.jar -R test_31.fa -T HaplotypeCaller -I aln_rg.bam -o aln_rg.vcf
~~~~

# Comparing VCF files

How many variants were called using BCFtools?

~~~~{.bash}
# convert BCF to VCF
bcftools convert -O v -o aln_consensus.vcf aln_consensus.bcf

# count
cat aln_consensus.vcf | grep -v "^#" | wc -l
9704
~~~~

How many variants using HaplotypeCaller?

~~~~{.bash}
cat aln_rg.vcf | grep -v "^#" | wc -l
9875
~~~~

## Decompose and normalise

Despite the VCF being a standard, there are still differences between VCF files. To ensure the VCF files are unified, we'll use the ```vt``` program to decompose and normalise the variants. For more information, refer to this [blog post](http://davetang.org/muse/2015/12/16/getting-acquainted-analysing-dna-sequencing-data/).

~~~~{.bash}
# download and compile
git clone https://github.com/atks/vt.git
cd vt
make
make test
cd ..

# decompose and normalise
vt/vt decompose -s aln_consensus.vcf | vt normalize -r test_31.fa - > aln_consensus.vt.vcf

# this step doesn't do anything because
# the GATK variant file is already decomposed and normalised
vt/vt decompose -s aln_rg.vcf | vt normalize -r test_31.fa - > aln_rg.vt.vcf
~~~~

## SnpSift

We can use [SnpSift](http://snpeff.sourceforge.net/SnpSift.html) to compare VCF files; I have a [blog post](http://davetang.org/muse/2015/08/26/vcf-concordance/) with more information.

~~~~{.bash}
# download
wget http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
~~~~

~~~~{.bash}
# run SnpSift
java -Xmx1g -jar \
snpEff/SnpSift.jar concordance \
-v aln_consensus.vt.vcf aln_rg.vt.vcf \
> concordance_by_variant.txt
~~~~

# Useful links

A very useful thread on SEQanswers on learning about the VCF format: <http://seqanswers.com/forums/showthread.php?t=9345>

Useful tutorial on VCFs files from the 1000 Genomes Project Page: <http://www.1000genomes.org/node/101>

The author of ANNOVAR writes about VCF files: <http://annovar.openbioinformatics.org/en/latest/articles/VCF/>

[Encoding Structural Variants in VCF](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/VCF%20(Variant%20Call%20Format)%20version%204.0/encoding-structural-variants) version 4.0

