Learning about VCF files
------------------------

To learn about VCF files, I created a pipeline that generates a random reference, which becomes mutated (SNVs are introduced), from which random reads are derived from, and subsequently mapped back to the original non-mutated random reference. The ```pipeline.groovy``` file contains the pipeline, which is processed by Bpipe. For more information take a look at my [blog post](http://davetang.org/muse/2015/06/04/paired-end-alignment-using-bpipe/).

Simply type ```make``` in the root directory of this repository to run the pipeline. You will need to install bwa, samtools, bcftools, and Bpipe.

Adjusting parameters
--------------------

All the variables are defined in ```pipeline.groovy```.

```
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
```

Viewing the BCF file
-------------------

```
bcftools view l100_n1000000_d400_31_1.bcf | grep -v "^#" | head
```

Converting to VCF
-----------------

```
bcftools convert -O v -o l100_n1000000_d400_31_1.vcf l100_n1000000_d400_31_1.bcf
```

Consensus caller
----------------

```
bcftools call -c -o l100_n1000000_d400_31_1_consensus.bcf -O b l100_n1000000_d400_31_1.bcf
```

Filtering for different types of mutations
------------------------------------------

SNPs

```
bcftools view -v snps l100_n1000000_d400_31_1.consensus.bcf | grep -v "^#" | head
1000000 336     .       A       G       221.999 .       DP=112;VDB=0.756462;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,102,0;MQ=60;FQ=-281.989    GT:PL   1/1:255,255,0
1000000 378     .       T       C       221.999 .       DP=101;VDB=0.706367;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,99,0;MQ=60;FQ=-281.989     GT:PL   1/1:255,255,0
1000000 1009    .       G       C       221.999 .       DP=203;VDB=0.263498;SGB=-0.693147;MQSB=9.69427e-13;MQ0F=0;AF1=1;AC1=2;DP4=0,0,94,101;MQ=53;FQ=-281.989  GT:PL   1/1:255,255,0
1000000 1207    .       T       G       221.999 .       DP=177;VDB=0.647852;SGB=-0.693147;MQSB=0.989795;MQ0F=0;AF1=1;AC1=2;DP4=0,0,94,79;MQ=60;FQ=-281.989      GT:PL   1/1:255,255,0
1000000 1281    .       C       A       219.999 .       DP=145;VDB=0.150661;SGB=-0.693147;MQSB=0.220993;MQ0F=0;AF1=1;AC1=2;DP4=0,0,64,79;MQ=48;FQ=-281.989      GT:PL   1/1:253,255,0
1000000 1405    .       A       T       221.999 .       DP=196;VDB=0.126483;SGB=-0.693147;MQSB=0.374041;MQ0F=0;AF1=1;AC1=2;DP4=0,0,104,89;MQ=40;FQ=-281.989     GT:PL   1/1:255,255,0
1000000 1669    .       G       C       221.999 .       DP=191;VDB=0.649092;SGB=-0.693147;MQSB=0.00462056;MQ0F=0;AF1=1;AC1=2;DP4=0,0,108,73;MQ=57;FQ=-281.989   GT:PL   1/1:255,255,0
1000000 1775    .       C       A       221.999 .       DP=225;VDB=0.413906;SGB=-0.693147;MQSB=1.35318e-35;MQ0F=0;AF1=1;AC1=2;DP4=0,0,101,115;MQ=46;FQ=-281.989 GT:PL   1/1:255,255,0
1000000 2036    .       T       A       221.999 .       DP=193;VDB=0.227246;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,83,98;MQ=60;FQ=-281.989     GT:PL   1/1:255,255,0
1000000 2180    .       G       C       221.999 .       DP=211;VDB=0.127508;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,97,105;MQ=60;FQ=-281.989    GT:PL   1/1:255,255,0
```

Indels

```
bcftools view -v indels l100_n1000000_d400_31_1.consensus.bcf | grep -v "^#" | head
1000000 58      .       AT      A       77.4563 .       INDEL;IDV=49;IMF=1;DP=49;VDB=1.20228e-08;SGB=-0.693136;MQ0F=0;AF1=1;AC1=2;DP4=0,0,35,0;MQ=29;FQ=-139.526        GT:PL   1/1:118,105,0
1000000 68      .       CTTTT   CTTT    72.4562 .       INDEL;IDV=68;IMF=0.971429;DP=70;VDB=0.0087156;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,56,0;MQ=37;FQ=-203.527   GT:PL   1/1:113,169,0
1000000 225     .       CTT     CT      186.458 .       INDEL;IDV=77;IMF=0.916667;DP=84;VDB=0.1415;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,82,0;MQ=60;FQ=-281.528      GT:PL   1/1:227,247,0
1000000 451     .       AGG     AGGG    214.458 .       INDEL;IDV=126;IMF=0.940298;DP=134;VDB=0.113628;SGB=-0.693147;MQSB=1.53453e-19;MQ0F=0;AF1=1;AC1=2;DP4=0,0,90,43;MQ=52;FQ=-289.528     GT:PL   1/1:255,255,0
1000000 915     .       G       GC      214.458 .       INDEL;IDV=175;IMF=0.883838;DP=198;VDB=0.977609;SGB=-0.693147;MQSB=2.36671e-25;MQ0F=0;AF1=1;AC1=2;DP4=0,0,93,103;MQ=49;FQ=-289.528    GT:PL   1/1:255,255,0
1000000 1062    .       ATT     AT      214.458 .       INDEL;IDV=183;IMF=0.928934;DP=197;VDB=0.369482;SGB=-0.693147;MQSB=0.000483862;MQ0F=0;AF1=1;AC1=2;DP4=0,0,98,94;MQ=56;FQ=-289.528     GT:PL   1/1:255,255,0
1000000 1278    .       TA      TAA     214.458 .       INDEL;IDV=143;IMF=0.953333;DP=150;VDB=0.27816;SGB=-0.693147;MQSB=0.154055;MQ0F=0;AF1=1;AC1=2;DP4=0,0,67,79;MQ=48;FQ=-289.528GT:PL    1/1:255,255,0
1000000 1328    .       AT      A       128.457 .       INDEL;IDV=163;IMF=0.953216;DP=171;VDB=8.3785e-26;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,38,21;MQ=29;FQ=-212.527    GT:PL    1/1:169,178,0
1000000 1380    .       TA      TAA     214.458 .       INDEL;IDV=167;IMF=0.976608;DP=171;VDB=9.35461e-09;SGB=-0.693147;MQSB=0.999867;MQ0F=0;AF1=1;AC1=2;DP4=0,0,74,68;MQ=30;FQ=-289.528     GT:PL   1/1:255,255,0
1000000 1449    .       GT      G       214.458 .       INDEL;IDV=191;IMF=0.927184;DP=206;VDB=0.137289;SGB=-0.693147;MQSB=0.408104;MQ0F=0;AF1=1;AC1=2;DP4=0,0,92,100;MQ=48;FQ=-289.528       GT:PL   1/1:255,255,0
```

Using GATK for calling variants
-------------------------------

The [HaplotypeCaller](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. Firstly, [download GATK](https://www.broadinstitute.org/gatk/download/) and then unzip the file:


```
tar -xjf GenomeAnalysisTK-3.4-0.tar.bz2
```

Firstly we need to setup Picard to prepare our reference fasta file:

```
git clone https://github.com/broadinstitute/picard.git
cd picard
git clone https://github.com/samtools/htsjdk.git
cd htsjdk
ant htsjdk-jar
cd ..
ant -lib lib/ant package-commands
```

Some necessary steps before running HaplotypeCaller:

```
java -jar picard/dist/picard.jar CreateSequenceDictionary R= test_31.fa O= test_31.dict
samtools faidx test_31.fa
java -jar picard/dist/picard.jar AddOrReplaceReadGroups \
INPUT=l100_n1000000_d400_31_1.bam \
OUTPUT=l100_n1000000_d400_31_1_rg.bam \
RGLB=test \
RGPL=illumina \
RGPU=test \
RGSM=test
#index
samtools index l100_n1000000_d400_31_1_rg.bam
```

Now to call variants:

```
java -Xmx4G -jar GenomeAnalysisTK.jar -R test_31.fa -T HaplotypeCaller -I l100_n1000000_d400_31_1_rg.bam -o l100_n1000000_d400_31_1_hc.vcf
```

Comparing VCF files
-------------------

Firstly download [VCFtools](http://sourceforge.net/projects/vcftools/files/)

```
tar -xzf vcftools_0.1.12b.tar.gz 
cd vcftools_0.1.12b
make
```

How many variants using SAMtools?

```
cat l100_n1000000_d400_31_1.consensus.vcf | grep -v "^#" | awk '$5!="." {print}' | wc -l
    9892
```

How many variants using HaplotypeCaller?

```
cat l100_n1000000_d400_31_1_hc.vcf | grep -v "^#" | awk '$5!="." {print}' | wc -l
    9875
```

Compare VCF files using VCFtools:

```
bcftools convert -O v -o l100_n1000000_d400_31_1.consensus.vcf l100_n1000000_d400_31_1.consensus.bcf
a=l100_n1000000_d400_31_1.consensus.vcf
b=l100_n1000000_d400_31_1_hc.vcf
vcftools_0.1.12b/bin/vcftools --vcf $a --diff $b --out compare

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf l100_n1000000_d400_31_1.consensus.vcf
	--out compare
	--diff l100_n1000000_d400_31_1_hc.vcf

After filtering, kept 1 out of 1 Individuals
Comparing individuals in VCF files...
N_combined_individuals:	2
N_individuals_common_to_both_files:	0
N_individuals_unique_to_file1:	1
N_individuals_unique_to_file2:	1
Comparing sites in VCF files...
Found 6572 sites common to both files.
Found 996385 sites only in main file.
Found 0 sites only in second file.
Found 3303 non-matching overlapping sites.
After filtering, kept 1006260 out of a possible 1006260 Sites
Run Time = 3.00 seconds
```
Useful links
------------

A very useful thread on SEQanswers on learning about the VCF format: <http://seqanswers.com/forums/showthread.php?t=9345>

Useful tutorial on VCFs files from the 1000 Genomes Project Page: <http://www.1000genomes.org/node/101>

The author of ANNOVAR writes about VCF files: <http://annovar.openbioinformatics.org/en/latest/articles/VCF/>

[Encoding Structural Variants in VCF](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/VCF%20(Variant%20Call%20Format)%20version%204.0/encoding-structural-variants) version 4.0

