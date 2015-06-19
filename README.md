Learning about VCF files
------------------------

To learn about VCF files, I created a pipeline that generates a random reference, which becomes mutated (SNVs are introduced), from which random reads are derived from, and subsequently mapped back to the original non-mutated random reference. The ```pipeline.groovy``` file contains the pipeline, which is processed by Bpipe. For more information take a look at my [blog post](http://davetang.org/muse/2015/06/04/paired-end-alignment-using-bpipe/).

Simply type make in the root directory of this repository to run the pipeline, assuming that bpipe is installed.

Adjusting parameters
--------------------

All the variables are defined in ```pipeline.groovy```.

```
SEED=31
REF_SIZE=1000000
REF="test_" + "$SEED" + ".fa"
REF_MUT="test_mutated.fa"
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
1000000 2       .       G       <X>     0       .       DP=1;I16=1,0,0,0,41,1681,0,0,60,3600,0,0,0,0,0,0;QS=1,0;MQ0F=0  PL      0,3,41
1000000 3       .       T       <X>     0       .       DP=2;I16=2,0,0,0,82,3362,0,0,120,7200,0,0,1,1,0,0;QS=1,0;MQ0F=0 PL      0,6,75
1000000 4       .       C       <X>     0       .       DP=2;I16=2,0,0,0,82,3362,0,0,120,7200,0,0,3,5,0,0;QS=1,0;MQ0F=0 PL      0,6,75
1000000 5       .       A       <X>     0       .       DP=6;I16=6,0,0,0,226,8546,0,0,360,21600,0,0,5,13,0,0;QS=1,0;MQ0F=0      PL      0,18,148
1000000 6       .       C       <X>     0       .       DP=7;I16=7,0,0,0,281,11311,0,0,420,25200,0,0,11,29,0,0;QS=1,0;MQ0F=0    PL      0,21,169
1000000 7       .       A       <X>     0       .       DP=7;I16=7,0,0,0,286,11686,0,0,420,25200,0,0,18,58,0,0;QS=1,0;MQ0F=0    PL      0,21,171
1000000 8       .       G       <X>     0       .       DP=10;I16=10,0,0,0,395,15655,0,0,600,36000,0,0,25,101,0,0;QS=1,0;MQ0F=0 PL      0,30,192
1000000 9       .       A       <X>     0       .       DP=11;I16=11,0,0,0,443,17899,0,0,660,39600,0,0,35,161,0,0;QS=1,0;MQ0F=0 PL      0,33,200
1000000 10      .       A       <X>     0       .       DP=12;I16=12,0,0,0,482,19402,0,0,720,43200,0,0,46,242,0,0;QS=1,0;MQ0F=0 PL      0,36,204
1000000 11      .       G       <X>     0       .       DP=12;I16=12,0,0,0,492,20172,0,0,720,43200,0,0,58,346,0,0;QS=1,0;MQ0F=0 PL      0,36,206
```

Converting to VCF
-----------------

```
bcftools convert -O v -o l100_n1000000_d400_31_1.vcf l100_n1000000_d400_31_1.bcf
cat l100_n1000000_d400_31_1.vcf | grep -v "^#" | head
1000000 2       .       G       <X>     0       .       DP=1;I16=1,0,0,0,41,1681,0,0,60,3600,0,0,0,0,0,0;QS=1,0;MQ0F=0  PL      0,3,41
1000000 3       .       T       <X>     0       .       DP=2;I16=2,0,0,0,82,3362,0,0,120,7200,0,0,1,1,0,0;QS=1,0;MQ0F=0 PL      0,6,75
1000000 4       .       C       <X>     0       .       DP=2;I16=2,0,0,0,82,3362,0,0,120,7200,0,0,3,5,0,0;QS=1,0;MQ0F=0 PL      0,6,75
1000000 5       .       A       <X>     0       .       DP=6;I16=6,0,0,0,226,8546,0,0,360,21600,0,0,5,13,0,0;QS=1,0;MQ0F=0      PL      0,18,148
1000000 6       .       C       <X>     0       .       DP=7;I16=7,0,0,0,281,11311,0,0,420,25200,0,0,11,29,0,0;QS=1,0;MQ0F=0    PL      0,21,169
1000000 7       .       A       <X>     0       .       DP=7;I16=7,0,0,0,286,11686,0,0,420,25200,0,0,18,58,0,0;QS=1,0;MQ0F=0    PL      0,21,171
1000000 8       .       G       <X>     0       .       DP=10;I16=10,0,0,0,395,15655,0,0,600,36000,0,0,25,101,0,0;QS=1,0;MQ0F=0 PL      0,30,192
1000000 9       .       A       <X>     0       .       DP=11;I16=11,0,0,0,443,17899,0,0,660,39600,0,0,35,161,0,0;QS=1,0;MQ0F=0 PL      0,33,200
1000000 10      .       A       <X>     0       .       DP=12;I16=12,0,0,0,482,19402,0,0,720,43200,0,0,46,242,0,0;QS=1,0;MQ0F=0 PL      0,36,204
1000000 11      .       G       <X>     0       .       DP=12;I16=12,0,0,0,492,20172,0,0,720,43200,0,0,58,346,0,0;QS=1,0;MQ0F=0 PL      0,36,206
```

Consensus caller
----------------

```
bcftools call -c -o l100_n1000000_d400_31_1_consensus.bcf -O b l100_n1000000_d400_31_1.bcf
```

Cross-check with introduced variants
------------------------------------

```
head test_31_mutated.txt
107     A       C
151     T       C
174     A       T
227     T       C
368     T       A
393     G       C
571     C       G
758     A       G
1023    A       G
1064    T       G

bcftools view -v snps l100_n1000000_d400_31_1_consensus.bcf | grep -v "^#" | head
1000000 107     .       A       C       221.999 .       DP=98;VDB=0.0346205;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,92,0;MQ=60;FQ=-281.989     GT:PL   1/1:255,255,0
1000000 151     .       T       C       221.999 .       DP=92;VDB=0.391205;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,85,0;MQ=60;FQ=-281.989      GT:PL   1/1:255,255,0
1000000 174     .       A       T       221.999 .       DP=88;VDB=0.741413;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,81,0;MQ=60;FQ=-270.989      GT:PL   1/1:255,244,0
1000000 227     .       T       C       221.999 .       DP=84;VDB=0.036411;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,79,0;MQ=60;FQ=-264.989      GT:PL   1/1:255,238,0
1000000 368     .       T       A       221.999 .       DP=99;VDB=0.626025;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,96,0;MQ=60;FQ=-281.989      GT:PL   1/1:255,255,0
1000000 393     .       G       C       221.999 .       DP=96;VDB=0.395084;SGB=-0.693147;MQ0F=0;AF1=1;AC1=2;DP4=0,0,92,0;MQ=60;FQ=-281.989      GT:PL   1/1:255,255,0
1000000 571     .       C       G       221.999 .       DP=186;VDB=0.582778;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,96,81;MQ=60;FQ=-281.989     GT:PL   1/1:255,255,0
1000000 758     .       A       G       221.999 .       DP=197;VDB=0.499294;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,92,100;MQ=60;FQ=-281.989    GT:PL   1/1:255,255,0
1000000 1023    .       A       G       221.999 .       DP=209;VDB=0.297086;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,104,96;MQ=60;FQ=-281.989    GT:PL   1/1:255,255,0
1000000 1064    .       T       G       221.999 .       DP=203;VDB=0.562939;SGB=-0.693147;MQSB=1;MQ0F=0;AF1=1;AC1=2;DP4=0,0,97,98;MQ=60;FQ=-281.989     GT:PL   1/1:255,255,0
```
