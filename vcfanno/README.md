Testing vcfanno
===============

[vcfanno](https://github.com/brentp/vcfanno) is a tool that allows you to annotate unannotated VCF files; this is extremely handy. Just download the [relevant binary](https://github.com/brentp/vcfanno/releases).

Creating a test VCF file
========================

Using SNPs from the ESP.

~~~~{.bash}
# 127M file
wget -c http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz
tar -xzf ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz
parallel bgzip ::: *.vcf
parallel tabix -p vcf ::: *.vcf.gz
bcftools merge -o ESP6500SI-V2-SSA137.all.vcf.gz -O z *.vcf.gz

bcftools view -v snps ~/data/esp/ESP6500SI-V2-SSA137.all.vcf.gz | grep "^#" > header
bcftools view -v snps ~/data/esp/ESP6500SI-V2-SSA137.all.vcf.gz | grep -v "^#"  | shuf | head -100 > variants

cat header > test.vcf && cat variants | sort -k1,1V -k2,2n >> test.vcf
rm header variants
bgzip test.vcf
tabix -p vcf test.vcf.gz
~~~~

Running vcfanno
===============

Requires a configuration file; I created one called `conf.toml`.

~~~~{.bash}
# 4.1G file
wget -c ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz
wget -c ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz.tbi

cat conf.toml 
[[annotation]]
file="ExAC.r0.3.1.sites.vep.vcf.gz"
fields = ["AF", "AC_Het", "AC_Hom", "VQSLOD"]
ops=["first", "first", "first", "first"]

time vcfanno -p 16 conf.toml test.vcf.gz > test.annotated.vcf

=============================================
vcfanno version 0.1.1-alpha [built with go1.8beta1]

see: https://github.com/brentp/vcfanno
=============================================
vcfanno.go:114: found 4 sources from 1 files
vcfanno.go:220: annotated 100 variants in 0.23 seconds (439.0 / second)

real    0m0.350s
user    0m0.844s
sys     0m0.148s

zcat test.vcf.gz | grep -v "^#" | head -1
1       11105045        rs369923268     C       T       .       PASS    DBSNP=dbSNP_138;EA_AC=1,8597;AA_AC=0,4404;TAC=1,13001;MAF=0.0116,0.0,0.0077;GTS=TT,TC,CC;EA_GTC=0,1,4298;AA_GTC=0,0,2202;GTC=0,1,6500;GL=MASP2;CP=0;CG=-5.8;AA=C;CA=.;EXOME_CHIP=no;GWAS_PUBMED=.;FG=NM_139208.2:intron,NM_006610.3:intron;HGVS_CDNA_VAR=NM_139208.2:c.545-35G>A,NM_006610.3:c.544+420G>A;HGVS_PROTEIN_VAR=.,.;CDS_SIZES=NM_139208.2:558,NM_006610.3:2061;GS=.,.;PH=.,.;EA_AGE=1.2+/-3.3;AA_AGE=.;GRCh38_POSITION=1:11044988;DP=29

cat test.annotated.vcf | grep -v "^#" | head -1
1       11105045        rs369923268     C       T       .       PASS    DBSNP=dbSNP_138;EA_AC=1,8597;AA_AC=0,4404;TAC=1,13001;MAF=0.0116,0.0,0.0077;GTS=TT,TC,CC;EA_GTC=0,1,4298;AA_GTC=0,0,2202;GTC=0,1,6500;GL=MASP2;CP=0;CG=-5.8;AA=C;CA=.;EXOME_CHIP=no;GWAS_PUBMED=.;FG=NM_139208.2:intron,NM_006610.3:intron;HGVS_CDNA_VAR=NM_139208.2:c.545-35G>A,NM_006610.3:c.544+420G>A;HGVS_PROTEIN_VAR=.,.;CDS_SIZES=NM_139208.2:558,NM_006610.3:2061;GS=.,.;PH=.,.;EA_AGE=1.2+/-3.3;AA_AGE=.;GRCh38_POSITION=1:11044988;DP=29;AF=0.0001652;AC_Het=20;AC_Hom=0;VQSLOD=1.88

bgzip test.annotated.vcf
tabix -p vcf test.annotated.vcf.gz

bcftools query -f 'AF=%AF\tAC_Het=%AC_Het\tAC_Hom=%AC_Hom\tVQSLOD=%VQSLOD\n' test.annotated.vcf.gz | head
AF=0.0001652    AC_Het=20       AC_Hom=0        VQSLOD=1.88
AF=8.853e-05    AC_Het=4        AC_Hom=0        VQSLOD=-0.4803
AF=0.0001894    AC_Het=23       AC_Hom=0        VQSLOD=4.48
AF=3.295e-05    AC_Het=4        AC_Hom=0        VQSLOD=1.54
AF=8.236e-06    AC_Het=1        AC_Hom=0        VQSLOD=-1.197
AF=1.647e-05    AC_Het=2        AC_Hom=0        VQSLOD=1.81
AF=8.237e-06    AC_Het=1        AC_Hom=0        VQSLOD=1.21
AF=0.014        AC_Het=675      AC_Hom=358      VQSLOD=86.21
AF=8.277e-06    AC_Het=1        AC_Hom=0        VQSLOD=-1.574
AF=0.001051,1.655e-05,8.277e-06 AC_Het=119,2,1,0,0,0    AC_Hom=1,0,0    VQSLOD=4.54
~~~~

