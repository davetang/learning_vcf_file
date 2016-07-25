#!/bin/bash

# wget https://www.cog-genomics.org/static/bin/plink160705/plink_linux_x86_64.zip
# unzip plink_linux_x86_64.zip

# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/browser/vcf_to_ped_converter/version_1.1/vcf_to_ped_convert.pl

vcf_to_ped_convert.pl -vcf ex2.vcf.gz -sample_panel_file ex2.panel -region 2 -base_format letter -population AUS

add_phenotype.pl 2.ped ex2.panel > blah
mv -f blah 2.ped

info_to_map.pl 2.info > 2.map

plink --file 2 --make-bed --out 2

plink --bfile 2 --assoc --out data --allow-no-sex

plink --bfile 2 --model --out data --allow-no-sex

