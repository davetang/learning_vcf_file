VCF to PED
----------

Example of generating a PED (and BED) file from a VCF file, allowing one to perform analyses using [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/).

~~~~{.bash}
# download PLINK
wget https://www.cog-genomics.org/static/bin/plink160705/plink_linux_x86_64.zip
unzip plink_linux_x86_64.zip

# download Perl script to convert VCF to PED
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/browser/vcf_to_ped_converter/version_1.1/vcf_to_ped_convert.pl

# convert VCF to PED
vcf_to_ped_convert.pl -vcf ex2.vcf.gz -sample_panel_file ex2.panel -region 2 -base_format letter -population AUS

# add phenotype to PED file
add_phenotype.pl 2.ped ex2.panel > blah
mv -f blah 2.ped

# create map file
info_to_map.pl 2.info > 2.map

# create BED from PED
plink --file 2 --make-bed --out 2

# basic descriptive summary
plink --bfile 2 --assoc --out data --allow-no-sex

# single SNP tests of association
plink --bfile 2 --model --out data --allow-no-sex
~~~~

## Further reading

[Basic statistical analysis in genetic case-control studies](http://www.nature.com/nprot/journal/v6/n2/abs/nprot.2010.182.html)

