#!/usr/bin/env Rscript

if("SNPRelate" %in% rownames(installed.packages()) == FALSE){
  stop('Please install the SNPRelate package from Biocondcutor')
}

args <- commandArgs(TRUE)
if (length(args) != 1){
  stop('Please input a single file name')
}

infile <- args[1]

if (!file.exists(infile)){
  msg <- paste('File does not exist: ', infile, sep='')
  stop(msg)
}

library(SNPRelate)

my_gds <- gsub(pattern="vcf$", replacement="gds", x=infile, ignore.case=TRUE)
my_base <- gsub(pattern=".vcf$", replacement="", x=infile, ignore.case=TRUE)

snpgdsVCF2GDS(infile, my_gds, method="copy.num.of.ref")

genofile <- snpgdsOpen(my_gds)

snpgdsGDS2PED(genofile, my_base)

quit()
