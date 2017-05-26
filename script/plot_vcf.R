#!/usr/bin/env Rscript
#
# Usage: plot_vcf.R <infile.vcf>
#

if("VariantAnnotation" %in% rownames(installed.packages()) == FALSE){
  stop('Please install the VariantAnnotation package from Bioconductor')
}

args <- commandArgs(TRUE)
if (length(args) != 1){
  stop('Please input a single argument')
}

infile <- args[1]

if (!file.exists(infile)){
  msg <- paste('File does not exist: ', infile, sep='')
  stop(msg)
}

library(VariantAnnotation)

outfile <- gsub('vcf', 'pdf', infile)
genome <- 'hg19'
vcf <- readVcf(infile, genome)
pdf(outfile)

my_qual <- qual(vcf)
hist(my_qual, main = 'QUAL')
my_filt <- filt(vcf)
barplot(table(my_filt), main = "FILTER", las = 2)

info_df <- as.data.frame(info(header(vcf)))
for (i in rownames(info_df)){
  print(i)
  my_command <- paste('info(vcf)$', i, sep = '')
  x <- 'x'
  assign(x, eval(parse(text=my_command)))
  x <- x[!is.na(x)]
  if (length(x) == 0){
    next
  }
  if (info_df[i, 2] == 'Integer'){
    barplot(table(unlist(x)), main=i, las=2, xlab='Count')
  } else if (info_df[i, 2] == 'Float'){
    hist(unlist(x), main=i, breaks=50, xlab='Breaks')
  }
}

dev.off()
quit()
