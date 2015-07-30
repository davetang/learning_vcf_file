#!/usr/bin/env Rscript

if("VariantAnnotation" %in% rownames(installed.packages()) == FALSE){
  stop('Please install the VariantAnnotation package from Biocondcutor')
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

info_df <- as.data.frame(info(header(vcf)))
for (i in rownames(info_df)){
  my_command <- paste('info(vcf)$', i, sep = '')
  x <- 'x'
  assign(x, eval(parse(text=my_command)))
  if (info_df[i, 2] == 'Integer'){
    x <- x[!is.na(x)]
    if (length(x)==0){
      next
    }
    barplot(table(unlist(x)), main=i, las=2, xlab='Count')
  } else if (info_df[i, 2] == 'Float'){
    x <- x[!is.na(x)]
    if (length(x)==0){
      next
    }
    hist(unlist(x), main=i, breaks=50, xlab='Breaks')
  }
}

dev.off()
quit()