#!/usr/bin/Rscript

################################## README #########################################
## Before running this script you need to copy the file "all_samples_metadata.txt" 
## to the data_processing directory. This file can be found at 
## https://figshare.com/s/fdd50024f3bb2094f70b

## I would further recommend running this script using R studio

library(phyloseq)
library(DESeq2)

## set directory for R analysis. 
## Execute script from same location as sequence_processing.sh script
setwd(paste(getwd(),"PRJNA331054_analysis/data_processing",sep="/"))

## It is very important to document and describe the exact data manipulation that 
## occurred before the differential testing step because it can be dangerous to alter
## or omit data without any clear documentation or reasoning. In our case we want to 
## perform differential testing on non transformed counts (ie not relative counts or
## or normalized in any way) but also we want to avoid OTUs that are trivially abundant
## or that have a trivial mean but a strong covariance
## The phyloseq object "cce12.raw.filt.2" has no count normalization but it has discarded OTUs that: 
##  1. have a mean percentage abundance of less than 0.005%.
##  2. are not seen more than 3 times in at least 20% of the samples.

## read cce12 phyloseq object.
cce12.raw.filt.2 <- readRDS("cce12.raw.filt.2.rds")

###################################################################################
## So there seems to be some controvery in the microbiome world about whether or 
## not to Rarefy data. Rarefying is basically a library size normalization technique 
## that randomly subsamples reads from each library/sample in order to make all libraries 
## contain the same number of reads as the smallest library. This is the default 
## normalization proceedure in QIIME as far as I can tell. A very recent paper 
## (https://peerj.com/preprints/1157/) states that Rarefying is a good option when they 
## are close to each other in read number (like ours) and you have a small sample 
## size (also like ours) then procedures utilizing generalized linear models (like DESeq2)
## to be more sensitive and powerful.

## Some researchers have gone as far as saying that Rarefying is “inexcusable”. I am not a
## statastician so I don't feel qualified to say what is or is not inaddmissable, but the paper
## McMurdie and Holmes is quite compelling. (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)
## It seems like the best option for us is to model microbiome counts with a Negative Binomial model
## as implemented in the popular package DESeq2 for detecting differential expression in RNAseq data. 
## Here we use DESeq2 to test whether any particular OTUs are significantly less or more abundant 
## in high iron versus low iron treatments.

## Import data with phyloseq, convert to DESeq2’s DESeqDataSet class
FeTreat.dds = phyloseq_to_deseq2(cce12.raw.filt.2, ~ Fe_Treatment)

## Calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(FeTreat.dds), 1, gm_mean)
FeTreat.dds = estimateSizeFactors(FeTreat.dds, geoMeans = geoMeans)

## Running DESeq2 on non transformed data
FeTreat.dds = DESeq(FeTreat.dds, fitType="local")
res = results(FeTreat.dds)
res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < 0.05 & res$baseMean > 100 ), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(cce12)[rownames(sigtab), ], "matrix"))

## Write the main text table 1
write.table(sigtab, file = "table_1.tsv", sep = "\t", row.names = T)
