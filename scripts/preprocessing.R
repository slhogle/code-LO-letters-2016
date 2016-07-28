#!/usr/bin/Rscript

################################## README #########################################
## Before running this script you need to copy the file "all_samples_metadata.txt" 
## to the data_processing directory. This file can be found at 
## https://figshare.com/s/fdd50024f3bb2094f70b

## I would further recommend running this script using R studio

library(phyloseq)
library(ggplot2)

packageVersion("phyloseq")
packageVersion("ggplot2")

## set directory for R analysis. 
## Execute script from same location as sequence_processing.sh script
setwd(paste(getwd(),"PRJNA331054_analysis/data_processing",sep="/"))

## Load the OTU matrix generated in the sequence_processing.sh script
otu.mat <- as.matrix(read.table("all_samples_otu_counts.txt", row.names=1, sep=" ", header=TRUE))

## Load and create a taxonomy matrix
tax.mat <- as.matrix(read.table("phyloseq_tax_table.txt", row.names=1, sep="\t", header=TRUE))

## Load and create a metadata data frame
meta.df <- as.data.frame(read.table("all_samples_metadata.txt", row.names=1, sep="\t", header=TRUE))

## Generate phyloseq object
## Here we’ll call some of phyloseq’s functions to genearate an object called cce12 that contains metadata, OTU, and taxonomy info.

OTU <- otu_table(otu.mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax.mat)
SAMPDAT <- sample_data(meta.df)
cce12 <- phyloseq(OTU, TAX, SAMPDAT)

## Save cce12 phyloseq object for later use
saveRDS(cce12, "cce12.rds")

## Generate OTU sums across all samples and sort
OTUsums <- data.frame(nreads = sort(taxa_sums(cce12), TRUE), sorted = 1:ntaxa(cce12))

## Now calculate the number of reads in each sample and make a table
samplesums <- data.frame(read.count = sort(sample_sums(cce12), TRUE))

## Most of the samples are pretty evenly sequenced with the exception of DNA08. 
# Now plot an OTU rank abundance curve using ggplot. 

## Generate rank abundance plot
ggplot(OTUsums, aes(x = sorted, y = nreads)) + 
  geom_line(stat = "identity", size = 1) + 
  scale_y_log10() + 
  xlab("OTU rank abundance") +
  ylab("Total number of reads") +
  ggtitle("Rank abundance plot")
ggsave("rank_abund.pdf")

## We can see that there are many OTUs with a very low abundance

## Preprocessing data
## Here we’ll do some preprocessing to remove spurious and
## nonabundant OTUs that we don’t want to affect downstream analysis
## Here we use the subset taxa function to select only the Domain 
## “Bacteria” while excluding “Cyanobacteria/chloroplast.” 
## cce12_st is a new phyloseq object
cce12_st <- subset_taxa(cce12, Phylum!="Cyanobacteria/Chloroplast")
cce12_st <- subset_taxa(cce12_st, Domain=="Bacteria")

## Now we will transform the reads for each OTU in each 
## sample to relative abundance - ie. reads for OTUx / total reads in sample x
cce12st.rel.abun <- transform_sample_counts(cce12_st, function(x) x / sum(x))
cce12st.rel.abun

## The cce12st.rel.abun has 832 OTUs

## Now we will remove OTUs that have a mean percentage abundance of less than 0.005%. 
cce12st.rel.abun.filter <- filter_taxa(cce12st.rel.abun, function(x) mean(x) > 5e-5, TRUE)
cce12st.rel.abun.filter

## The resulting phyloseq object has 304 OTUs

## Collect OTUs from cce12.rel.abun.filter and use them to filter the raw OTUs from cce12
cce12.raw.filt <- prune_taxa(taxa_names(cce12st.rel.abun.filter), cce12)

## Remove taxa not seen more than 3 times in at least 20% of the samples. 
## This protects against OTUs with small mean & trivially large covariance
cce12.raw.filt.2 <- filter_taxa(cce12.raw.filt, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

## Save this data object for use with differential abundance testing
saveRDS(cce12.raw.filt.2, "cce12.raw.filt.2.rds")

## Generate OTU sums across all samples in the new filtered set and sort
OTUsums1 <- data.frame(nreads = sort(taxa_sums(cce12.raw.filt.2), TRUE), sorted = 1:ntaxa(cce12.raw.filt.2))

## Generate new rank abundance plot based on filtered OTUs
ggplot(OTUsums1, aes(x = sorted, y = nreads)) + 
  geom_line(stat = "identity", size = 1) + 
  scale_y_log10() + 
  xlab("OTU rank abundance") +
  ylab("Total number of reads") +
  ggtitle("Rank abundance plot after OTU filtering")
ggsave("rank_abund_filter.pdf")

## Again make a summary table of read counts for each library
## Some very annoying formatting to make a nice table

names <- c("DNA07", "DNA08", "DNA09", "DNA10", "DNA11", "DNA12")
cond <- c("High Fe A", "High Fe B", "High Fe C", "Low Fe A", "Low Fe B", "Low Fe C")
df.cond <- data.frame(names, cond)

d <- samplesums
names <- rownames(d)
rownames(d) <- NULL
samplesums1 <- cbind(names, d)

filt1 <- data.frame(read.count.postfilter = sort(sample_sums(cce12.raw.filt), TRUE))
e <- filt1
names <- rownames(e)
rownames(e) <- NULL
filt.1 <- cbind(names, e)

filt2 <- data.frame(read.count.postfilter.2 = sort(sample_sums(cce12.raw.filt.2), TRUE))
f <- filt2
names <- rownames(f)
rownames(f) <- NULL
filt.2 <- cbind(names, f)

df.list <- list(df.cond, samplesums1, filt.1, filt.2)
seq.stats <- Reduce(function(...) merge(..., all=T), list(df.cond, samplesums1, filt.1, filt.2))

## write a table consisting of sequencing stats
write.table(seq.stats, sep = "\t", row.names = FALSE, file = "supplemental_table_5.tsv")

## Now we’ll standardize abundances to the median sequencing depth of each sample
## First calculate the median for the sample sums in the filtered dataset
total.median <- median(sample_sums(cce12.raw.filt.2))

## Now create a function that standardized the relative 
## abundance of each OTU to the median sequencing depth
stdize.fun <- function(x, t=total.median) round(t * (x / sum(x)))

## Create the standardized result
cce12norm <- transform_sample_counts(cce12.raw.filt.2, stdize.fun)

## save the final formatted phyloseq object for later use
saveRDS(cce12norm, "cce12norm.rds")
