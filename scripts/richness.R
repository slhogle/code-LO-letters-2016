#!/usr/bin/Rscript

################################## README #########################################
## Before running this script you need to copy the file "all_samples_metadata.txt" 
## to the data_processing directory. This file can be found at 
## https://figshare.com/s/fdd50024f3bb2094f70b

## I would further recommend running this script using R studio

library(phyloseq)
library(ggplot2)

## set directory for R analysis. 
## Execute script from same location as sequence_processing.sh script
setwd(paste(getwd(),"PRJNA331054_analysis/data_processing",sep="/"))

## read cce12 phyloseq object
cce12 <- readRDS("cce12.rds")

## Richness/alpha diversity of data partitioned by Fe_Treatment
## Here we’ll make some richness plots of our data separated by high iron and low 
## iron treatments. We use raw data instead of the normalized data here because many 
## richness estimates are based on singletons and doubletons. We need to leave them 
## in order to get meaningful estimates.

## Get numerical values for richness estimates and write output as a table
cce12_st <- subset_taxa(cce12, Phylum!="Cyanobacteria/Chloroplast")
cce12_st <- subset_taxa(cce12_st, Domain=="Bacteria")

richness <- estimate_richness(cce12_st)
write.table(richness, file = "sample_richness_estimates.tsv", sep = "\t", row.names = T)

p = plot_richness(cce12_st, x = "Fe_Treatment", color = "Sample_ID")
p + geom_point(data = p$data, aes(x = Fe_Treatment, y = value, color = Sample_ID), size = 4) + ggtitle("Richness indices")
ggsave("richness.pdf")

## Right now, we actually only care about richness plots of Chao1 and “observed” diversity so we’ll just select those
## This is part A of supplemental figure 5
psub = plot_richness(cce12_st, x = "Fe_Treatment", color = "Sample_ID", measures = c("Observed", "Chao1"))
psub + geom_point(data = psub$data, aes(x = Fe_Treatment, y = value, color = Sample_ID), size = 4) + ggtitle("Subset richness indices")
ggsave("supp_fig5_A.pdf")
