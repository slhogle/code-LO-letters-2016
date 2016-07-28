#!/usr/bin/Rscript

################################## README #########################################
## Before running this script you need to copy the file "all_samples_metadata.txt" 
## to the data_processing directory. This file can be found at 
## https://figshare.com/s/fdd50024f3bb2094f70b

## I would further recommend running this script using R studio

library(phyloseq)
library(ggplot2)
library(reshape2)

## set directory for R analysis. 
## Execute script from same location as sequence_processing.sh script
setwd(paste(getwd(),"PRJNA331054_analysis/data_processing",sep="/"))

## read cce12 phyloseq object
cce12norm <- readRDS("cce12norm.rds")

###################################################################################
## Here we'll calculate some average read counts for high and Low Fe samples
## that we'll be using to make relative abundance bar plots

## Subset datasets for Hi and Low Fe
HiFe <- subset_samples(cce12norm, Fe_Treatment=="Hi_Fe")
LoFe <- subset_samples(cce12norm, Fe_Treatment=="Lo_Fe")

## Calculate mean number of normalized reads in High and Low Fe samples
HiFe_mean <- sum(sample_sums(HiFe))/3
LoFe_mean <- sum(sample_sums(LoFe))/3

## Calculate mean number of normalized reads per OTU in each Hi Fe sample
HiFe_OTU_mean <- taxa_sums(HiFe)/3
LoFe_OTU_mean <- taxa_sums(LoFe)/3

## Merge Hi and Low Fe OTU means into single table
Fesums <- merge(HiFe_OTU_mean, LoFe_OTU_mean, by="row.names", all=TRUE)
colnames(Fesums) <- c("OTU", "hiFe", "loFe")
rownames(Fesums) <- Fesums[,1]
Fesums <- subset(Fesums, select = -OTU)

###################################################################################

## Now let’s make some bar plots to get a fuller picture as to how the top OTUs are 
## distributed throughout our Fe treatments. Here we’ll select the top 20 most 
## abundant OTUs from the normalized dataset “cce12norm” and partition them by Family

topOTUs = names(sort(taxa_sums(cce12norm), TRUE)[1:20])
cce12norm20 = prune_taxa(topOTUs, cce12norm)

## Transform to percent relative abundance. Reads are now in percentage of 
## OTU X relative to the mean number of normalized OTUs across all samples
## The actual function just divides by mean of HiFe_mean and LoFe_mean
## which are both practically the same value thus creating a mean normalized 
## value across all six samples. Since the samples were already normalized this
## doesnt really change anything, otherwise OTUS from each sample should be normalized 
## to each individual sample.
cce12norm20.rel <- transform_sample_counts(cce12norm20 , function(x) (x / (mean(HiFe_mean, LoFe_mean)))*100 )

## Combine OTUs at rank family - makes for prettier plotting
glommed <- tax_glom(cce12norm20.rel, taxrank="Family")

## Supplemental Fig. 4: Bar plot showing the relative 
## abundance of the 20 most abundant OTUs binned by taxonomic family
plot_bar(glommed, "Sample_ID", fill = "Family", facet_grid = ~Family) + 
  theme_bw() +
  xlab("Sample") +
  ylab("% Relative Abundance") +
  ggtitle("relative abundance of the 20 most abundant OTUs binned by taxonomic family")
ggsave("supp_fig4.pdf", width = 20, height = 10)

###################################################################################
## Now we'll make a plot for OTUs that are disproportionately abundant in Hi Fe
## samples relative to Low Fe samples

## Select only OTUs from Fesums dataframe that are in greater abundance than 
## 1% of the mean number of normalized reads across all Hi Fe samples and that
## are greater than 1.5x the mean number of normalized reads across Low Fe
## samples
selotus <- row.names(subset(Fesums, hiFe>(0.01*HiFe_mean) & hiFe > 1.5*loFe))

## Reduce the phyloseq object down to the selected OTUS
cce12norm.select <- prune_taxa(selotus, cce12norm)

## Transform to percent relative abundance. Reads are now in percentage of 
## OTU X relative to the mean number of normalized OTUs across all samples
## The actual function just divides by mean of HiFe_mean and LoFe_mean
## which are both practically the same value thus creating a mean normalized 
## value across all six samples. Since the samples were already normalized this
## doesnt really change anything, otherwise OTUS from each sample should be normalized 
## to each individual sample.
cce12norm.select.rel <- transform_sample_counts(cce12norm.select , function(x) (x / (mean(HiFe_mean, LoFe_mean)))*100 )

## Main text Fig. 3: Bar plot of OTUs with a mean relative abundance greater 
## than 1% of all reads in each replicate and where the High Fe mean is at 
## least 1.5 times greater than the Low Fe mean. OTUs are colored by taxonomic order. 
plot_bar(cce12norm.select.rel, "Sample_ID", fill = "Order", facet_grid = ~OTU) + 
  theme_bw() +
  xlab("Sample") +
  ylab("% Relative Abundance") +
  ggtitle("relative abundance of the OTUs overrepresented in High Fe treatments")
ggsave("Figure3.pdf", width = 20, height = 10)

###################################################################################
## To make bar plots of all OTUs combined into a particular taxonomic rank you 
## need to use the tax_glom function.

## Combine OTUs at the taxonomic rank of Class. Results in 23 OTUs
glommed <- tax_glom(cce12norm, taxrank="Class")
                    
## Create variable “OTUS_remove” that contains the names of the 19 least 
## abundant OTUs in the glommed dataset
OTUS_remove <- names(sort(taxa_sums(glommed), FALSE)[1:18])

## Remove the 19 least abundant Class “OTUs” from the glommed dataset
glommed.red <- merge_taxa(glommed, OTUS_remove, 1)

## When you make the bar plot using the phyloseq wrapper “plot_bar” it puts the 
## components not in the same order for each sample. I don’t actually understand 
## how this wrapper is implemented in phyloseq so I'll just set it up with 
## vanilla ggplot
OTUs.class <- data.frame(otu_table(glommed.red), stringsAsFactors=FALSE)
tax.class <- data.frame(tax_table(glommed.red), stringsAsFactors=FALSE)

## Merge these data frames into one
merged.class.plot <- merge(OTUs.class, tax.class, by="row.names", all=TRUE)

## We’ll now get rid of unecessary columns
merged.class.plot.red <- subset(merged.class.plot, select=-c(Row.names, Domain, Phylum, Order, Family, Genus))
merged.class.plot.red[4,7] <- "Other"
rownames(merged.class.plot.red) <- merged.class.plot.red[,7]
merged.class.plot.red$Class <- NULL

## Change the dataframe so it's in relative abundance
merged.class.plot.red.rel <- merged.class.plot.red/colSums(merged.class.plot.red)
colnames(merged.class.plot.red.rel)[7] <- "Class"

## make a column that is the same as the rownames (for reshaping)
merged.class.plot.red.rel[,7] <- rownames(merged.class.plot.red.rel)

## Now we’ll put df into long format and select id.vars as “Class” using the 
## reshape package
merged.class.plot.red.rel.m <- melt(merged.class.plot.red.rel, id.vars="Class")

## Finally now we’ll make the bar plot with the correct ordering of the 
## different classes, which looks a lot nicer.
ggplot(merged.class.plot.red.rel.m, aes(variable, value, fill=Class)) + 
  geom_bar(stat = "identity") + 
  xlab("Sample") +
  ylab("Relative Abundance") +
  ggtitle("Relative Abundance at Class Level")
ggsave("supp_fig5_B.pdf")

