#!/usr/bin/Rscript

################################## README #########################################
## Before running this script you need to copy the file "all_samples_metadata.txt" 
## to the data_processing directory. This file can be found at 
## https://figshare.com/s/fdd50024f3bb2094f70b

## I would further recommend running this script using R studio

library(phyloseq)
library(ggplot2)
library(vegan)
library(grid) ## to load the arrows function for drawing arrows on plots

## set directory for R analysis. 
## Execute script from same location as sequence_processing.sh script
setwd(paste(getwd(),"PRJNA331054_analysis/data_processing",sep="/"))

## read cce12 phyloseq object.
cce12norm <- readRDS("cce12norm.rds")

## Ordination of data using phyloseq, vegan, ggplot and some other stuff…
## Finally we’ll make an ordination based on our community data and fit continuous environmental 
## variable to it as factors

## First define a function veganotu for extracting an OTU table from a phyloseq object 
## and coercing it into a form vegan likes

veganotu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

## Convert phyloseq metadata into vanilla R data frames
keepvariables = which(sapply(sample_data(cce12norm), is.numeric))
physeqsd = data.frame(sample_data(cce12norm))[keepvariables]
physeqsd = subset(physeqsd, select = -Fe_level)
treatment = data.frame(sample_data(cce12norm))["Fe_Treatment"]

## Using the normalized OTU matrix in the metaMDS function
set.seed(4019)
vare.mds <- metaMDS(veganotu(cce12norm), trace = TRUE, 
                    plot = FALSE, autotransform = TRUE, 
                    noshare = FALSE, wascores = FALSE, 
                    trymax = 1000)

## Use the score function to load sample sites into vectors. We also use plot_ordination
## function with the justDF option =TRUE in order to load the specific coordinates of 
## just the OTUs.

scrs.site <- as.data.frame(scores(vare.mds, display = "sites"))
scrs.site <- cbind(scrs.site, treatment = treatment)
scrs.spp <- plot_ordination(cce12norm, vare.mds, type = "species", justDF = T) # To get species coordinates

## Now we use the envfit function to fit vectors to our continuous environmental data
set.seed(152)
ef <- envfit(vare.mds, physeqsd, permu = 999)

## Use the score function again, but this time to load direction of environmental vectors. Environmental 
## vectors here are weighted by R2 value
scrs.vct <- as.data.frame(scores(ef, display = "vectors"))

## ef is a list of lists. to get pvals you need to selected the pvals list from the vectors list
vct.pvals <- ef[['vectors']]['pvals']

## combine both vectors and pvals
scrs.vct <- cbind(scrs.vct, vct.pvals)

## select only rows with a P value < 0.1
scrs.vct <- scrs.vct[scrs.vct$pvals < 0.1,  ]

## Make the plot
p <- ggplot(scrs.site, aes(x = NMDS1, y = NMDS2))
p + geom_point(data = scrs.site, aes(shape = Fe_Treatment), size = 4) + 
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = scrs.vct, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "black") +
  geom_text(data = scrs.vct, aes(x = NMDS1, y = NMDS2, label = rownames(scrs.vct)), size = 3) +
  geom_point(data = scrs.spp, mapping = aes(x = NMDS1, y = NMDS2, color = Class), size = 1.5)
