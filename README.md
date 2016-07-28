# Introduction
Scripts for reproducing analyses from Hogle et al. L&amp;O Letters. This project investigated the production of different classes of iron-binding ligands, iron concentrations, macronutrient concentrations, and phytoplankton and bacterioplankton assemblage composition in iron amended microcosm incubations conducted in oligotrophic waters collected off the Southern California Bight.
--------
# Data and Code
Scripts for reproducing analysis are found in the scripts folder. Metadata is located in the metadata folder. These files are needed to follow the analysis in the shell and R scripts.
--------
# Required Software
## NCBI SRA Toolkit v2.7.0
* [Download](https://github.com/ncbi/sra-tools/wiki/Downloads)  
* [Manual](http://ncbi.github.io/sra-tools/)  
* Must add SRA Toolkit/bin folder to path  

## USEARCH v8.1.1861
* [Download](http://www.drive5.com/usearch/download.html)  
* [Manual](http://www.drive5.com/usearch/manual/)  
* Must add usearch6.1.1861... executable to path  

## RDP classifier v2.12
* [Download](https://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/)  
* [Manual](https://github.com/rdpstaff/classifier)  
* Must add classifier Jar file to path  

## R v3.3.0
* [Download](https://www.r-project.org/)  
* [Manual](https://cran.r-project.org/manuals.html)  

## R Libraries
* phyloseq  
* ggplot2  
* reshape2  
* grid  
* vegan  
* DESeq2  
--------
# How to reproduce analysis...
After installing the above software and adding it to your path run the scripts in the following order

```shell
bash sequence_processing.sh
R preprocessing.R
R richness.R
R bar_plots.R
R diff_abundance.R
R ordination

All scripts are heavily documented so open them in a text editor (or R studio) to see what is really going on  
  
The sequence_processing.sh script should create a directory hierarchy for the project that is self explanatory. It will also connect remotely to the NCBI SRA to download the sequecing data. For the sake of completeness, you can find the links to the NCBI SRA here:  
  
* [Bioproject](http://www.ncbi.nlm.nih.gov/bioproject/PRJNA331054)  
* [DNA07](http://www.ncbi.nlm.nih.gov/sra/SRX1973167)  
* [DNA08](http://www.ncbi.nlm.nih.gov/sra/SRX1973168)  
* [DNA09](http://www.ncbi.nlm.nih.gov/sra/SRX1973169)  
* [DNA10](http://www.ncbi.nlm.nih.gov/sra/SRX1973170)  
* [DNA11](http://www.ncbi.nlm.nih.gov/sra/SRX1973171)  
* [DNA12](http://www.ncbi.nlm.nih.gov/sra/SRX1973172)  
