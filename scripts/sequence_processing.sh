#!/bin/bash

################### Generate Folder Hierarchy #######################
cwd=$(pwd)

mkdir ${cwd}/PRJNA331054_analysis
mkdir ${cwd}/PRJNA331054_analysis/sequence_data
mkdir ${cwd}/PRJNA331054_analysis/sequence_data/paired_fastq
mkdir ${cwd}/PRJNA331054_analysis/sequence_data/overlapped_fastq
mkdir ${cwd}/PRJNA331054_analysis/sequence_data/filtered_fasta
mkdir ${cwd}/PRJNA331054_analysis/sequence_data/combined_fasta
mkdir ${cwd}/PRJNA331054_analysis/sequence_data/dereplicated_combined_fasta
mkdir ${cwd}/PRJNA331054_analysis/sequence_data/OTUs
mkdir ${cwd}/PRJNA331054_analysis/sequence_data/Taxonomy
mkdir ${cwd}/PRJNA331054_analysis/sequence_data/OTU_counts
## make directory for R analysis
mkdir ${cwd}/PRJNA331054_analysis/data_processing

base=${cwd}"/PRJNA331054_analysis"
raw=${base}"/sequence_data"
paired=${raw}"/paired_fastq"
overlap=${raw}"/overlapped_fastq"
filtered=${raw}"/filtered_fasta"
combo=${raw}"/combined_fasta"
derep=${raw}"/dereplicated_combined_fasta"
otu=${raw}"/OTUs"
tax=${raw}"/Taxonomy"
counts=${raw}"/OTU_counts"

cd ${base}

###################### Download fastq files ##########################
## requires NCBI SRA Toolkit
## Download here: https://github.com/ncbi/sra-tools/wiki/Downloads
## Manual Here: http://ncbi.github.io/sra-tools/
## SRA Toolkit version used here: 2.7.0
## Must add SRA Toolkit/bin folder to path
#######################################################################

## Notes: Unprocessed reads from this study are archived
## in the NCBI Sequence Read Archive (SRA). The SRA identifiers are:
## SRR3948168 --> DNA07 High Fe A
## SRR3948169 --> DNA08 High Fe B
## SRR3948170 --> DNA09 High Fe C
## SRR3948171 --> DNA10 Low Fe A
## SRR3948172 --> DNA11 Low Fe B
## SRR3948173 --> DNA12 Low Fe C

fastq-dump -I --split-files SRR3948168
mv SRR3948168_1.fastq ${paired}/DNA07_1.fastq
mv SRR3948168_2.fastq ${paired}/DNA07_2.fastq

fastq-dump -I --split-files SRR3948169
mv SRR3948169_1.fastq ${paired}/DNA08_1.fastq
mv SRR3948169_2.fastq ${paired}/DNA08_2.fastq

fastq-dump -I --split-files SRR3948170
mv SRR3948170_1.fastq ${paired}/DNA09_1.fastq
mv SRR3948170_2.fastq ${paired}/DNA09_2.fastq

fastq-dump -I --split-files SRR3948171
mv SRR3948171_1.fastq ${paired}/DNA10_1.fastq
mv SRR3948171_2.fastq ${paired}/DNA10_2.fastq

fastq-dump -I --split-files SRR3948172
mv SRR3948172_1.fastq ${paired}/DNA11_1.fastq
mv SRR3948172_2.fastq ${paired}/DNA11_2.fastq

fastq-dump -I --split-files SRR3948173
mv SRR3948173_1.fastq ${paired}/DNA12_1.fastq
mv SRR3948173_2.fastq ${paired}/DNA12_2.fastq

#################### Read quality control ####################
## Requires usearch
## Download here: http://www.drive5.com/usearch/download.html
## Manual Here: http://www.drive5.com/usearch/manual/
## usearch version used here: 8.1.1861
## Must add usearch6.1.1861... executable to path
##############################################################

## NOTES: Edgar says that for paired reads length truncation is probably
## unnecessary. Filter out reads with more than 1.0 total expected errors 
## for all bases in the read. There is more information about the expected
## error prediction in USEARCH here and here. This is more conservative threshold 
## than in many other comparable OTU processing pipelines. Edgar recommends using
## the default threshold of 1.0. This command also relabels each sequence 
## as Filt1, Filt2, Filt3…

declare -a arr=("DNA07" "DNA08" "DNA09" "DNA10" "DNA11" "DNA12")

# Merge paired end reads, relabel each read as FILENAME.1,2,3…
for i in "${arr[@]}"; do usearch8.1.1861_i86linux32 -fastq_mergepairs ${paired}/"$i"_1.fastq -reverse ${paired}/"$i"_2.fastq -relabel @ -fastqout ${overlap}/"$i"_overlapped.fastq; done

# Filter out reads with more than 1.0 total expected errors for all bases in the read
for i in "${arr[@]}"; do usearch8.1.1861_i86linux32 -fastq_filter ${overlap}/"$i"_overlapped.fastq -fastq_maxee 1.0 -relabel Filt -fastaout ${filtered}/"$i"_filtered.fasta; done

# Combine individual samples into single file
cat ${filtered}/*.fasta > ${combo}/combined_overlapped_samples.fasta

# Read dereplication
## NOTES: The input sequences to ’-cluster_otus’ must be a set of unique sequences
## sorted in order of decreasing abundance with size annotations in the labels.
## The ’-derep_fulllength’ command finds the unique sequences and adds size 
## annotaitions (ie roughly how many identical sequences there were 
## for each uniq sequence)

usearch8.1.1861_i86linux32 -derep_fulllength ${combo}/combined_overlapped_samples.fasta -relabel Uniq -sizeout -fastaout ${derep}/combined_overlapped_samples_derep.fasta

#################### OTU clustering ###########################
## Requires usearch
## Download here: http://www.drive5.com/usearch/download.html
## Manual Here: http://www.drive5.com/usearch/manual/
## usearch version used here: 8.1.1861
## Must add usearch6.1.1861... executable to path
################################################################

## Notes: Additionally, in USEARCH v8.1.1861 there is an additional
## chimera filtering step included in the -clusterotus command. Edgar
## states that it may be possible to find extra chimera sequences by
## running -uchimeref to remove predicted chimeric sequences using a
## certified reference database. In this case we use the UCHIME gold
## standard reference. It seems that the extra chimera filtering step
## is mostly unecessary, so we won’t be using it in this tutorial.
## Check the older pipeline for the commands for uchime.

usearch8.1.1861_i86linux32 -cluster_otus ${derep}/combined_overlapped_samples_derep.fasta -minsize 2 -relabel OTU -otus ${otu}/combined_overlapped_samples_derep_OTUs.fasta

#################### Map Reads to OTUs #########################
## Requires usearch
## Download here: http://www.drive5.com/usearch/download.html
## Manual Here: http://www.drive5.com/usearch/manual/
## usearch version used here: 8.1.1861
## Must add usearch6.1.1861... executable to path
################################################################

for i in "${arr[@]}"; do usearch8.1.1861_i86linux32 -usearch_global ${overlap}/"$i"_overlapped.fastq -db ${otu}/combined_overlapped_samples_derep_OTUs.fasta -strand plus -id 0.97 -otutabout ${counts}/"$i"_otutab.txt; done

#################### Taxonomy Assignment #######################
## Requires the RDP classifier
## Download precompiled version here: https://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/
## Manual Here: https://github.com/rdpstaff/classifier
## RDP classifier version used here: 2.12
## Must add classifier Jar file to path
################################################################

classifier.jar classify -c 0.8 -f fixrank -o ${tax}/combined_overlapped_samples_derep_OTUs_classified.txt -h ${tax}/combined_overlapped_samples_derep_OTUs_hierarchy.txt ${otu}/combined_overlapped_samples_derep_OTUs.fasta


############# Clean up files for inport to R ################

## Get rid of annoying space in output file
for i in "${arr[@]}"; do sed -e "s/ /_/g" ${counts}/"$i"_otutab.txt > ${counts}/"$i"_format.txt; done

## join the first two files by the first OTU column
## --header means that input files have header
## -a1 -a2 means to keep all fields not present in the other file
## -o auto means to create as many entries as in the first line (filles missing entries)
join --header -a1 -a2 -e 0 -o auto <(sort -k1,1 ${counts}/DNA07_format.txt) <(sort -k1,1 ${counts}/DNA08_format.txt) > ${counts}/tmp.tmp

for f in DNA09_format.txt DNA10_format.txt DNA11_format.txt DNA12_format.txt 
do                              
    join --header -a1 -a2 -e 0 -o auto <(sort -k1,1 ${counts}/tmp.tmp) <(sort -k1,1 ${counts}/$f) > tmpf           
    mv tmpf tmp.tmp                  
done

## ensures the header is at the first line
sort -k1,1r ${counts}/tmp.tmp > ${counts}/tmp1.tmp
sed -e "s/#//g" ${counts}/tmp1.tmp > ${counts}/all_samples_otu_counts.txt

## cleanup
rm ${counts}/tmp.tmp 
rm ${counts}/tmp1.tmp 
rm ${counts}/*_format.txt

## copy to R script processing directory
cp ${counts}/all_samples_otu_counts.txt ${base}/data_processing


## Format the Taxonomy file to a way that phyloseq likes it...
cut -f1,3,6,9,12,15,18 ${tax}/combined_overlapped_samples_derep_OTUs_classified.txt > ${tax}/tmp

## remove quotations
sed -e "s/\"//g" ${tax}/tmp > ${tax}/tmp1

## add header line
sed '1 i\OTU"\t"Domain"\t"Phylum"\t"Class"\t"Order"\t"Family"\t"Genus' ${tax}/tmp1 > ${tax}/phyloseq_tax_table.txt

## copy to R script processing directory
cp ${tax}/phyloseq_tax_table.txt ${base}/data_processing