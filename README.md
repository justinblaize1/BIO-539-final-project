# BIO-539-final-project
Git hub repository for BIO 539 final project

To perform big gene analysis first we start biomaRt search in R. Script is detailed below.

#Load library
library(biomaRt)

Warning: package 'biomaRt' was built under R version 3.6.3
#List databases, click on myMarts in Environment to view options
myMarts <- listMarts()

#Select database, for assignment CHANGE to Ensembl Variation
ensembl=useMart("ENSEMBL_MART_ENSEMBL")

#List datasets, click on myDatasets in Environment to view options
myDatasets<-listDatasets(ensembl)

#Select dataset, for assignment CHANGE to short Human variants
ensembl=useMart("ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl")

#List filters, click on myFilters in Environment to view options
myFilters <- listFilters(ensembl)

#Make a vector for your filters
filter1 <- 'ensembl_gene_id'

#Use/read in differential expression dataset
values1 <- read.csv("Differential_expression_analysis_table_B_vs_E.csv")

#List attributes, click on myAttributes in Environment to view options
myAttributes <- listAttributes(ensembl)

#Make a variable for your attributes
att1 <- c('ensembl_gene_id','chromosome_name','start_position','end_position',
          'transcript_start','transcript_end','transcript_length')

searchResults <-getBM(att1,
                      filters= filter1,
                      values= values1, mart=ensembl)

#Filter for max transcript length
searchResults <- searchResults %>% 
  group_by(ensembl_gene_id) %>% 
  filter(transcript_length == max(transcript_length)) %>% ungroup()

#Change name of ensemble gene id column to match the biomart file you just generated
colnames(values1) [1] <- "ensembl_gene_id"

#Merge two datasets using merge function
B_vs_E_comp <- merge(values1, searchResults, by = "ensembl_gene_id", all.x = TRUE)

#Create a column for gene length
B_vs_E_comp <- B_vs_E_comp %>% mutate(c=end_position-start_position)


#Rename column 
colnames(B_vs_E_comp) [20] <- "gene_length"


After performing biomaRt search in R. I filtered for big genes and made volcano plots as seen below.

## For big genes

library(tidyverse)
library(RColorBrewer)
library(ggrepel)

A_vs_D_1000 <- A_vs_D_comp %>% filter(gene_length > 220649)

A_vs_D_1000$diffexpressed <- 'no'
A_vs_D_1000$diffexpressed[A_vs_D_1000$log2FoldChange > 1 & A_vs_D_1000$padj < 0.05] <- 'UP'
A_vs_D_1000$diffexpressed[A_vs_D_1000$log2FoldChange < -1 & A_vs_D_1000$padj < 0.05] <- 'DOWN'

top5genes_1000 <- head(A_vs_D_1000[order(A_vs_D_1000$padj), 'Gene.name'], 5)

A_vs_D_1000$delabel <- ifelse(A_vs_D_1000$Gene.name %in% top5genes_1000, A_vs_D_1000$Gene.name, NA)
  
theme_set(theme_classic(base_size = 12) + theme(axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'), axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'), plot.title = element_text(hjust = -0.2), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))) 

ggplot(data = A_vs_D_1000, aes(x= log2FoldChange, y= -log10(padj), col = diffexpressed, label = delabel)) + geom_point(size = 2) + scale_color_manual(values = c("green", "gray", "red"), labels = c("Downregulated (60)", "Not Significant (831)", "Upregulated (109)")) + coord_cartesian(ylim = c(0, 200), xlim = c(-10,10)) + labs(color = "Differentially Expressed", x = expression("log"[2]*"FoldChange"), y = expression("-log"[10]*"padj")) + ggtitle('FA-D2 vs. FA-D2 + FANCD2 1000 largest genes') + geom_text_repel(max.overlaps = Inf, fontface = "italic")

## For small genes

A_vs_D_low_1000 <- A_vs_D_comp %>% filter(gene_length < 1999)

A_vs_D_low_1000$diffexpressed <- 'no'
A_vs_D_low_1000$diffexpressed[A_vs_D_low_1000$log2FoldChange > 1 & A_vs_D_low_1000$padj < 0.05] <- 'UP'
A_vs_D_low_1000$diffexpressed[A_vs_D_low_1000$log2FoldChange < -1 & A_vs_D_low_1000$padj < 0.05] <- 'DOWN'

top5genes_low_1000 <- head(A_vs_D_low_1000[order(A_vs_D_low_1000$padj), 'Gene.name'], 5)

A_vs_D_low_1000$delabel <- ifelse(A_vs_D_low_1000$Gene.name %in% top5genes_low_1000, A_vs_D_low_1000$Gene.name, NA)
  
theme_set(theme_classic(base_size = 12) + theme(axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'), axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'), plot.title = element_text(hjust = -0.2), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))) 

ggplot(data = A_vs_D_low_1000, aes(x= log2FoldChange, y= -log10(padj), col = diffexpressed, label = delabel)) + geom_point(size = 2) + scale_color_manual(values = c("green", "gray", "red"), labels = c("Downregulated (58)", "Not Significant (905)", "Upregulated (38)")) + coord_cartesian(ylim = c(0, 200), xlim = c(-10,10)) + labs(color = "Differentially Expressed", x = expression("log"[2]*"FoldChange"), y = expression("-log"[10]*"padj")) + ggtitle('FA-D2 vs. FA-D2 + FANCD2 1000 smallest genes') + geom_text_repel(max.overlaps = Inf, fontface = "italic")


## For big genes under replication stress




