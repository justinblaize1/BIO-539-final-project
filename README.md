# BIO-539-final-project
Git hub repository for BIO 539 final project

To perform big gene analysis first we start biomaRt search.

#Load library
library(biomaRt)

## Warning: package 'biomaRt' was built under R version 3.6.3
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

