#'

## Snakemake ---------------------------
if(exists("snakemake")){
  kinhub_file <- snakemake@input$kinhub
  quickgo_file <- snakemake@input$go
  output_file <- snakemake@output$filter
}else{
  kinhub_file <- "data/misc/kinase_list_kinhub.txt"
  quickgo_file <- "data/misc/QuickGO_kinase_20231214.tsv"
  output_file <- "data/misc/kinase_list.tsv"
}

## Libraries ---------------------------
library(tidyverse)
library(biomaRt)

## Load data ---------------------------
kinhub <- read_tsv(kinhub_file, col_types = cols())
quickgo <- read_tsv(quickgo_file, col_types = cols())

kinases <- unique(c(kinhub$`HGNC Name`, quickgo$SYMBOL))
kinases <- kinases[!is.na(kinases)]

# Add uniprot id of kinases to identify autophosphorylation
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_name <- getBM(attributes = c('external_gene_name'),
                   values = kinases,
                   mart = mart)

kinases <- kinases[kinases %in% gene_name$external_gene_name]

write_csv(data.frame(kinase = kinases), output_file)
