#'

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$raw
  kinhub_file <- snakemake@input$kinhub
  quickgo_file <- snakemake@input$go
  output_file <- snakemake@output$filter
}else{
  input_file <- "results/prior/raw/GSknown.tsv"
  kinhub_file <- "data/kinase_list_kinhub.txt"
  quickgo_file <- "data/QuickGO_kinase_20231214.tsv"
  output_file <- "results/prior/GSknown.tsv"
}

## Libraries ---------------------------
library(tidyverse)
library(biomaRt)

## Load data ---------------------------
net <- read_tsv(input_file, col_types = cols())

kinhub <- read_tsv(kinhub_file, col_types = cols())
quickgo <- read_tsv(quickgo_file, col_types = cols())

kinases <- unique(c(kinhub$`HGNC Name`, quickgo$SYMBOL))
kinases <- kinases[!is.na(kinases)]

# Add uniprot id of kinases to identify autophosphorylation
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

gene_name <- getBM(attributes = c('external_gene_name'),
                   values = kinases,
                   mart = mart)

kinases <- kinases[kinases %in% gene_name$external_gene_name]


## Check kinases not covered in data ---------------------------
net %>%
  filter(!source %in% kinases) %>%
  pull(source) %>%
  unique()

net_filtered <- net %>% filter(source %in% kinases)

write_tsv(net_filtered, output_file)
