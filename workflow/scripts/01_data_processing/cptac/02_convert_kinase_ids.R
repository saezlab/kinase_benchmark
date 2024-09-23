if(exists("snakemake")){
  prior_files <- snakemake@input$prior_files
  output_file <- snakemake@output$tsv
}else{
  prior_files <- list.files("results/00_prior", pattern = "tsv", recursive = T, full.names = T)
  output_file <- "resources/kinase_mapping.tsv"
}

## Libraries ---------------------------
library(tidyverse)
library(biomaRt)

## Load kinases in prior ---------------------------
prior <- map_dfr(prior_files, read_tsv)
kinases <- unique(prior$source)

## Convert Ensemble gene IDs to Gene names
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

res <- getBM(attributes = c('external_gene_name',
                            'ensembl_gene_id',
                            'uniprotswissprot'),
             values = kinases,
             mart = mart)

kin_map <- res %>%
  dplyr::filter(external_gene_name %in% kinases) %>%
  distinct() %>%
  group_by(external_gene_name, ensembl_gene_id) %>%
  filter(!(n() > 1 & uniprotswissprot == "")) %>%
  ungroup()

## Save kinase ID map ---------------------------
write_tsv(kin_map, output_file)
