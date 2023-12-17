if(exists("snakemake")){
  ppsp_file <- snakemake@input$ppsp
  file_datasets <- snakemake@input$file_dataset
  output_file <- snakemake@output$tsv
}else{
  ppsp_file <- "results/prior/phosphositeplus.tsv"
  output_file <- "results/cptac/prior/phosphositeplus.tsv"
  file_datasets <- list.files("data/CPTAC_phospho/final", full.names = T)
}

## Libraries ---------------------------
library(biomaRt)
library(tidyverse)
library(org.Hs.eg.db)

## Construct kinase-substrate interaction network ---------------------------
network <- read_tsv(ppsp_file, col_types = cols())

## Prepare CPTAC ---------------------------
# Map targets to pps in data
pps <- map_dfr(file_datasets, function(file){
  df <- readRDS(file)
  data.frame(site = rownames(df))
})

pps <- pps %>%
  distinct(.keep_all = TRUE)

identifiers <- map_chr(str_split(pps$site,  "\\|"), 1)
identifiers <- map_chr(str_split(identifiers,  "\\."), 1)

pps_df <- data.frame(site = pps,
                     ensembl_gene_id = identifiers,
                     position = map_chr(str_split(pps$site,  "\\|"), 3),
                     surrounding = map_chr(str_split(pps$site,  "\\|"), 4))

## Convert Ensemble gene IDs to Gene names
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

res <- getBM(attributes = c('ensembl_gene_id',
                            'external_gene_name'),
             values = identifiers,
             mart = mart)

target_df <- full_join(pps_df, res, by = "ensembl_gene_id")

## Merge with network ---------------------------
if(!any(is.na(network$sequence))){
  target_df <- target_df %>%
    mutate(target_site = paste0(external_gene_name, "_", surrounding))

  network <- network %>%
    mutate(target_site = paste0(target_protein, "_", sequence))

  mapped_prior <- left_join(network %>%
                            dplyr::select(source, target_site),
                          target_df, by = "target_site", relationship = "many-to-many")

  prior_df <- mapped_prior %>%
    mutate(target = case_when(
      !is.na(site) ~ site,
      is.na(site) ~ target_site
    )) %>%
    mutate(target = case_when(
      str_detect(pattern = source, string = target_site) ~ paste0(target, "|auto"), #mark autophosphorylation
      !str_detect(pattern = source, string = target_site) ~ target
    )) %>%
    dplyr::mutate(mor = 1) %>%
    dplyr::select(source, target, mor) %>%
    distinct()
} else {
  target_df <- target_df %>%
    mutate(target_site = paste0(external_gene_name, "_", position))

  network <- network %>%
    mutate(target_site = paste0(target_protein, "_", position))

  mapped_prior <- left_join(network %>%
                              dplyr::select(source, target_site),
                            target_df, by = "target_site", relationship = "many-to-many")

  prior_df <- mapped_prior %>%
    mutate(target = case_when(
      !is.na(site) ~ site,
      is.na(site) ~ target_site
    )) %>%
    mutate(target = case_when(
      str_detect(pattern = source, string = target_site) ~ paste0(target, "|auto"), #mark autophosphorylation
      !str_detect(pattern = source, string = target_site) ~ target
    )) %>%
    dplyr::mutate(mor = 1) %>%
    dplyr::select(source, target, mor) %>%
    distinct()
}


## Save processed phosphositeplus
write_tsv(prior_df, output_file)
