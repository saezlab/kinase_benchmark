if(exists("snakemake")){
  ppsp_file <- snakemake@input$ppsp
  dataset <- snakemake@input$file_dataset
  output_file <- snakemake@output$tsv
}else{
  ppsp_file <- "results/prior/phosphositeplus.tsv"
  output_file <- "results/hernandez/prior/phosphositeplus.tsv"
  dataset <- "results/hernandez/processed_data/benchmark_data.csv"
}

## Libraries ---------------------------
library(biomaRt)
library(tidyverse)
library(org.Hs.eg.db)

## Construct kinase-substrate interaction network ---------------------------
network <- read_tsv(ppsp_file, col_types = cols())

## Prepare hernandez ---------------------------
hernandez_df <- read_csv(dataset, col_types = cols())
identifiers <- hernandez_df %>%
  dplyr::select(ID) %>%
  dplyr::mutate(gene_name = map_chr(str_split(ID, "\\|"), 1))  %>%
  dplyr::mutate(position = map_chr(str_split(ID, "\\|"), 2)) %>%
  dplyr::mutate(ensg = map_chr(str_split(ID, "\\|"), 3)) %>%
  dplyr::mutate(ensp = map_chr(str_split(ID, "\\|"), 4)) %>%
  dplyr::mutate(target_site = paste0(gene_name, "_", position)) %>%
  dplyr::rename("site" = ID)

## Merge with network ---------------------------
target_df <- identifiers

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

## Save processed phosphositeplus
write_tsv(prior_df, output_file)
