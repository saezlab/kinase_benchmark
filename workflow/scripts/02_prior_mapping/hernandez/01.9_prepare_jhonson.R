if(exists("snakemake")){
  ppsp_file <- snakemake@input$ppsp
  dataset <- snakemake@input$file_dataset
  output_file <- snakemake@output$tsv
}else{
  ppsp_file <- "data/prior/simplified_jhonson_with_psite.csv"
  output_file <- "results/hernandez/prior/jhonson.tsv"
  dataset <- "results/hernandez/processed_data/benchmark_data.csv"
}

## Libraries ---------------------------
library(biomaRt)
library(tidyverse)
library(org.Hs.eg.db)

## Construct kinase-substrate interaction network ---------------------------
net <- read_csv(ppsp_file, col_types = cols())
net <- net %>%
  mutate(source_tmp = case_when(
    str_detect(source, "_") ~ target,
    str_detect(target, "_") ~ source
  )) %>%
  mutate(target_tmp = case_when(
    str_detect(source, "_") ~ source,
    str_detect(target, "_") ~ target
  )) %>%
  dplyr::select(source_tmp, target_tmp, source_protein, target_protein)

net <- net %>%
  dplyr::rename("target_site" = target_tmp) %>%
  dplyr::rename("uniprot_gn_id" = source_tmp) %>%
  mutate(autophosphorylation = case_when(
    source_protein == target_protein ~ "auto",
    !source_protein == target_protein ~ "no"
  ))

## Convert Ensemble gene IDs to Gene names
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

res <- getBM(attributes = c('uniprot_gn_id',
                            'external_gene_name'),
             values = unique(net$uniprot_gn_id),
             mart = mart)

net_df <- left_join(net, res, by = "uniprot_gn_id", relationship = "many-to-many") %>%
  dplyr::rename("source" = external_gene_name) %>%
  dplyr::filter(!is.na(source))

res <- getBM(attributes = c('uniprot_gn_id',
                            'external_gene_name'),
             values = unique(net$target_protein),
             mart = mart) %>%
  dplyr::rename("target_protein" = uniprot_gn_id)

net_df <- left_join(net_df, res, by = "target_protein", relationship = "many-to-many") %>%
  dplyr::rename("target_name" = external_gene_name) %>%
  dplyr::filter(!is.na(target_name)) %>%
  mutate(position = map_chr(str_split(target_site, "_"), 2)) %>%
  mutate(target_site = paste0(target_name, "_", position))

## Prepare data ---------------------------
hernandez_df <- read_csv(dataset)
identifiers <- hernandez_df %>%
  dplyr::select(ID) %>%
  dplyr::mutate(gene_name = map_chr(str_split(ID, "\\|"), 1))  %>%
  dplyr::mutate(position = map_chr(str_split(ID, "\\|"), 2)) %>%
  dplyr::mutate(ensg = map_chr(str_split(ID, "\\|"), 3)) %>%
  dplyr::mutate(ensp = map_chr(str_split(ID, "\\|"), 4)) %>%
  dplyr::mutate(target_site = paste0(gene_name, "_", position)) %>%
  dplyr::rename("site" = ID)


## Merge with network ---------------------------
ppsp_prior <- left_join(net_df %>%
                          dplyr::select(source, target_site, autophosphorylation),
                        identifiers, by = "target_site", relationship = "many-to-many")

ppsp_prior_df <- ppsp_prior %>%
  mutate(target = case_when(
    !is.na(site) ~ site,
    is.na(site) ~ target_site
  )) %>%
  mutate(target = case_when(
    autophosphorylation == "auto" ~ paste0(target, "|auto"), #mark autophosphorylation
    autophosphorylation == "no" ~ target
  )) %>%
  dplyr::mutate(mor = 1) %>%
  dplyr::select(source, target, mor)

ppsp_prior_df <- ppsp_prior_df %>%
  distinct(.keep_all = TRUE)

## Save processed phosphositeplus
write_tsv(ppsp_prior_df, output_file)
