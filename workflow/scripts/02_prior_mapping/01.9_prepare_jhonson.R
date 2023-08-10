if(exists("snakemake")){
  ppsp_file <- snakemake@input$ppsp
  file_datasets <- snakemake@input$file_dataset
  decryptm_dataset <- snakemake@input$decryptm
  output_file <- snakemake@output$tsv
  output_file_decryptm <- snakemake@output$out_decryptm
}else{
  ppsp_file <- "data/prior/simplified_jhonson_with_psite.csv"
  output_file <- "results/prior/jhonson.tsv"
  file_datasets <- list.files("data/CPTAC_phospho", full.names = T)
  decryptm_dataset <- "results/decryptm/processed_data/R2_pEC50.csv"
  output_file_decryptm <- "results/decryptm/prior/jhonson.tsv"
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


## Prepare CPTAC ---------------------------
# Map targets to pps in data
pps <- map_dfr(file_datasets, function(file){
  df <- read_tsv(file, col_types = cols())
  data.frame(site = df$site)
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
res <- getBM(attributes = c('ensembl_gene_id',
                            'uniprot_gn_id'),
             values = identifiers,
             mart = mart)

target_df <- full_join(pps_df, res, by = "ensembl_gene_id") %>%
  mutate(target_site = paste0(uniprot_gn_id, "_", position))


## Prepare decryptm ---------------------------
decryptm_df <- read_csv(decryptm_dataset)
decryptm_identifiers <- decryptm_df %>%
  dplyr::select(pps_id) %>%
  mutate(gene_name = map_chr(str_split(pps_id, "\\|"), 1))  %>%
  mutate(position = map_chr(str_split(pps_id, "\\|"), 3)) %>%
  mutate(target_site = paste0(gene_name, "_", position)) %>%
  dplyr::rename("site" = pps_id)

## Merge with network ---------------------------
## CPTAC ---------------------------
ppsp_prior <- left_join(net_df %>%
                          dplyr::select(source, target_site, autophosphorylation),
                        target_df, by = "target_site", relationship = "many-to-many")

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

## Save processed phosphositeplus
write_tsv(ppsp_prior_df, output_file)

## decryptm ---------------------------
ppsp_prior <- left_join(net_df %>%
                          dplyr::select(source, target_site, autophosphorylation),
                        decryptm_identifiers, by = "target_site", relationship = "many-to-many")

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

## Save processed phosphositeplus
write_tsv(ppsp_prior_df, output_file_decryptm)
