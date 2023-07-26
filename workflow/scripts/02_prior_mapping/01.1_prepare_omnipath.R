if(exists("snakemake")){
  decryptm_dataset <- snakemake@input$decryptm
  file_datasets <- snakemake@input$file_dataset
  output_file_decryptm <- snakemake@output$out_decryptm
  output_file <- snakemake@output$tsv
}else{
  file_datasets <- list.files("data/CPTAC_phospho", full.names = T)
  output_file <- "results/prior/omnipath.tsv"
  decryptm_dataset <- "results/decryptm/processed_data/R2_pEC50.csv"
  output_file_decryptm <- "results/decryptm/prior/omnipath.tsv"
  if(!require("OmnipathR")) remotes::install_github('saezlab/OmnipathR', upgrade='never')
}

## Libraries ---------------------------
library(OmnipathR)
library(tidyverse)
library(biomaRt)

## Construct kinase-substrate interaction network ---------------------------
omnipath_ptm <- OmnipathR::get_signed_ptms()
omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation", "phosphorylation"), ]

# Filter out ProtMapper
omnipath_ptm_filtered <- omnipath_ptm %>%
  dplyr::filter(!(stringr::str_detect(omnipath_ptm$source, "ProtMapper") & n_resources == 1))
omnipath_ptm_filtered <- omnipath_ptm_filtered %>%
  mutate(target_site = paste0(substrate_genesymbol, "_", residue_type, residue_offset))

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
                  position = map_chr(str_split(pps$site,  "\\|"), 3))

## Convert Ensemble gene IDs to Gene names
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

res <- getBM(attributes = c('ensembl_gene_id',
                            'external_gene_name'),
             values = identifiers,
             mart = mart)

target_df <- full_join(pps_df, res, by = "ensembl_gene_id") %>%
  mutate(target_site = paste0(external_gene_name, "_", position))

## Prepare decryptm ---------------------------
decryptm_df <- read_csv(decryptm_dataset)
decryptm_identifiers <- decryptm_df %>%
  dplyr::select(pps_id) %>%
  dplyr::mutate(gene_name = map_chr(str_split(pps_id, "\\|"), 2))  %>%
  dplyr::mutate(position = map_chr(str_split(pps_id, "\\|"), 3)) %>%
  dplyr::mutate(target_site = paste0(gene_name, "_", position)) %>%
  dplyr::rename("site" = pps_id)


## Merge with network ---------------------------
## CPTAC ---------------------------
omnipath_prior <- left_join(omnipath_ptm_filtered %>% dplyr::select(target_site, enzyme_genesymbol, modification), target_df, by = "target_site", relationship = "many-to-many")

omnipath_prior_df <- omnipath_prior %>%
  mutate(target = case_when(
    !is.na(site) ~ site,
    is.na(site) ~ target_site
  )) %>%
  mutate(target = case_when(
    str_detect(pattern = enzyme_genesymbol, string = target_site) ~ paste0(target, "|auto"), #mark autophosphorylation
    !str_detect(pattern = enzyme_genesymbol, string = target_site) ~ target
  )) %>%
  mutate(mor = case_when(
    modification == "phosphorylation" ~ 1,
    modification == "dephosphorylation" ~ -1
  )) %>%
  dplyr::select(enzyme_genesymbol, target, mor) %>%
  drop_na() %>%
  dplyr::rename("source" = enzyme_genesymbol) %>%
  distinct()

# Remove edges with duplicated sign information
edges <- paste(omnipath_prior_df$source, omnipath_prior_df$target, sep = "_")
dup_edges <- edges[duplicated(edges)]

omnipath_prior_df <- omnipath_prior_df %>%
  dplyr::mutate(KS = paste(source, target, sep = "_")) %>%
  dplyr::filter(!KS %in% dup_edges) %>%
  dplyr::select(-KS)

## Save processed Omnipath
write_tsv(omnipath_prior_df, output_file)

## decryptm ---------------------------
omnipath_prior <- left_join(omnipath_ptm_filtered %>% dplyr::select(target_site, enzyme_genesymbol, modification), decryptm_identifiers, by = "target_site", relationship = "many-to-many")

omnipath_prior_df <- omnipath_prior %>%
  mutate(target = case_when(
    !is.na(site) ~ site,
    is.na(site) ~ target_site
  )) %>%
  mutate(target = case_when(
    str_detect(pattern = enzyme_genesymbol, string = target_site) ~ paste0(target, "|auto"), #mark autophosphorylation
    !str_detect(pattern = enzyme_genesymbol, string = target_site) ~ target
  )) %>%
  mutate(mor = case_when(
    modification == "phosphorylation" ~ 1,
    modification == "dephosphorylation" ~ -1
  )) %>%
  dplyr::select(enzyme_genesymbol, target, mor) %>%
  drop_na() %>%
  dplyr::rename("source" = enzyme_genesymbol) %>%
  distinct()

# Remove edges with duplicated sign information
edges <- paste(omnipath_prior_df$source, omnipath_prior_df$target, sep = "_")
dup_edges <- edges[duplicated(edges)]

omnipath_prior_df <- omnipath_prior_df %>%
  dplyr::mutate(KS = paste(source, target, sep = "_")) %>%
  dplyr::filter(!KS %in% dup_edges) %>%
  dplyr::select(-KS)

## Save processed Omnipath
write_tsv(omnipath_prior_df, output_file_decryptm)


