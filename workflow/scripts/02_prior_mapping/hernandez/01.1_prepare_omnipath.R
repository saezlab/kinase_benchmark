if(exists("snakemake")){
  dataset <- snakemake@input$file_dataset
  output_file <- snakemake@output$tsv
}else{
  dataset <- "results/hernandez/processed_data/benchmark_data.csv"
  output_file <- "results/hernandez/prior/omnipath.tsv"
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

# Filter out KEA3
omnipath_ptm_filtered <- omnipath_ptm_filtered %>%
  dplyr::filter(!(stringr::str_detect(omnipath_ptm_filtered$source, "KEA") & n_resources == 1))

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
omnipath_prior <- left_join(omnipath_ptm_filtered %>%
                              dplyr::select(target_site, enzyme_genesymbol, modification),
                            identifiers,
                            by = "target_site", relationship = "many-to-many")

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


