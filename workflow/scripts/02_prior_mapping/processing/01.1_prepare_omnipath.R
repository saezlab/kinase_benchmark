if(exists("snakemake")){
  output_file <- snakemake@output$tsv
}else{
  output_file <- "results/prior/raw/omnipath.tsv"
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

## Merge with network ---------------------------
omnipath_prior_df <- omnipath_ptm_filtered %>%
  dplyr::select(target_site, enzyme_genesymbol, modification) %>%
  rename("source" = enzyme_genesymbol, "target" = target_site) %>%
  mutate(target_protein = map_chr(str_split(target, "_"), 1)) %>%
  mutate(position = map_chr(str_split(target, "_"), 2)) %>%
  mutate(mor = case_when(
    modification == "phosphorylation" ~ 1,
    modification == "dephosphorylation" ~ -1
  )) %>%
  dplyr::select(source, target, target_protein, position, mor) %>%
  distinct() %>%
  mutate(sequence = NA)

# Remove edges with duplicated sign information
edges <- paste(omnipath_prior_df$source, omnipath_prior_df$target, sep = "_")
dup_edges <- edges[duplicated(edges)]

omnipath_prior_df <- omnipath_prior_df %>%
  dplyr::mutate(KS = paste(source, target, sep = "_")) %>%
  dplyr::filter(!KS %in% dup_edges) %>%
  dplyr::filter(mor == 1) %>%
  dplyr::select(-KS)

## Save processed Omnipath
write_tsv(omnipath_prior_df, output_file)


