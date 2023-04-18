if(!require("OmnipathR")) remotes::install_github('saezlab/OmnipathR', upgrade='never')
library(OmnipathR)
library(tidyverse)
library(biomaRt)

### Construct kinase-substrate interaction network
omnipath_ptm <- OmnipathR::get_signed_ptms()
omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation", "phosphorylation"), ]

# Filter out ProtMapper
omnipath_ptm_filtered <- omnipath_ptm %>%
  dplyr::filter(!(stringr::str_detect(omnipath_ptm$source, "ProtMapper") & n_resources == 1))
omnipath_ptm_filtered <- omnipath_ptm_filtered %>%
  mutate(target_site = paste0(substrate_genesymbol, "_", residue_type, residue_offset))

files_phospho <- list.files("data/CPTAC_phospho", full.names = T)
pps <- map_dfr(files_phospho, function(file){
  df <- read_tsv(file)
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

omnipath_prior <- left_join(target_df, omnipath_ptm_filtered %>% dplyr::select(target_site, enzyme_genesymbol, modification), by = "target_site", relationship = "many-to-many")

omnipath_prior_df <- omnipath_prior %>%
  drop_na() %>%
  mutate(mor = case_when(
    modification == "phosphorylation" ~ 1,
    modification == "dephosphorylation" ~ -1
  )) %>%
  dplyr::select(site, enzyme_genesymbol, mor) %>%
  dplyr::rename("Kinase" = enzyme_genesymbol, "Site2" = site)

write_tsv(omnipath_prior_df, "data/prior/omnipath.tsv")




