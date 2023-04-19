if(exists("snakemake")){
  ppsp_file <- snakemake@input$ppsp
  output_file <- snakemake@output$tsv
  file_datasets <- snakemake@input$file_dataset
}else{
  ppsp_file <- "data/prior/phosphositeplus"
  output_file <- "results/prior/phosphositeplus.tsv"
  file_datasets <- list.files("data/CPTAC_phospho", full.names = T)
}

## Libraries ---------------------------
library(biomaRt)
library(tidyverse)

## Construct kinase-substrate interaction network ---------------------------
phosphositeplus <- read_tsv(ppsp_file, skip = 2)
phosphositeplus_human <- phosphositeplus %>%
  dplyr::filter(KIN_ORGANISM == "human" & SUB_ORGANISM == "human") %>%
  mutate(surrounding = toupper(`SITE_+/-7_AA`)) %>%
  mutate(target_site = paste0(SUB_GENE, "_", surrounding))

# Map targets to pps in data
pps <- map_dfr(file_datasets, function(file){
  df <- read_tsv(file)
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
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

res <- getBM(attributes = c('ensembl_gene_id',
                            'external_gene_name'),
             values = identifiers,
             mart = mart)

target_df <- full_join(pps_df, res, by = "ensembl_gene_id") %>%
  mutate(target_site = paste0(external_gene_name, "_", surrounding))

## merge with phosphositeplus
ppsp_prior <- left_join(target_df, phosphositeplus_human %>% dplyr::select(KINASE, target_site), by = "target_site", relationship = "many-to-many")

ppsp_prior_df <- ppsp_prior %>%
  drop_na() %>%
  dplyr::mutate(mor = 1) %>%
  dplyr::select(KINASE, site, mor) %>%
  dplyr::rename("source" = KINASE, "target" = site) %>%
  distinct()


## Save processed Omnipath ---------------------------
write_tsv(ppsp_prior_df, output_file)

