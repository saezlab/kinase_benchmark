if(exists("snakemake")){
  ptmsig_file <- snakemake@input$ptmsig
  file_datasets <- snakemake@input$file_dataset
  output_file <- snakemake@output$tsv
}else{
  ptmsig_file <- "data/prior/ptm.sig.db.all.uniprot.human.v1.9.0.gmt"
  output_file <- "results/prior/ptmsigdb.tsv"
  file_datasets <- list.files("data/CPTAC_phospho", full.names = T)
}

## Libraries ---------------------------
library(biomaRt)
library(qusage)
library(tidyverse)

## Construct kinase-substrate interaction network ---------------------------
PTMsigDB <- qusage::read.gmt(ptmsig_file)
PTMsigDB.kinase <- PTMsigDB[str_detect(names(PTMsigDB), "KINASE")]

PTMsig_df <- map_dfr(names(PTMsigDB.kinase), function(kin){
  data.frame(source = kin,
             target_site = PTMsigDB.kinase[[kin]])
})

PTMsig_df <- PTMsig_df %>%
  mutate(source = str_remove(source, "KINASE-PSP_")) %>%
  mutate(target_site = str_remove(target_site, "-p;u"))

## Change kinases to common gene names---------------------------
PTMsig_df$source <- map_chr(str_split(paste(PTMsig_df$source,
                                            PTMsig_df$source,
                                            sep = "/"),
                                      "/"),
                            2)

# Add uniprot id of kinases to identify autophosphorylation
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
res_kin <- getBM(attributes = c('external_gene_name',
                                'uniprot_gn_id'),
                 values = PTMsig_df$source,
                 mart = mart)  %>%
  dplyr::rename("source" = external_gene_name) %>%
  dplyr::rename("source_uniprot" = uniprot_gn_id)

PTMsig_df <- left_join(PTMsig_df, res_kin, by = "source", relationship = "many-to-many")

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

target_df <- full_join(pps_df, res, by = "ensembl_gene_id", relationship = "many-to-many") %>%
  mutate(target_site = paste0(uniprot_gn_id, ";", position))


## Merge with network ---------------------------
## CPTAC ---------------------------
PTMsig_prior <- left_join(PTMsig_df,
                        target_df, by = "target_site", relationship = "many-to-many")

PTMsig_prior_df <- PTMsig_prior %>%
  mutate(target = case_when(
    !is.na(site) ~ site,
    is.na(site) ~ target_site
  )) %>%
  mutate(target = case_when(
    str_detect(pattern = source_uniprot, string = target_site) ~ paste0(target, "|auto"), #mark autophosphorylation
    !str_detect(pattern = source_uniprot, string = target_site) ~ target
  )) %>%
  dplyr::mutate(mor = 1) %>%
  dplyr::select(source, target, mor)


# Remove duplicated edges (if present)
PTMsig_prior_df <- PTMsig_prior_df %>%
  distinct()

## Save processed PTMsigDB
write_tsv(PTMsig_prior_df, output_file)
