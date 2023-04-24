if(exists("snakemake")){
  ikip_file <- snakemake@input$ikip
  output_file <- snakemake@output$tsv
  file_datasets <- snakemake@input$file_dataset
}else{
  ikip_file <- "data/prior/iKiP-DB-Table.tsv"
  output_file <- "results/prior/iKiPdb.tsv"
  file_datasets <- list.files("data/CPTAC_phospho", full.names = T)
}

## Libraries ---------------------------
library(biomaRt)
library(tidyverse)

## Construct kinase-substrate interaction network ---------------------------
ikipdb <- read_tsv(ikip_file)

ikipdb_df <- map_dfr(unique(ikipdb$Kinase_name), function(kinase){
    ids <- ikipdb %>%
      dplyr::filter(Kinase_name == kinase) %>%
      pull(ids) %>%
      str_split(pattern = ";") %>%
      unlist()

    sites <- ikipdb %>%
      dplyr::filter(Kinase_name == kinase) %>%
      pull(sites) %>%
      str_split(pattern = ";") %>%
      unlist()

    data.frame(kinase = kinase,
               id = ids,
               sequence = sites) %>%
      mutate(target_protein = map_chr(str_split(id, "\\."), 1)) %>%
      mutate(target_protein = map_chr(str_split(target_protein, "-"), 1)) %>%
      mutate(aa = map_chr(str_split(id, "\\."), 2)) %>%
      mutate(position = map_chr(str_split(id, "\\."), 3)) %>%
      mutate(target_site = paste0(target_protein, "_", sequence))
  })


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
                            'uniprot_gn_id'),
             values = identifiers,
             mart = mart)

target_df <- full_join(pps_df, res, by = "ensembl_gene_id", relationship = "many-to-many") %>%
  mutate(target_site = paste0(uniprot_gn_id, "_", surrounding))

## merge with phosphositeplus
ikipd_prior <- left_join(ikipdb_df %>%
                          dplyr::select(kinase, target_site),
                        target_df, by = "target_site", relationship = "many-to-many")

ikipd_prior_df <- ikipd_prior %>%
  mutate(target = case_when(
    !is.na(site) ~ site,
    is.na(site) ~ target_site
  )) %>%
  dplyr::mutate(mor = 1) %>%
  dplyr::select(kinase, target, mor) %>%
  dplyr::rename("source" = kinase)


## Rename kinases in prior knowledge ---------------------------
## Change kinases to common gene names
ikipd_prior_df <- ikipd_prior_df %>%
  dplyr::mutate(source = recode(source,
                                "ACK" = "TNK2",
                                "AMPKA1" = "PRKAA1",
                                "AMPKA2" = "PRKAA2",
                                "CAMK1A" = "CAMK1",
                                "CDC2" = "CDK1",
                                "CGK2" = "PRKG2",
                                "LYNA" = "LYN",
                                "LYNB" = "LYN",
                                "MLCK" = "MYLK",
                                "NPM1:ALK" = "NA",
                                "PRKCB1" = "PRKCB",
                                "PRKCB2" = "PRKCB",
                                "RIPK5" = "DSTYK",
                                "TSSK1" = "TSSK1B",
                                "ZAK" = "MAP3K20"
                                )) %>%
  dplyr::filter(!source == "NA") %>%
  distinct()

## Save processed phosphositeplus ---------------------------
write_tsv(ikipd_prior_df, output_file)
