if(exists("snakemake")){
  ikip_file <- snakemake@input$ikip
  output_file <- snakemake@output$tsv
}else{
  ikip_file <- "data/prior/iKiP-DB-Table.tsv"
  output_file <- "results/prior/raw/iKiPdb.tsv"
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
      mutate(target_site = paste0(target_protein, "_", aa, position))
  })

## Rename kinases in prior knowledge ---------------------------
## Change kinases to common gene names
ikipdb_df <- ikipdb_df %>%
  dplyr::mutate(kinase = recode(kinase,
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
                                "ZAK" = "MAP3K20",
                                "ICK" = "CILK1"
  )) %>%
  dplyr::filter(!kinase == "NA") %>%
  distinct()

# Add uniprot id of kinases to identify autophosphorylation
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

res_kin <- getBM(attributes = c('uniprot_gn_id',
                                'external_gene_name'),
                 values = ikipdb_df$target_protein,
                 mart = mart)  %>%
  dplyr::rename("target_protein" = uniprot_gn_id) %>%
  dplyr::rename("target_gene_name" = external_gene_name)

ikipdb_df <- left_join(ikipdb_df, res_kin, by = "target_protein", relationship = "many-to-many") %>%
  filter(!is.na(target_gene_name))
ikipdb_df <- ikipdb_df %>%
  mutate(target_site = paste0(target_gene_name, "_", aa, position))

## merge with ikipdb ---------------------------
ikipd_prior_df <- ikipdb_df %>%
  dplyr::select(kinase, sequence, target_site) %>%
  dplyr::rename("source" = kinase, "target" = target_site) %>%
  mutate(target_protein = map_chr(str_split(target, "_"), 1)) %>%
  mutate(position = map_chr(str_split(target, "_"), 2)) %>%
  mutate(mor = 1) %>%
  dplyr::select(source, target, target_protein, position, mor, sequence) %>%
  distinct()

## Save processed ikipdb
write_tsv(ikipd_prior_df, output_file)
