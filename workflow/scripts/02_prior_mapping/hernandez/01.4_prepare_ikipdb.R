if(exists("snakemake")){
  ikip_file <- snakemake@input$ikip
  dataset <- snakemake@input$file_dataset
  output_file <- snakemake@output$tsv
}else{
  ikip_file <- "data/prior/iKiP-DB-Table.tsv"
  output_file <- "results/hernandez/prior/iKiPdb.tsv"
  dataset <- "results/hernandez/processed_data/benchmark_data.csv"
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
                                "ZAK" = "MAP3K20"
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


## merge with ikipdb ---------------------------
ikipd_prior <- left_join(ikipdb_df,
                         identifiers, by = "target_site", relationship = "many-to-many")

ikipd_prior_df <- ikipd_prior %>%
  mutate(target = case_when(
    !is.na(site) ~ site,
    is.na(site) ~ target_site
  )) %>%
  mutate(target = case_when(
    str_detect(pattern = kinase, string = target_site) ~ paste0(target, "|auto"), #mark autophosphorylation
    !str_detect(pattern = kinase, string = target_site) ~ target
  )) %>%
  dplyr::mutate(mor = 1) %>%
  dplyr::select(kinase, target, mor) %>%
  dplyr::rename("source" = kinase) %>%
  distinct()

rm_auto_duplicated <- ikipd_prior_df %>%
  mutate(target_auto = paste(source, str_remove(target, "\\|auto"), sep = "_")) %>%
  filter(str_detect(target, "auto")) %>%
  pull(target_auto)

ikipd_prior_df <- ikipd_prior_df %>%
  mutate(tmp = paste(source, target, sep = "_")) %>%
  filter(!tmp %in% rm_auto_duplicated) %>%
  dplyr::select(-tmp)

## Save processed ikipdb
write_tsv(ikipd_prior_df, output_file)
