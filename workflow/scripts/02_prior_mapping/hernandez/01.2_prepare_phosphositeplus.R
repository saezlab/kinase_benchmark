if(exists("snakemake")){
  ppsp_file <- snakemake@input$ppsp
  dataset <- snakemake@input$file_dataset
  output_file <- snakemake@output$tsv
}else{
  ppsp_file <- "data/prior/phosphositeplus"
  output_file <- "results/hernandez/prior/phosphositeplus.tsv"
  dataset <- "results/hernandez/processed_data/benchmark_data.csv"
}

## Libraries ---------------------------
library(biomaRt)
library(tidyverse)
library(org.Hs.eg.db)

## Construct kinase-substrate interaction network ---------------------------
phosphositeplus <- read_tsv(ppsp_file, skip = 2)
phosphositeplus_human <- phosphositeplus %>%
  dplyr::filter(KIN_ORGANISM == "human" & SUB_ORGANISM == "human") %>%
  mutate(surrounding = toupper(`SITE_+/-7_AA`)) %>%
  mutate(target_site = paste0(SUB_GENE, "_", SUB_MOD_RSD))

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
ppsp_prior <- left_join(phosphositeplus_human %>%
                          dplyr::select(KINASE, target_site),
                        identifiers, by = "target_site", relationship = "many-to-many")

ppsp_prior_df <- ppsp_prior %>%
  mutate(target = case_when(
    !is.na(site) ~ site,
    is.na(site) ~ target_site
  ))  %>%
  dplyr::mutate(mor = 1) %>%
  dplyr::select(KINASE, target, mor) %>%
  dplyr::rename("source" = KINASE)

## Change kinases to common gene names
source_keys <- unique(ppsp_prior_df$source)
translate_df <- AnnotationDbi::select(org.Hs.eg.db, keys=toupper(source_keys),
                                      columns=c("ENTREZID","SYMBOL"), keytype="ALIAS")

## Manual renaming for kinases where the alias was not found
translate_df <- translate_df %>%
  dplyr::distinct(ALIAS, .keep_all = T)  %>%
  add_column(source = source_keys) %>%
  mutate(symbol_manual = recode(ALIAS,
                                "AMPKG2" = "PRKAG2",
                                "AMPKA1" = "PRKAA1",
                                "AMPKA2" = "PRKAA2",
                                "CK1E" = "CSNK1E",
                                "CAMK1A" = "CAMK1",
                                "BCR-ABL1" = "ABL1",
                                "YES" = "YES1",
                                "PKCB ISO2" = "PRKCB",
                                "CAMLCK" = "MYLK3",
                                "MARK3 ISO3" = "MARK3",
                                "CK1G2" = "CSNK1G2",
                                "PKCT" = "PRKCQ",
                                "CK1G1" = "CSNK1G1",
                                "CK1D" = "CSNK1D",
                                "MNK1 ISO2" = "MKNK1",
                                "PDHK4" = "PDK4",
                                "AURC" = "AURKC",
                                "MPSK1" = "STK16",
                                "P70S6KB" = "RPS6KB2",
                                "JNK2 ISO2" = "MAPK9",
                                "PKG1 ISO2" = "PRKG1",
                                "SKMLCK" = "MYLK2",
                                "DMPK1" = "DMPK",
                                "CAMK2D ISO8" = "CAMK2D",
                                "JNK1 ISO2" = "MAPK8",
                                "RET ISO3" = "RET",
                                "P38G" = "MAPK14",
                                "P38A" = "MAPK14",
                                "ALPHAK3" = "ALPK3",
                                "GSK3B ISO2" = "GSK3B",
                                "PKACA ISO2" = "PRKACA",
                                "PKCH" = "PRKCH",
                                "CAMK1B" = "PNCK",
                                "PKM ISO2" = "PKM",
                                "P90RSK" = "RPS6KA1",
                                "NPM-ALK" = "NA",
                                "AMPKB1" = "PRKAB1",
                                "CK1A" = "CSNK1A1",
                                "PDHK3" = "PDK3",
                                "P70S6K" = "RPS6KB1",
                                "P70S6K ISO2" = "RPS6KB1",
                                "P38D" = "MAPK14",
                                "CK1A2" = "CSNK1A1L",
                                "AURB" = "AURKB",
                                "KHS2" = "MAP4K3",
                                "CDK11A ISO10" = "CDK11A",
                                "PDHK1" = "PDK1",
                                "SMMLCK" = "MYLK",
                                "PKCZ" = "PRKCZ",
                                "CK1G3" = "CSNK1G3",
                                "BRSK1 ISO2" = "BRSK1"
  )) %>%
  mutate(symbol = case_when(
    !is.na(SYMBOL) ~ SYMBOL,
    is.na(SYMBOL) ~ symbol_manual
  )) %>%
  dplyr::filter(!symbol == "NA")


## Rename kinases in prior knowledge
ppsp_prior_df <- left_join(ppsp_prior_df,
                           translate_df %>%
                             dplyr::select(source, symbol), by = "source") %>%
  dplyr::filter(!is.na(symbol)) %>%
  dplyr::mutate(source = symbol) %>%
  dplyr::select(-symbol) %>%
  mutate(target = case_when(
    str_detect(pattern = source, string = target) ~ paste0(target, "|auto"), #mark autophosphorylation
    !str_detect(pattern = source, string = target) ~ target
  )) %>%
  distinct()


## Save processed phosphositeplus
write_tsv(ppsp_prior_df, output_file)
