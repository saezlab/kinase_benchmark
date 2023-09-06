if(exists("snakemake")){
  networkin_file <- snakemake@input$networkin_file
  dataset <- snakemake@input$file_dataset
  output_file <- snakemake@output$tsv
  networkin_score <- snakemake@params$score
}else{
  networkin_file <- "data/prior/networkin_human_predictions_3.1.tsv"
  output_file <- "results/hernandez/prior/networkin.tsv"
  dataset <- "results/hernandez/processed_data/benchmark_data.csv"
  networkin_score <- 5
}
networkin_score <- as.numeric(networkin_score)

## Libraries ---------------------------
library(tidyverse)

## load NetworKIN ---------------------------
# load NetworKin prediction scores; convert kinase names to gene symbol and filter to kinases in KS table
nwkin <- read.table(networkin_file, sep = "\t",  stringsAsFactors = F, header= T, comment.char = "")

# rename kinases
nwkin$id <- toupper(nwkin$id)
nwkin <- nwkin %>%
  mutate(id = recode(id,
                     "ABL" = "ABL1",
                     "AURORAA" = "AURKA",
                     "AURORAA" = "AURKA",
                     "AURORAB" = "AURKB",
                     "AURORAC" = "AURKC",
                     "BRK" = "PTK6",
                     "CAMKIALPHA" = "CAMK1",
                     "CAMKIIALPHA" = "CAMK2A",
                     "CAMKIIBETA" = "CAMK2B",
                     "CAMKIIDELTA" = "CAMK2D",
                     "CAMKIIGAMMA" = "CAMK2G",
                     "CAMKIV" = "CAMK4",
                     "CAMKIDELTA" = "CAMK1D",
                     "CK1ALPHA" = "CSNK1A1",
                     "CK1DELTA" = "CSNK1D",
                     "CK1EPSILON" = "CSNK1E",
                     "CK1GAMMA1" = "CSNK1G1",
                     "CK1GAMMA2" = "CSNK1G2",
                     "CK1GAMMA3" = "CSNK1G3",
                     "CK2ALPHA" = "CSNK2A1",
                     "CK2A2" = "CSNK2A2",
                     "DNAPK" = "PRKDC",
                     "GSK3ALPHA" = "GSK3A",
                     "GSK3BETA" = "GSK3B",
                     "IKKALPHA" = "CHUK",
                     "IKKBETA" = "IKBKB",
                     "LKB1" = "STK11",
                     "MRCKA" = "CDC42BPA",
                     "MST2" = "STK3",
                     "MRCKA" = "CDC42BPA",
                     "PDGFRBETA" = "PDGFRB",
                     "PDGFRALPHA" = "PDGFRA",
                     "PDHK1" = "PDK1",
                     "PDHK2" = "PDK2",
                     "PDHK3" = "PDK3",
                     "PKBALPHA" = "AKT1",
                     "PKBBETA" = "AKT2",
                     "PKBGAMMA" = "AKT3",
                     "PKCALPHA" = "PRKCA",
                     "PKCEPSILON" = "PRKCE",
                     "PKCETA" = "PRKCH",
                     "PKCGAMMA" = "PRKCG",
                     "PKCIOTA" = "PRKCI",
                     "PKCTHETA" = "PRKCQ",
                     "PKG1CGKI" = "PRKG1",
                     "TRKA" = "NTRK1",
                     "TRKB" = "NTRK2",
                     "YES" = "YES1",
                     "PKAALPHA" = "PRKACA",
                     "PKABETA" = "PRKACB",
                     "PKAGAMMA" = "PRKACG"))

# Select targets above certain NetworKIN score
nwkin_filtered <- nwkin[nwkin$networkin_score >= networkin_score, ]
nwkin_filtered$AA <- unlist(str_extract_all(nwkin_filtered$sequence, '[[:lower:]]+'))

nwkin_filtered <- nwkin_filtered %>%
  mutate(target = map_chr(str_split(X.substrate, " "), 1)) %>%
  mutate(target_site = paste0(target, "_", toupper(AA), position))

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
nwkin_df <- left_join(nwkin_filtered, identifiers, by = "target_site", relationship = "many-to-many") %>%
  mutate(target = case_when(
    is.na(site) ~ target_site,
    !is.na(site) ~ site
  )) %>%
  dplyr::rename("source" = id) %>%
  dplyr::select(source, target) %>%
  mutate(target = case_when(
    str_detect(target, source) ~ paste0(target, "|auto"), #mark autophosphorylation
    !str_detect(target, source) ~ target
  )) %>%
  add_column(mor = 1) %>%
  distinct()

## Save prior
write_tsv(nwkin_df, output_file)


