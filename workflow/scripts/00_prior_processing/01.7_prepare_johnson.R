if(exists("snakemake")){
  ppsp_file <- snakemake@input$ppsp
  tyr_file <- snakemake@input$tyr
  par <- snakemake@params$perc
  output_file <- snakemake@output$tsv
}else{
  ppsp_file <- "data/kinase_libraries/johnson_library/pps_fifteenmer_ser_thr_percent.tsv"
  tyr_file <- "data/kinase_libraries/johnson_library/pps_fifteenmer_tyr_percent.tsv"
  par <- 15
  output_file <- "results/prior/raw/johnson99.tsv"
}
par <- as.numeric(par)
par_quant <- par >= 100
if (par_quant){
  par <- as.numeric(par)/100
}


## Libraries ---------------------------
library(tidyverse)
library(org.Hs.eg.db)

## Construct kinase-substrate interaction network ---------------------------
net <- read_tsv(ppsp_file, col_types = cols())
net_long <- net %>%
  mutate(id = paste(pps, fifteenmer, sep = "|")) %>%
  dplyr::select(-ID, -pps, -fifteenmer, -phos_res, -`SITE_+/-7_AA`) %>%
  pivot_longer(!id, names_to = "source", values_to = "percentile")

net_tyr <- read_tsv(tyr_file, col_types = cols())
tyr_long <- net_tyr %>%
  mutate(id = paste(pps, fifteenmer, sep = "|")) %>%
  dplyr::select(-ID, -pps, -fifteenmer, -phos_res, -`SITE_+/-7_AA`) %>%
  pivot_longer(!id, names_to = "source", values_to = "percentile")

full_net <- rbind(net_long, tyr_long) %>%
  add_column(mor = 1)

if (par_quant){
  johnson_df <- full_net %>%
    filter(percentile >= par) %>%
    mutate(target = map_chr(str_split(id, "\\|"), 1)) %>%
    mutate(target_protein = map_chr(str_split(target, "_"), 1)) %>%
    mutate(position = map_chr(str_split(target, "_"), 2)) %>%
    mutate(sequence = map_chr(str_split(id, "\\|"), 2)) %>%
    dplyr::select(source, target, target_protein, position, mor, sequence) %>%
    distinct(.keep_all = TRUE)
} else {
  johnson_df <- full_net %>%
    group_by(id) %>%
    arrange(desc(percentile)) %>%
    dplyr::slice(1:par) %>%
    ungroup() %>%
    mutate(target = map_chr(str_split(id, "\\|"), 1)) %>%
    mutate(target_protein = map_chr(str_split(target, "_"), 1)) %>%
    mutate(position = map_chr(str_split(target, "_"), 2)) %>%
    mutate(sequence = map_chr(str_split(id, "\\|"), 2)) %>%
    dplyr::select(source, target, target_protein, position, mor, sequence) %>%
    distinct(.keep_all = TRUE)
}


## Change kinases to common gene names
source_keys <- unique(johnson_df$source)
translate_df <- AnnotationDbi::select(org.Hs.eg.db, keys=toupper(source_keys),
                                      columns=c("ENTREZID","SYMBOL"), keytype="ALIAS")

## Manual renaming for kinases where the alias was not found
translate_df <- translate_df %>%
  dplyr::distinct(ALIAS, .keep_all = T)  %>%
  add_column(source = source_keys)  %>%
  mutate(symbol = case_when(
    !is.na(SYMBOL) ~ SYMBOL,
    is.na(SYMBOL) ~ ALIAS
  )) %>%
  dplyr::filter(!symbol == "NA")

johnson_df <- left_join(johnson_df,
                           translate_df %>%
                             dplyr::select(source, symbol), by = "source") %>%
  dplyr::filter(!is.na(symbol)) %>%
  dplyr::mutate(source = symbol) %>%
  dplyr::select(-symbol) %>%
  distinct()

johnson_df <- johnson_df %>%
  mutate(source = recode(source,
                         "ALPHAK3" = "ALPK3",
                         "CK1G2" = "CSNK1G2",
                         "CK1D" = "CSNK1D",
                         "MPSK1" = "STK16",
                         "CK1A" = "CSNK1A1",
                         "SKMLCK" = "MYLK2",
                         "PDHK4" = "PDK4",
                         "CK1A2" = "CSNK1A1L",
                         "CK1E" = "CSNK1E",
                         "CK1G1" = "CSNK1G1",
                         "STLK3" = "STK39",
                         "AURB" = "AURKB",
                         "CAMK1B" = "PNCK",
                         "CK1G3" = "CSNK1G3",
                         "P70S6K" = "RPS6KB1",
                         "P70S6KB" = "RPS6KB2",
                         "P90RSK" = "RPS6KA1",
                         "PKACG" = "PRKACG",
                         "KHS2" = "MAP4K3",
                         "TAO3" = "TAOK3",
                         "P38A" = "MAPK14",
                         "P38D" = "MAPK14",
                         "P38G" = "MAPK14",
                         "PDHK1" = "PDK1",
                         "TSSK1A" = "TSSK1B",
                         "CAMK1A" = "CAMK1",
                         "PKCH" = "PRKCH",
                         "AURC" = "AURKC",
                         "CAMLCK" = "MYLK3",
                         "SMMLCK" = "MYLK",
                         "AMPKA1" = "PRKAA1",
                         "AMPKA2" = "PRKAA2",
                         "IRE2" = "ERN2",
                         "PKCT" = "PRKCQ",
                         "PKCZ" = "PRKCZ",
                         "DMPK1" = "DMPK",
                         "YES" = "YES1"))

## Save processed phosphositeplus
write_tsv(johnson_df, output_file)
