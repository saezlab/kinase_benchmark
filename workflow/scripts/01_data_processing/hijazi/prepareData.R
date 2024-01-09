#'

## Snakemake ---------------------------
if(exists("snakemake")){
  raw_HL60 <- snakemake@input$HL60
  raw_MCF7 <- snakemake@input$MCF7
  p_HL60 <- snakemake@input$HL60p
  p_MCF7 <- snakemake@input$MCF7p
  meta_file <- snakemake@input$meta
  targets_file <- snakemake@input$targets
  bench_output <- snakemake@output$mat
  meta_prior_output <- snakemake@output$meta_out
  meta_discoverx_output <- snakemake@output$meta_discover
}else{
  raw_HL60 <- "data/hijazi/HL60_fc.tsv"
  p_HL60 <- "data/hijazi/HL60_pval.tsv"
  raw_MCF7 <- "data/hijazi/MCF7_fc.tsv"
  p_MCF7 <- "data/hijazi/MCF7_pval.tsv"
  meta_file <- "data/hijazi/inhibitor_selectivity.tsv"
  targets_file <- "data/hijazi/targets_hijazi.tsv"
  bench_output <- "results/hijazi/01_processed_data/benchmark_data.csv"
  meta_prior_output <- "results/hijazi/01_processed_data/benchmark_metadataPrior.csv"
  meta_discoverx_output <- "results/hijazi/01_processed_data/benchmark_metadataDiscoverX.csv"
}

## Libraries ---------------------------
library(tidyverse)

## Prepare data ---------------------------
# LogFC of benchmark data HL60
HL60_df <- read_tsv(raw_HL60, col_types = cols())
HL60_pval_df <- read_tsv(p_HL60, col_types = cols())

HL60p_long <- HL60_pval_df %>%
  pivot_longer(!sh.index.sites, names_to = "experiment", values_to = "pVal") %>%
  mutate(experiment = str_remove(experiment, "\\.p\\.value"))

HL60FC_long <-HL60_df %>%
  pivot_longer(!sh.index.sites, names_to = "experiment", values_to = "logFC") %>%
  mutate(experiment = str_remove(experiment, "\\.fold"))

HL60_mono <- full_join(HL60FC_long, HL60p_long, by = c("sh.index.sites", "experiment")) %>%
  tidyr::separate_rows(sh.index.sites, sep = ";") %>%
  filter(!sh.index.sites == "")

HL60_mono_df <- HL60_mono %>%
  dplyr::group_by(sh.index.sites, experiment) %>%
  filter(pVal == min(pVal)) %>%
  ungroup() %>%
  dplyr::select(-pVal) %>%
  distinct(sh.index.sites, experiment, .keep_all = T) %>%
  pivot_wider(names_from = "experiment", values_from = "logFC") %>%
  mutate(protein = map_chr(str_split(sh.index.sites, "\\("), 1), .after = sh.index.sites) %>%
  mutate(aa = map_chr(str_split(sh.index.sites, "\\("), 2) %>% str_remove("\\)"), .after = protein) %>%
  mutate(id = paste0(protein, "_", aa, "|", protein, "|", aa), .before = sh.index.sites) %>%
  dplyr::select(-protein, -aa, -sh.index.sites)


# LogFC of benchmark data MCF7
MCF7_df <- read_tsv(raw_MCF7, col_types = cols())
MCF7_pval_df <- read_tsv(p_MCF7, col_types = cols())

MCF7p_long <- MCF7_pval_df %>%
  pivot_longer(!sh.index.sites, names_to = "experiment", values_to = "pVal") %>%
  mutate(experiment = str_remove(experiment, "\\.p\\.value"))

MCF7FC_long <-MCF7_df %>%
  pivot_longer(!sh.index.sites, names_to = "experiment", values_to = "logFC") %>%
  mutate(experiment = str_remove(experiment, "\\.fold"))

MCF7_mono <- full_join(MCF7FC_long, MCF7p_long, by = c("sh.index.sites", "experiment")) %>%
  tidyr::separate_rows(sh.index.sites, sep = ";") %>%
  filter(!sh.index.sites == "")

MCF7_mono_df <- MCF7_mono %>%
  dplyr::group_by(sh.index.sites, experiment) %>%
  filter(pVal == min(pVal)) %>%
  ungroup() %>%
  dplyr::select(-pVal) %>%
  distinct(sh.index.sites, experiment, .keep_all = T) %>%
  pivot_wider(names_from = "experiment", values_from = "logFC") %>%
  mutate(protein = map_chr(str_split(sh.index.sites, "\\("), 1), .after = sh.index.sites) %>%
  mutate(aa = map_chr(str_split(sh.index.sites, "\\("), 2) %>% str_remove("\\)"), .after = protein) %>%
  mutate(id = paste0(protein, "_", aa, "|", protein, "|", aa), .before = sh.index.sites) %>%
  dplyr::select(-protein, -aa, -sh.index.sites)

phospho_df <- full_join(HL60_mono_df, MCF7_mono_df, by = "id") %>%
  mutate(id = str_remove_all(id, ";")) %>%
  column_to_rownames("id")

phospho_df[phospho_df == 0] <- NA
phospho_df <- phospho_df[!rowSums(is.na(phospho_df)) == ncol(phospho_df),]


# MetaData
metaData <- data.frame(id = colnames(phospho_df)) %>%
  mutate(cell_line = map_chr(str_split(id, "\\."), 1)) %>%
  mutate(drug = map_chr(str_split(id, "\\."), 2)) %>%
  add_column(sign = -1)

targets_manual <- read_tsv(targets_file, col_types = cols()) %>%
  mutate(target = case_when(
    is.na(Target_TTD) & !is.na(Manual) ~ Manual,
    !is.na(Target_TTD) & is.na(Manual) ~ Target_TTD,
    is.na(Target_TTD) & is.na(Manual) ~ NA
  )) %>%
  filter(!is.na(target))

metaData_prior <- left_join(metaData, targets_manual %>% dplyr::select(drug, target), by = "drug") %>%
  filter(!is.na(target))


## Select targets based on drug selectivity data
drug_selectivity <- read_tsv(meta_file, col_types = cols())

drug_selectivity_df <- drug_selectivity %>%
  pivot_longer(!kinase, names_to = "drug_id", values_to = "selectivity") %>%
  mutate(dataset = map_chr(str_split(drug_id, "\\.\\."), 1)) %>%
  mutate(drug = map_chr(str_split(drug_id, "\\.\\."), 2)) %>%
  mutate(drug = str_remove_all(drug, "-")) %>%
  mutate(drug = recode(drug,
                       "Silmitasertib" = "CX4945",
                       "Pictilisib" = "GDC0941",
                       "Ribociclib" = "LEE011",
                       "Abemaciclib" = "LY2835219",
                       "Amuvatinib" = "MP470"))

targets_discoverX <- drug_selectivity_df %>%
  filter(dataset == "DiscoverX") %>%
  filter(selectivity < 50) %>%
  distinct(kinase, drug) %>%
  group_by(drug) %>%
  summarise(target = paste(kinase, collapse = ";"))

metaData_discoverX <- left_join(metaData, targets_discoverX, by = "drug") %>%
  filter(!is.na(target))

# Save data
write_csv(phospho_df %>% rownames_to_column("ID"), bench_output)
write_csv(metaData_prior, meta_prior_output)
write_csv(metaData_discoverX, meta_discoverx_output)

