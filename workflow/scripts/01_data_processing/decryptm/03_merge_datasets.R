#' Merging decryptm metadata and matrix from different dataset

## Snakemake ---------------------------
if(exists("snakemake")){
  EC50_files <- snakemake@input$EC50_files
  meta_files <- snakemake@input$meta_files
  targets <- snakemake@input$targets
  output_meta <- snakemake@output$output_meta
  output_EC50 <- snakemake@output$output_EC50
  output_R2_EC50 <- snakemake@output$output_R2_EC50
  output_drugs <- snakemake@output$output_drugs
}else{
  EC50_files <- list.files("results/decryptm/phosphoproteome/", pattern = "^EC50", full.names = T)
  meta_files <- list.files("results/decryptm/phosphoproteome/", pattern = "^meta", full.names = T)
  targets <- "data/decryptm/drug_targets.csv"
  output_meta <- "results/decryptm/processed_data/meta_data.csv"
  output_EC50 <- "results/decryptm/processed_data/pEC50.csv"
  output_R2_EC50 <- "results/decryptm/processed_data/R2_pEC50.csv"
  output_drugs <- "results/decryptm/processed_data/overview_drugs.csv"
}

## Libraries ---------------------------
library(tidyverse)

## Merge metadata ---------------------------
meta <- map_dfr(meta_files, read_csv)
meta <- meta %>%
  mutate(sample = recode(sample,
                         "ddPTM_A549_Dasatinib_60min_R1_Rep_Dasatinib" = "ddPTM_A549_Dasatinib_60min_R1_Rep",
                         "ddPTM_A549_PD325901_60min_R1_Rep_PD325901" = "ddPTM_A549_PD325901_60min_R1_Rep",
                         "ddPTM_A549_Refametinib_60min_R1_Rep_Refametinib" = "ddPTM_A549_Refametinib_60min_R1_Rep",
                         "ddPTM_A549_Staursporin_60min_R1_Rep_Staursporin" = "ddPTM_A549_Staursporin_60min_R1_Rep",))

target_df <- read_csv(targets) %>%
  rename("drug" = "Drug_name")

drug_targets <- target_df %>%
  filter(!MoA == "None") %>%
  select(-Target_type) %>%
  group_by(drug, Drug_type, MoA) %>%
  summarise(Target = paste(UniProtID, collapse = ";"),
            Target_name = paste(Target_name, collapse = ";")) %>%
  mutate(drug = recode(drug,
                       "MK-2206" = "MK2206",
                       "AZD-4547" = "AZD4547",
                       "AZD-8055" = "AZD8055",
                       "PD-325901" = "PD325901",
                       "Nintendanib" = "Nintedanib",
                       "Staurosporine" = "Staursporin"))

meta <- left_join(meta, drug_targets, relationship = "many-to-many")

# Find molar mass for drugs which are not yet in mol/l
meta %>%
  filter(dose_unit == "ng/ml") %>%
  pull(drug) %>%
  unique()

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7493232/
# Pertuzumab: 148556 g/mol
# Trastuzumab: 148531 g/mol

meta <- meta %>%
  mutate(molecular_weight = case_when(
    drug == "Pertuzumab" ~ 148556,
    drug == "Trastuzumab" ~ 148531
  )) %>%
  mutate(dose_scale_factor = case_when(
    dose_unit == "ng/ml" ~ 1e-6 / molecular_weight,
    dose_unit == "M" ~ 1
  ))


## Merge metadata and EC50 matrices ---------------------------
EC50_df <- map_dfr(EC50_files, read_csv)

EC50_corr_df <- EC50_df %>%
  left_join(meta %>% select(sample, dose_scale_factor), by = "sample") %>%
  mutate(EC50 = EC50*dose_scale_factor) %>%
  mutate(pEC50 = -log10(EC50)) %>%
  mutate(R2_pEC50 = R2*pEC50) %>%
  mutate(regulation = sign(Curve.effect.size)) %>%
  mutate(pEC50_sign = pEC50*regulation) %>%
  mutate(R2_pEC50_sign = R2_pEC50*regulation) %>%
  select(-c(dose_scale_factor, regulation))

pEC50_wide <- EC50_corr_df %>%
  select(pps_id, sample, pEC50_sign) %>%
  pivot_wider(values_from = pEC50_sign, names_from = sample)

R2_pEC50_wide <- EC50_corr_df %>%
  select(pps_id, sample, R2_pEC50_sign) %>%
  pivot_wider(values_from = R2_pEC50_sign, names_from = sample)

colnames(pEC50_wide)[!colnames(pEC50_wide) %in% meta$sample]
meta$sample[!meta$sample %in% colnames(pEC50_wide)]

meta <- meta %>%
  filter(sample %in% colnames(pEC50_wide))

table_drugs <- table(meta$drug, meta$cells) %>%
  data.frame() %>%
  pivot_wider(names_from = Var2, values_from = Freq) %>%
  rename("drug" = Var1)

## Save_data ---------------------------
write_csv(meta, output_meta)
write_csv(pEC50_wide, output_EC50)
write_csv(R2_pEC50_wide, output_R2_EC50)
write_csv(table_drugs, output_drugs)
