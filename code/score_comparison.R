# Copyright (c) Sophia Müller-Dott [2023]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script I will compare the mean results of me and Eric

library(tidyverse)

## ---------------------------
brca_eric_combined <- readRDS("data/activity_score_eric/kinase_activity_scores_brca_v6_combined.Rds")
brca_sophia_combined <- readRDS("output/cptac/brca_combined.rds")

brca_eric_gS <- readRDS("data/activity_score_eric/kinase_activity_scores_brca_v6.Rds")
brca_sophia_gS <- readRDS("output/cptac/brca_goldStandard.rds")

hnscc_eric_combined <- readRDS("data/activity_score_eric/kinase_activity_scores_hnscc_v6_combined.Rds")
hnscc_sophia_combined <- readRDS("output/cptac/hnscc_combined.rds")

hnscc_eric_gS <- readRDS("data/activity_score_eric/kinase_activity_scores_brca_v6.Rds")
hnscc_sophia_gS <- readRDS("output/cptac/brca_goldStandard.rds")


## brca combined ---------------------------
mean_brca_eric_combined <- brca_eric_combined$mean %>%
  as.data.frame() %>%
  rownames_to_column("kinase") %>%
  pivot_longer(!kinase, values_to = "activity_eric", names_to = "sample")

mean_brca_sophia_combined <- brca_sophia_combined$wmean %>%
  as.data.frame() %>%
  rownames_to_column("kinase") %>%
  pivot_longer(!kinase, values_to = "activity_sophia", names_to = "sample")

merged_brca_combined <- full_join(mean_brca_eric_combined, mean_brca_sophia_combined, by = c("kinase", "sample"))
#View(merged_brca_combined)
merged_brca_combined$activity_eric[1]
merged_brca_combined$activity_sophia[1]

## hnscc combined ---------------------------
mean_hnscc_eric_combined <- hnscc_eric_combined$mean %>%
  as.data.frame() %>%
  rownames_to_column("kinase") %>%
  pivot_longer(!kinase, values_to = "activity_eric", names_to = "sample")

mean_hnscc_sophia_combined <- hnscc_sophia_combined$wmean %>%
  as.data.frame() %>%
  rownames_to_column("kinase") %>%
  pivot_longer(!kinase, values_to = "activity_sophia", names_to = "sample")

merged_hnscc_combined <- full_join(mean_hnscc_eric_combined, mean_hnscc_sophia_combined, by = c("kinase", "sample"))
#View(merged_hnscc_combined)
merged_hnscc_combined$activity_eric[1]
merged_hnscc_combined$activity_sophia[1]

## brca gold standard ---------------------------
mean_brca_eric_gS <- brca_eric_gS$mean %>%
  as.data.frame() %>%
  rownames_to_column("kinase") %>%
  pivot_longer(!kinase, values_to = "activity_eric", names_to = "sample")

mean_brca_sophia_gS <- brca_sophia_gS$wmean %>%
  as.data.frame() %>%
  rownames_to_column("kinase") %>%
  pivot_longer(!kinase, values_to = "activity_sophia", names_to = "sample")

merged_brca_gS <- full_join(mean_brca_eric_gS, mean_brca_sophia_gS, by = c("kinase", "sample"))
#View(merged_brca_gS)
merged_brca_gS$activity_eric[1]
merged_brca_gS$activity_sophia[1]

## hnscc gold standard ---------------------------
mean_hnscc_eric_gS <- hnscc_eric_gS$mean %>%
  as.data.frame() %>%
  rownames_to_column("kinase") %>%
  pivot_longer(!kinase, values_to = "activity_eric", names_to = "sample")

mean_hnscc_sophia_gS <- hnscc_sophia_gS$wmean %>%
  as.data.frame() %>%
  rownames_to_column("kinase") %>%
  pivot_longer(!kinase, values_to = "activity_sophia", names_to = "sample")

merged_hnscc_gS <- full_join(mean_hnscc_eric_gS, mean_hnscc_sophia_gS, by = c("kinase", "sample"))
#View(merged_hnscc_gS)
merged_hnscc_gS$activity_eric[1]
merged_hnscc_gS$activity_sophia[1]


## Gold Standard comparison ---------------------
brca_goldStandard <- read_tsv("data/brca_gold_standard_targets.tsv") %>%
  mutate(edge = paste(Kinase, Site2, sep = ":"))
hnscc_goldStandard <- read_tsv("data/hnscc_gold_standard_targets.tsv") %>%
  mutate(edge = paste(Kinase, Site2, sep = ":"))
all_goldStandard <- read_tsv("data/prior/allCT_gold_standard_targets.tsv") %>%
  mutate(edge = paste(Kinase, Site2, sep = ":"))

brca <- read_tsv("data/CPTAC_phospho/brca_phospho_data_median_centered.tsv")


# Check difference between combined goldStandard and cancer specififc gold Standard
sum(brca$site %in% brca_goldStandard$Site2)
sum(brca$site %in% all_goldStandard$Site2)

all_goldStandard
any(!hnscc_goldStandard$edge %in% all_goldStandard$edge)

# Let's test this for one specific kinase
kin <- "AKT1" #all_goldStandard$Kinase[1] #ABL1
brca_goldStandard %>% filter(Kinase == "ABL1") %>%
  filter(Site2 %in% brca$site)
targets_kin <- all_goldStandard %>% filter(Kinase == kin) %>%
  pull(Site2)

brca_kin <- brca %>%
  filter(site %in% targets_kin) %>%
  column_to_rownames("site")

samp <- colnames(brca_kin)[1] #X11BR047
brca_kin_samp <- brca_kin[samp] %>%
  drop_na() %>%
  colMeans()

brca_kin_samp
merged_brca_gS %>%
  filter(sample == samp) %>%
  filter(kinase == kin)




