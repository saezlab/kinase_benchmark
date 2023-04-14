# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we estimate kinase
#' activities with all combinations of
#' KSN and method

library(tidyverse)
library(decoupleR)
library(tidyselect)
source("code/support_functions.R")
map(list.files("code", pattern = "run", full.names = T), source)

## Load data---------------------------
### phosphoproteomics
ccrcc <- read_csv(file = "data/CPTAC_estefani_clean_zip/ccrcc/ccrcc_phospho_tvalues.csv") %>%
  dplyr::select(ID, NATvsTUM_t) %>%
  rename("ccrcc_tVal" = "NATvsTUM_t")
luad <- read_csv(file = "data/CPTAC_estefani_clean_zip/luad/luad_phospho_tvalues.csv") %>%
  dplyr::select(ID, NATvsTUM_t) %>%
  rename("luad_tVal" = "NATvsTUM_t")
ucec <- read_csv(file = "data/CPTAC_estefani_clean_zip/ucec/ucec_phospho_tvalues.csv") %>%
  dplyr::select(ID, NATvsTUM_t) %>%
  rename("ucec_tVal" = "NATvsTUM_t")

mat <- reduce(list(ccrcc, luad, ucec), full_join) %>%
  filter(rowSums(is.na(.)) == 0) %>%
  column_to_rownames("ID")


### networks
network_files <- list.files("data/networks", full.names = T)
networks <- map(network_files, read_csv)
names(networks) <- str_remove(list.files("data/networks"), ".csv")

## ---------------------------
kin_scores <- map_dfr(names(networks), function(network_i){
  print(network_i)
  network <- networks[[network_i]] %>%
    distinct()

  #prepare network
  network <- intersect_regulons(mat, network, .source = "source", .target = "target", minsize = 5)
  cor.source <- check_corr(network) %>% filter(correlation > 0.9) %>% pull(source.2)
  network <- network %>% filter(!source %in% cor.source) %>% ungroup()

  KARP <- run_KARP(mat, network)
  RoKAI_z <- run_zscore_RoKAI(mat, network)
  KSEA_z <- run_zscore_KSEA(mat, network)
  INKA <- run_INKA(mat, network)
  Rokai_lm <- run_lm_rokai(mat, network)
  decoupler <- decouple(mat = as.matrix(mat), network = network)
  fgsea <- run_fgsea(mat = as.matrix(mat), network = network)
  gsva <- run_gsva(mat = as.matrix(mat), network = network)
  wmean <- run_wmean(mat = as.matrix(mat), network = network)
  viper <- run_viper(mat = as.matrix(mat), network = network)

  decoupler <- decoupler %>%
    select(c(source, condition, score, statistic)) %>%
    rename("method" = "statistic")
  fgsea <- fgsea %>%
    select(c(source, condition, score, statistic)) %>%
    rename("method" = "statistic")
  gsva <- gsva %>%
    select(c(source, condition, score, statistic)) %>%
    rename("method" = "statistic")
  wmean <- wmean %>%
    select(c(source, condition, score, statistic)) %>%
    rename("method" = "statistic")
  viper <- viper %>%
    select(c(source, condition, score, statistic)) %>%
    rename("method" = "statistic")

  rbind(KARP, RoKAI_z, KSEA_z, INKA, Rokai_lm, decoupler, fgsea, gsva, wmean, viper) %>%
    add_column(network = network_i)
})

gmt_f <- qusage::read.gmt("data/ptm.sig.db.all.uniprot.human.v1.9.0.gmt")
data_df <- map_dfr(names(gmt_f), function(name_i){
  data.frame(source = name_i, target = gmt_f[[name_i]], mor = 1)
})

id_uniprot <- map_chr(str_split(data_df$target, ";"), 1)
pos <- str_remove(map_chr(str_split(data_df$target, ";"), 2), "-p")

map_df <- OmnipathR::uniprot_id_mapping_table(id_uniprot, uniprot, genesymbol)

data_df <- cbind(data_df, From = id_uniprot, pos)

data_df <- full_join(map_df, data_df, by = "From", multiple = "all")
final_df <- data_df %>%
  mutate(target = paste0(data_df$To, "_", data_df$pos)) %>%
  dplyr::select(source, target, mor) %>%
  distinct(source, target, mor, .keep_all = TRUE)

fgsea <- run_fgsea(mat = as.matrix(mat["ccrcc_tVal"]), network = final_df)
res1 <- fgsea %>%
  filter(condition == "ccrcc_tVal") %>%
  arrange(desc(score)) %>%
  filter(statistic == "norm_fgsea") %>%
  rename("id" = source)

all <- full_join(res1, res2, by = "id") %>%
  filter(is.numeric(score))

all <- all %>%
  mutate(t = as.numeric(all$t)) %>%
  drop_na() %>%
  filter(is.numeric(t))

cor(all$score, all$t)

kin_scores <- kin_scores %>%
  filter(condition == "ccrcc_tVal") %>%
  filter(network == "networkin_5")

kin_scores <- kin_scores %>%
  add_column(combination = paste(kin_scores$method, kin_scores$network, kin_scores$condition, sep = "_"))

scores <- kin_scores %>%
  select(source, score, combination) %>%
  pivot_wider(names_from = combination, values_from = score) %>%
  column_to_rownames("source")

cor_scores <- cor(scores,  method = "pearson", use = "complete.obs")
idx_rows <- !rowSums(is.na(cor_scores)) >= (ncol(cor_scores)-1)
idx_cols <- !colSums(is.na(cor_scores)) >= (nrow(cor_scores)-1)

cor_scores <- cor_scores[idx_rows, idx_cols]

cor_scores_list <- map(colnames(cor_scores), function(x){
  idx <- cor_scores[,x] > 0.98
  rownames(cor_scores)[idx]
})
names(cor_scores_list) <- colnames(cor_scores)


ann_col <- kin_scores %>%
  select(combination, method, network) %>%
  distinct() %>%
  column_to_rownames("combination")

paletteLength <- 100
myColor <- colorRampPalette(c("steelblue3", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(cor_scores), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(cor_scores)/paletteLength, max(cor_scores), length.out=floor(paletteLength/2)))

networkcolors <- RColorBrewer::brewer.pal(n = length(unique(ann_col$network)), name = "Dark2")
names(networkcolors) <- unique(ann_col$network)

methodcolors <- RColorBrewer::brewer.pal(n = length(unique(ann_col$method)), name = "Set3")
names(methodcolors) <- unique(ann_col$method)
mycolors <- list(method = methodcolors,
                 network = networkcolors)



# Plot the heatmap
pheatmap::pheatmap(cor_scores, display_numbers = round(cor_scores, digits = 2), annotation_col = ann_col, annotation_row = ann_col,
                   color=myColor, breaks=myBreaks,
                   cutree_cols = 4,
                   cutree_rows = 4,
                   show_rownames = F, show_colnames = F,
                   annotation_colors = mycolors)




