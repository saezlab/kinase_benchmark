# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will check the difference in the kinase activity estimation
#' 1. prepare output files to merge them
#' 2. correlation
#' 3. Jaccard

library(tidyverse)
library(pheatmap)
library(ggplotify)
source("code/support_functions.R")

## Prepare output ---------------------------

## RoKAI
rokai_act_files <- list.files("output/CPTAC/rokai", pattern = "output", full.names = T)
rokai_names <- str_remove(str_remove(rokai_act_files, "output/CPTAC/rokai/kinase_output_"), ".csv")
rokai_act <- map(rokai_act_files, read_csv)
names(rokai_act) <- rokai_names

## KSEA
KSEA_act_files <- list.files("output/CPTAC/KSEA", pattern = "output", full.names = T)
KSEA_names <- str_remove(str_remove(KSEA_act_files, "output/CPTAC/KSEA/KSEA_output_"), ".csv")
KSEA_act <- map(KSEA_act_files, read_csv)
names(KSEA_act) <- KSEA_names

## KEA3
KEA3_act_files <- list.files("output/CPTAC/KEA3", pattern = "top_rank", full.names = T)
KEA3_names <- str_remove(str_remove(KEA3_act_files, "output/CPTAC/KEA3/KEA3_output_"), "_top_rank.tsv")
KEA3_act <- map(KEA3_act_files, read_tsv)
names(KEA3_act) <- KEA3_names

## IKAP
IKAP_act_files <- list.files("output/CPTAC/IKAP", pattern = "kinase", full.names = T)
IKAP_names <- str_remove(str_remove(IKAP_act_files, "output/CPTAC/IKAP/kinase_act_"), ".csv")
IKAP_names <- str_remove(IKAP_names, "output/CPTAC/IKAP/kinase_")

IKAP_act <- map(unique(IKAP_names), function(cohort){
  files <- list.files("output/CPTAC/IKAP", pattern = cohort, full.names = T)
  act_df <- cbind(read_csv(files[3], col_names = F), read_csv(files[2], col_names = F))
  names(act_df) <- c("Gene", "Act")
  act_df
})

names(IKAP_act) <- unique(IKAP_names)

## PTM-SEA
PTM_act_files <- list.files("output/CPTAC/PTM-SEA/2022-08-15", pattern = "ssGSEA-scores.gct", full.names = T, recursive = TRUE)
PTM_names <- map_chr(str_split(PTM_act_files, "_"), function(x) x[3])
PTM_act <- map(PTM_act_files, function(x) read.delim(file = x, skip=2))
names(PTM_act) <- PTM_names


## Scatter plot ---------------------------
act_cells <- map(names(IKAP_act), function(cohort){
  rokai <- rokai_act[[cohort]] %>%
    filter(!duplicated(Gene)) %>%
    dplyr::select(Gene, ZScore)
  KSEA <- KSEA_act[[cohort]] %>%
    dplyr::select(Kinase.Gene, z.score) %>%
    dplyr::rename("Gene" = "Kinase.Gene")
  IKAP <- IKAP_act[[cohort]]
  PTM <- PTM_act[[cohort]]
  PTM <- PTM[str_detect(PTM_act[[cohort]]$id, "KINASE"), c("id", "t")] %>%
    rename("Gene" = "id")
  PTM$Gene <- map_chr(str_split(PTM$Gene, "_"), function(x) x[2])
  PTM$Gene <-  map_chr(str_split(PTM$Gene, "/"), function(x) x[length(x)])

  activity_df <- list(rokai, KSEA, IKAP, PTM) %>% reduce(full_join, by = "Gene") %>%
    column_to_rownames("Gene")
  names(activity_df) <- c("ROKAI", "KSEA", "IKAP", "PTM")

  activity_df
})
names(act_cells) <- names(IKAP_act)


map(names(act_cells), function(cohort){
  ggplot(act_cells[[cohort]], aes(x=KSEA, y=ROKAI, color = "ROKAI")) +
    geom_point() +
    geom_point(aes(x=KSEA, y=IKAP, color = "IKAP")) +
    geom_point(aes(x=KSEA, y=PTM, color = "PTM-SEA")) +
    geom_smooth(aes(x=KSEA, y=ROKAI),method=lm, se=FALSE, color = "#619CFF", fullrange=TRUE) +
    geom_smooth(aes(x=KSEA, y=IKAP), method=lm, se=FALSE, color = "#F8766D", fullrange=TRUE) +
    geom_smooth(aes(x=KSEA, y=PTM), method=lm, se=FALSE, color = "#00BA38", fullrange=TRUE) +
    xlab("KSEA activity scores") +
    ylab("kinase activity scores") +
    labs(color='Methods') +
    theme_minimal() +
    theme(text = element_text(size = 26))
})

ccrcc <- act_cells[[cohort]]
colnames(ccrcc)
map(colnames(ccrcc), function(test){
  tmp <- ccrcc %>%
    rownames_to_column("Gene")

  new <- tmp[c("Gene", test)]
  new <- new[!is.na(new[2]),]
  new[2] <- rank(new$ROKAI)
})
rokai <- ccrcc %>%
  rownames_to_column("Gene")

rank <- data.frame(ROKAI = rank(abs(act_cells[[cohort]]$ROKAI)*-1, na.last = "keep"),
                   KSEA = rank(abs(act_cells[[cohort]]$KSEA)*-1, na.last = "keep"),
                   IKAP = rank(abs(act_cells[[cohort]]$IKAP)*-1, na.last = "keep"),
                   PTM = rank(abs(act_cells[[cohort]]$PTM)*-1, na.last = "keep"))

rownames(rank) <- rownames(act_cells[[cohort]])


ggplot(rank, aes(x=KSEA, y=ROKAI, color = "ROKAI")) +
  geom_point() +
  geom_point(aes(x=KSEA, y=IKAP, color = "IKAP")) +
  geom_point(aes(x=KSEA, y=PTM, color = "PTM-SEA")) +
  geom_smooth(aes(x=KSEA, y=ROKAI),method=lm, se=FALSE, color = "#619CFF", fullrange=TRUE) +
  geom_smooth(aes(x=KSEA, y=IKAP), method=lm, se=FALSE, color = "#F8766D", fullrange=TRUE) +
  geom_smooth(aes(x=KSEA, y=PTM), method=lm, se=FALSE, color = "#00BA38", fullrange=TRUE) +
  xlab("KSEA activity ranks") +
  ylab("kinase activity ranks") +
  labs(color='Methods') +
  theme_minimal() +
  theme(text = element_text(size = 26))

rokai_act
KSEA_act
IKAP_act
PTM_act

KEA3_act

## Correlation ---------------------------
act_cells <- map(names(IKAP_act), function(cohort){
  rokai <- rokai_act[[cohort]] %>%
    filter(!duplicated(Gene)) %>%
    dplyr::select(Gene, ZScore)
  KSEA <- KSEA_act[[cohort]] %>%
    dplyr::select(Kinase.Gene, z.score) %>%
    dplyr::rename("Gene" = "Kinase.Gene")
  IKAP <- IKAP_act[[cohort]]
  PTM <- PTM_act[[cohort]]
  PTM <- PTM[str_detect(PTM_act[[cohort]]$id, "KINASE"), c("id", "t")] %>%
    rename("Gene" = "id")
  PTM$Gene <- map_chr(str_split(PTM$Gene, "_"), function(x) x[2])
  PTM$Gene <-  map_chr(str_split(PTM$Gene, "/"), function(x) x[length(x)])

  activity_df <- list(rokai, KSEA, IKAP, PTM) %>% reduce(full_join, by = "Gene") %>%
    column_to_rownames("Gene")
  names(activity_df) <- c("ROKAI", "KSEA", "IKAP", "PTM")

  activity_df
})

names(act_cells) <- names(IKAP_act)

corr_cells <- map(act_cells, function(act_df){
  cor(act_df, method = "pearson", use = "complete.obs")
})

names(corr_cells) <- names(act_cells)
## Jaccard Index ---------------------------
n <- 5
top_kinase_cells <- map(names(IKAP_act), function(cohort){
  rokai <- rokai_act[[cohort]] %>%
    arrange(desc(abs(ZScore))) %>%
    slice(1:n) %>%
    pull(Gene)

  KSEA <- KSEA_act[[cohort]] %>%
    arrange(desc(abs(z.score))) %>%
    slice(1:n) %>%
    pull(Kinase.Gene)

  IKAP <- IKAP_act[[cohort]] %>%
    arrange(desc(abs(Act))) %>%
    slice(1:n) %>%
    pull(Gene)

  KEA3 <- KEA3_act[[cohort]] %>%
    slice(1:n) %>%
    pull(Protein)
  KEA3 <- str_remove(string = KEA3, pattern = "\\*")

  PTM <- PTM_act[[cohort]]
  PTM <- PTM[str_detect(PTM_act[[cohort]]$id, "KINASE"), c("id", "t")] %>%
    rename("Gene" = "id")
  PTM$Gene <- map_chr(str_split(PTM$Gene, "_"), function(x) x[2])
  PTM$Gene <-  map_chr(str_split(PTM$Gene, "/"), function(x) x[length(x)])
  PTM <- PTM %>%
    arrange(desc(abs(t))) %>%
    slice(1:n) %>% pull(Gene)

 list(rokai = rokai,
      KSEA = KSEA,
      IKAP = IKAP,
      KEA3 = KEA3,
      PTM = PTM)
})

names(top_kinase_cells) <- names(IKAP_act)

x <- top_kinase_cells$ccrcc
x <- x[c("rokai", "KSEA", "IKAP", "PTM")]
names(x) <- c("RoKAI", "KSEA", "IKAP", "PTM-SEA")
library(ggvenn)
ggvenn(
  x,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 5, text_size = 5
)

jacc_cells <- map(top_kinase_cells, calcJaccard)
jacc_cells$ccrcc
top_kinase_cells$ccrcc

map(jacc_cells, get_mat_plot)

map(corr_cells, function(x){get_mat_plot(x, main = 'Pearson correlation', palette = "Spectral")})

mean_corr_matrix <- matrix(ncol = ncol(corr_cells$ccrcc),
                           nrow = nrow(corr_cells$ccrcc))

for (i in 1:ncol(corr_cells$ccrcc)){
  for (j in 1:nrow(corr_cells$ccrcc)){
    mean_corr_matrix[j,i] <-  mean(c(corr_cells$ccrcc[j,i],
                                  corr_cells$luad[j,i],
                                  corr_cells$ucec[j,i]))

  }

}

colnames(mean_corr_matrix) <- colnames(corr_cells$ccrcc)
rownames(mean_corr_matrix) <- rownames(corr_cells$ccrcc)
get_mat_plot(mean_corr_matrix, main = 'Pearson correlation', palette = "Greens")
