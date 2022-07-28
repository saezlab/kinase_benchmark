# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will check the difference in the kinase activity estimation
#' 1. prepare output files to merge them
#' 2. correlation
#' 3. Jaccard

library(tidyverse)

## Prepare output ---------------------------

## RoKAI
rokai_act_files <- list.files("output/rokai", pattern = "output", full.names = T)
rokai_names <- str_remove(str_remove(rokai_act_files, "output/rokai/rokai_output_"), ".csv")
rokai_act <- map(rokai_act_files, read_csv)
names(rokai_act) <- rokai_names

## KSEA
KSEA_act_files <- list.files("output/KSEA", pattern = "output", full.names = T)
KSEA_names <- str_remove(str_remove(KSEA_act_files, "output/KSEA/KSEA_output_"), ".csv")
KSEA_act <- map(KSEA_act_files, read_csv)
names(KSEA_act) <- KSEA_names

## KEA3
KEA3_act_files <- list.files("output/KEA3", pattern = "top_rank", full.names = T)
KEA3_names <- str_remove(str_remove(KEA3_act_files, "output/KEA3/KEA3_output_"), "_top_rank.tsv")
KEA3_act <- map(KEA3_act_files, read_tsv)
names(KEA3_act) <- KEA3_names

## IKAP
IKAP_act_files <- list.files("output/IKAP", pattern = "kinase", full.names = T)
IKAP_names <- str_remove(str_remove(IKAP_act_files, "output/IKAP/kinase_act_"), ".csv")
IKAP_names <- str_remove(IKAP_names, "output/IKAP/kinase_")

IKAP_act <- map(unique(IKAP_names), function(cell_line){
  files <- list.files("output/IKAP", pattern = cell_line, full.names = T)
  act_df <- cbind(read_csv(files[3], col_names = F), read_csv(files[2], col_names = F))
  names(act_df) <- c("Gene", "Act")
  act_df
})

names(IKAP_act) <- unique(IKAP_names)

## Correlation ---------------------------
act_cells <- map(names(IKAP_act), function(cell_line){
  rokai <- rokai_act[[cell_line]] %>%
    select(Gene, ZScore)
  KSEA <- KSEA_act[[cell_line]] %>%
    select(Kinase.Gene, z.score) %>%
    rename("Gene" = "Kinase.Gene")
  IKAP <- IKAP_act[[cell_line]]

  activity_df <- list(rokai, KSEA, IKAP) %>% reduce(full_join, by = "Gene") %>%
    column_to_rownames("Gene")
  names(activity_df) <- c("ROKAI", "KSEA", "IKAP")

  activity_df
})

names(act_cells) <- names(IKAP_act)

corr_cells <- map(act_cells, function(act_df){
  cor(act_df, method = "pearson", use = "complete.obs")
})


## Jaccard Index ---------------------------
top_kinase_cells <- map(names(IKAP_act), function(cell_line){
  rokai <- rokai_act[[cell_line]] %>%
    arrange(desc(abs(ZScore))) %>%
    slice(1:20) %>%
    pull(Gene)

  KSEA <- KSEA_act[[cell_line]] %>%
    arrange(desc(abs(z.score))) %>%
    slice(1:20) %>%
    pull(Kinase.Gene)

  IKAP <- IKAP_act[[cell_line]] %>%
    arrange(desc(abs(Act))) %>%
    slice(1:20) %>%
    pull(Gene)

  KEA3 <- KEA3_act[[cell_line]] %>%
    slice(1:20) %>%
    pull(Protein)
  KEA3 <- str_remove(string = KEA3, pattern = "\\*")

 list(rokai = rokai,
      KSEA = KSEA,
      IKAP = IKAP,
      KEA3 = KEA3)
})

names(corr_cells) <- names(IKAP_act)

jacc_cells <- map(top_kinase_cells, calcJaccard)

map(jacc_cells, get_mat_plot)

map(corr_cells, function(x){get_mat_plot(x, main = 'Pearson correlation', palette = "Greens")})
