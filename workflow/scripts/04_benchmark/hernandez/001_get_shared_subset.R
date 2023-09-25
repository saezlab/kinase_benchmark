#' Process data for benchmark

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input$rds
  output_file <- snakemake@output$output
  methods <- snakemake@params$methods
}else{
  input_file <- list.files("results/hernandez/final_scores/scaled", pattern = "GPSppsp", full.names = T)
  methods <- c("KARP", "KS", "KSEA_z", "PC1", "RoKAI_z", "UQ", "Wilcox", "fgsea", "mean", "median", "mlm", "norm_fgsea","norm_wmean","number_of_targets","rokai_lm","ulm","viper","wsum", "ptmsea")
}

## Libraries ---------------------------
library(tidyverse)

prior <- str_remove(str_remove(input_file, "results/hernandez/final_scores/scaled/"), ".rds")
scores <- map(input_file, readRDS)

names(scores) <- prior

scores_long <- map_dfr(names(scores), function(score_i){
  map_dfr(methods, function(method_i){
    scores[[score_i]][[method_i]] %>%
      rownames_to_column("kinase") %>%
      pivot_longer(!kinase, names_to = "sample", values_to = "score") %>%
      add_column(prior = score_i) %>%
      add_column(method = method_i) %>%
      filter(!is.na(score))
  })
})

th_all <- length(unique(scores_long$method)) * length(unique(scores_long$prior))

msk <- scores_long %>%
  group_by(kinase, sample) %>%
  summarise(n = n()) %>%
  filter(n >= th_all) %>%
  mutate(tmp = paste(kinase, sample, sep = ":")) %>%
  pull(tmp)

scores_filtered <- scores_long %>%
  mutate(tmp = paste(kinase, sample, sep = ":")) %>%
  filter(tmp %in% msk) %>%
  select(-tmp)

map(prior, function(prior_i){
  scores_prior <- scores_filtered %>%
    filter(prior == prior_i) %>%
    group_by(method)

  score_list <- scores_prior %>%
    group_split(.keep = F)

  names(score_list) <- group_keys(scores_prior)$method

  score_list_wide <- map(names(score_list), function(met){
    score_list[[met]] %>%
      select(kinase, sample, score) %>%
      pivot_wider(names_from = sample, values_from = score) %>%
      column_to_rownames("kinase")
  })
  names(score_list_wide) <- names(score_list)

  saveRDS(score_list_wide, paste0("results/hernandez/final_scores/subset/", prior_i, ".rds"))
})
