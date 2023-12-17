#' Process data for benchmark

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file_raw <- snakemake@input$rds_raw
  input_file_scaled <- snakemake@input$rds_scaled
  output_file <- snakemake@output$output
  output_file_scaled <- snakemake@output$output_scaled
  methods <- snakemake@params$methods
}else{
  input_file <- list.files("results/hernandez/final_scores/subset", full.names = T)
  methods <- c("KARP", "KS", "KSEA_z", "PC1", "RoKAI_z", "Wilcox", "fgsea", "mean", "median", "mlm" ,"norm_wmean","number_of_targets","rokai_lm","ulm","viper","wsum", "ptmsea")
}

## Libraries ---------------------------
library(tidyverse)
library(corrplot)

methods <- methods[!methods == "number_of_targets"]

act_list <- map(input_file, readRDS)
names(act_list) <- str_remove(str_remove(input_file, ".rds"), "results/hernandez/final_scores/subset/")
act_list <- map(act_list, function(x) x[methods])
act_list <- act_list[c("iKiPdb", "networkin", "omnipath", "GPS", "phosphositeplus", "ptmsigdb")]

method_comparison <- map(names(act_list), function(prior_idx){
  scores <- act_list[[prior_idx]]

  scores_all <- map(names(scores), function(method_idx){
    score_long <- scores[[method_idx]] %>%
      rownames_to_column("kinase") %>%
      pivot_longer(!kinase, names_to = "sample", values_to = "score") %>%
      mutate(id = paste(kinase, sample, sep = ":")) %>%
      dplyr::select(id, score)

    colnames(score_long)[2] <- method_idx
    score_long
  }) %>%
    purrr::reduce(full_join, by = "id") %>%
    column_to_rownames("id")

  cor(scores_all, use = "complete.obs")
})

names(method_comparison) <- names(act_list)

Y <- do.call(cbind, method_comparison)
Y_array <- array(Y, dim=c(dim(method_comparison[[1]]), length(method_comparison)))

mean_correlation <- colMeans(aperm(Y_array, c(3, 1, 2)), na.rm = TRUE)
colnames(mean_correlation) <- colnames(method_comparison[[1]])
rownames(mean_correlation) <- rownames(method_comparison[[1]])

pdf("results/manuscript_figures/figure_2/corrplot_methods.pdf", height = 3, width = 3.5)
corrplot(mean_correlation, col.lim=c(0, 1), is.corr = FALSE, type = 'lower', tl.col = 'black', order = 'hclust', cl.pos = 'r')
dev.off()

## Prior comparison ---------------------------
prior_comparison <- map(names(act_list[[1]]), function(method_idx){
  scores <- map(act_list, function(x) x[[method_idx]])

  scores_all <- map(names(scores), function(method_idx){
    score_long <- scores[[method_idx]] %>%
      rownames_to_column("kinase") %>%
      pivot_longer(!kinase, names_to = "sample", values_to = "score") %>%
      mutate(id = paste(kinase, sample, sep = ":")) %>%
      dplyr::select(id, score)

    colnames(score_long)[2] <- method_idx
    score_long
  }) %>%
    purrr::reduce(full_join, by = "id") %>%
    column_to_rownames("id")

  cor(scores_all, use = "complete.obs")
})

names(prior_comparison) <- names(act_list[[1]])


Y <- do.call(cbind, prior_comparison)
Y_array <- array(Y, dim=c(dim(prior_comparison[[1]]), length(prior_comparison)))

mean_correlation <- colMeans(aperm(Y_array, c(3, 1, 2)), na.rm = TRUE)
colnames(mean_correlation) <- colnames(prior_comparison[[1]])
rownames(mean_correlation) <- rownames(prior_comparison[[1]])

pdf("results/manuscript_figures/figure_2/corrplot_priors.pdf", height = 3, width = 3.5)
corrplot(mean_correlation, col.lim=c(0, 1), is.corr = FALSE, type = 'lower', tl.col = 'black', order = 'hclust',  cl.pos = 'r')
dev.off()
