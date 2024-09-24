if(exists("snakemake")){
  act_files <- snakemake@input$act_files
  corr_type <-  snakemake@params$corr_type
  correlation_methods <- snakemake@output$corrMethods
  correlation_priors <-  snakemake@output$corrPriors
  plot_methods <- snakemake@output$plotMethods
  plot_priors <-  snakemake@output$plotPriors
}else{
  act_files <- list.files("results/02_activity_scores/merged/scores", pattern = "rds", recursive = T, full.names = T)
  corr_type <- "spearman"
  correlation_methods <- "results/04_exploration/merged/correlation/correlation_methods_spearman.rds"
  correlation_priors <- "results/04_exploration/merged/correlation/correlation_priors_spearman.rds"
  plot_methods <- "results/04_exploration/merged/correlation/plots/correlation_methods_spearman.rds"
  plot_priors <- "results/04_exploration/merged/correlation/plots/correlation_priors_spearman.rds"
}

## Libraries ---------------------------
library(tidyverse)
library(corrplot)

## Load and merge scores ---------------------------
act_list <- map(act_files, readRDS)
names(act_list) <- map_chr(str_split(act_files, "/"), 5) %>% str_remove(".rds")

## Correlation ---------------------------
# Correlation between methods for each prior
methods <- names(act_list[[1]])

method_comparison <- map(names(act_list), function(prior_idx){
  scores <- act_list[[prior_idx]]
  scores <- scores[names(scores) %in% methods]

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

  cor(scores_all, method = corr_type, use = "complete.obs")
})

names(method_comparison) <- names(act_list)

# Correlation between priors for each method
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

  cor(scores_all, method = corr_type, use = "complete.obs")
})

names(prior_comparison) <- names(act_list[[1]])
 
## Plot jaccard ---------------------------
method_mean <- map_dfr(names(method_comparison), function(prior_id){
    mat <- method_comparison[[prior_id]]
    mat %>% 
    data.frame() %>%
    rownames_to_column("method1") %>%
      pivot_longer(!method1, names_to = "method2", values_to = "jaccard")
}) %>%
  group_by(method1, method2) %>%
  summarise(jaccard = mean(jaccard)) %>%
  ungroup() %>%
  pivot_wider(names_from = method2, values_from = jaccard) %>%
  column_to_rownames("method1") 
  
method_mean <- method_mean %>%
  as.matrix()


pdf(plot_methods, height = 3, width = 3.5)
corrplot(method_mean,
         col.lim=c(min(method_mean), 1),
         is.corr = FALSE, type = 'lower',
         tl.col = 'black', order = 'hclust', cl.pos = 'r')
dev.off()

prior_mean <- map_dfr(names(prior_comparison), function(prior_id){
    mat <- prior_comparison[[prior_id]]
    mat %>% 
    data.frame() %>% 
    rownames_to_column("prior1") %>%
      pivot_longer(!prior1, names_to = "prior2", values_to = "jaccard")
}) %>%
  group_by(prior1, prior2) %>%
  summarise(jaccard = mean(jaccard)) %>%
  ungroup() %>%
  pivot_wider(names_from = prior2, values_from = jaccard) %>%
  column_to_rownames("prior1") 
  
prior_mean <- prior_mean %>%
  as.matrix()

pdf(plot_priors, height = 3, width = 3.5)
corrplot(prior_mean,
         col.lim=c(min(prior_mean), 1),
         is.corr = FALSE, type = 'lower',
         tl.col = 'black', order = 'hclust', cl.pos = 'r')
dev.off()

## Save files ---------------------------
saveRDS(method_comparison, correlation_methods)
saveRDS(prior_comparison, correlation_priors)