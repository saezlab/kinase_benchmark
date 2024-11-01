if(exists("snakemake")){
  act_files <- snakemake@input$act_files
  jaccard_i <-  snakemake@params$jaccard_i
  jacc_type <-  snakemake@params$jacc_type
  jaccard_methods <- snakemake@output$jaccMethods
  jaccard_priors <-  snakemake@output$jaccPriors
  plot_methods <- snakemake@output$plotMethods
  plot_priors <-  snakemake@output$plotPriors
}else{
  act_files <- list.files("results/02_activity_scores/merged/scores", pattern = "rds", recursive = T, full.names = T)
  act_files <- act_files[c(8, 4, 10, 9, 13, 14, 3, 7)]
  jaccard_i <- 10
  jacc_type <- "up"
  jaccard_methods <- "results/04_exploration/merged/jaccard/jaccard_methods_up_10.rds"
  jaccard_priors <- "results/04_exploration/merged/jaccard/jaccard_priors_up_10.rds"
  plot_methods <- "results/04_exploration/merged/jaccard/plots/jaccard_methods_up_10.pdf"
  plot_priors <- "results/04_exploration/merged/jaccard/plots/jaccard_priors_up_10.pdf"
}
jaccard_i <- as.numeric(jaccard_i)

## Libraries ---------------------------
library(tidyverse)
library(corrplot)

## Load and merge scores ---------------------------
act_list <- map(act_files, readRDS)
names(act_list) <- map_chr(str_split(act_files, "/"), 5) %>% str_remove(".rds")

methods <- names(act_list[[1]])

## Jaccard index ---------------------------
jaccard_index <- function(set1, set2) {
  length(intersect(set1, set2)) / length(union(set1, set2))
}

# Jaccard between methods for each prior
method_jaccard <- map(names(act_list), function(prior_idx){
  scores <- act_list[[prior_idx]]
  scores <- scores[names(scores) %in% methods]

  kin_input <- map_dfr(names(scores), function(method_idx){
    score_long <- scores[[method_idx]] %>%
      rownames_to_column("kinase") %>%
      pivot_longer(!kinase, names_to = "sample", values_to = "score") %>%
      filter(!is.na(score)) %>%
      group_by(sample)
    if (jacc_type == "up"){
        score_long <- score_long %>%
          arrange(desc(score))
    } else if (jacc_type == "down"){
        score_long <- score_long %>%
          arrange(score)
    } else if (jacc_type == "abs"){
        score_long <- score_long %>%
          arrange(desc(abs(score)))
    }
    score_long <- score_long %>%
      slice(1:jaccard_i) %>%
      add_column(method = method_idx) %>%
      ungroup()
  })

  # Calculate Jaccard Index for each sample
  jaccard_per_sample <- kin_input %>%
    group_by(sample) %>%
    summarise(
        jaccard_values = list({
        # Find all unique methods
        methods <- unique(kin_input$method)
        # Generate all combinations of method pairs
        combs <- expand.grid(methods, methods, stringsAsFactors = FALSE)
        # Calculate Jaccard index for each combination of methods
        sapply(1:nrow(combs), function(i) {
            set1 <- kinase[method == combs$Var1[i]]
            set2 <- kinase[method == combs$Var2[i]]
            jaccard_index(set1, set2)
        })
        }),
        method_pairs = list(expand.grid(unique(method), unique(method), stringsAsFactors = FALSE))
    ) %>%
    unnest(c(jaccard_values, method_pairs)) %>%
    ungroup()

  mean_jaccard_per_pair <- jaccard_per_sample %>%
    group_by(Var1, Var2) %>%
    summarise(mean_jaccard = mean(jaccard_values)) %>%
    ungroup()

  mean_jaccard_per_pair %>%
  pivot_wider(names_from = Var2, values_from = mean_jaccard) %>%
  column_to_rownames("Var1")
})

names(method_jaccard) <- names(act_list)


prior_jaccard <- map(names(act_list[[1]]), function(method_idx){
  scores <- map(act_list, function(x) x[[method_idx]])

  kin_input <- map_dfr(names(scores), function(prior_idx){
    score_long <- scores[[prior_idx]] %>%
      rownames_to_column("kinase") %>%
      pivot_longer(!kinase, names_to = "sample", values_to = "score") %>%
      filter(!is.na(score)) %>%
      group_by(sample)
    if (jacc_type == "up"){
        score_long <- score_long %>%
          arrange(desc(score))
    } else if (jacc_type == "down"){
        score_long <- score_long %>%
          arrange(score)
    } else if (jacc_type == "abs"){
        score_long <- score_long %>%
          arrange(desc(abs(score)))
    }
    score_long <- score_long %>%
      slice(1:jaccard_i) %>%
      add_column(method = prior_idx) %>%
      ungroup()
  })

  # Calculate Jaccard Index for each sample
  jaccard_per_sample <- kin_input %>%
    group_by(sample) %>%
    summarise(
        jaccard_values = list({
        # Find all unique methods
        methods <- unique(kin_input$method)
        # Generate all combinations of method pairs
        combs <- expand.grid(methods, methods, stringsAsFactors = FALSE)
        # Calculate Jaccard index for each combination of methods
        sapply(1:nrow(combs), function(i) {
            set1 <- kinase[method == combs$Var1[i]]
            set2 <- kinase[method == combs$Var2[i]]
            jaccard_index(set1, set2)
        })
        }),
        method_pairs = list(expand.grid(unique(method), unique(method), stringsAsFactors = FALSE))
    ) %>%
    unnest(c(jaccard_values, method_pairs)) %>%
    ungroup()

  mean_jaccard_per_pair <- jaccard_per_sample %>%
    group_by(Var1, Var2) %>%
    summarise(mean_jaccard = mean(jaccard_values)) %>%
    ungroup()

  mean_jaccard_per_pair %>%
    pivot_wider(names_from = Var2, values_from = mean_jaccard) %>%
    column_to_rownames("Var1")
})

names(prior_jaccard) <- names(act_list[[1]])

## Plot jaccard ---------------------------
method_mean <- map_dfr(names(method_jaccard), function(prior_id){
    mat <- method_jaccard[[prior_id]]
    mat %>%
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

prior_mean <- map_dfr(names(prior_jaccard), function(prior_id){
    mat <- prior_jaccard[[prior_id]]
    mat %>%
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
saveRDS(method_jaccard, jaccard_methods)
saveRDS(prior_jaccard, jaccard_priors)
