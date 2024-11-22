#'

## Snakemake ---------------------------
if(exists("snakemake")){
  rank_files <- snakemake@input$rank
  prior_id <- snakemake@wildcards$PKN
  method_id <- snakemake@wildcards$hernandez_methods
  performance_plot <- snakemake@output$plot
}else{
  rank_files <-  "results/03_benchmark/merged/02_mean_rank/omnipath/zscore-omnipath.csv"                        
  performance_plot <- "results/04_exploration/merged/rank/comparison_rank_zscore_omnipath.pdf"
  prior_id <- "omnipath"
  method_id <- "zscore"
}

## Libraries ---------------------------
library(tidyverse)

## Load scaled rank ---------------------------
rank_df <- read_csv(rank_files, col_types = cols())

rank_df <- rank_df %>%
  dplyr::filter(prior == prior_id) %>%
  dplyr::filter(method == method_id)

kinase_pert <- rank_df %>%
  filter(!is.na(rank)) %>%
  pull(targets) %>%
  unique()

rank_overview <- map_dfr(kinase_pert, function(kin_id){
    exp_target <- rank_df %>%
      dplyr::filter(targets == kin_id) %>%
      pull(sample)
    perturbed_rank <- rank_df %>%
      dplyr::filter(targets == kin_id) %>%
      pull(scaled_rank)

    nonpert_df <- rank_df %>%
      dplyr::filter(!sample %in% exp_target)

     nonperturbed_rank <- map_dbl(1:nrow(nonpert_df), function(row_idx){
        kin <- nonpert_df$all_kinases_act[row_idx] %>%
          str_split(";") %>%
          unlist() 
        np_rank <- which(kin == kin_id)/length(kin)

        if (length(np_rank) == 0) {
            np_rank <- NA
            }
        np_rank
      })
    data.frame(kinase = kin_id, 
               mean_pert = mean(perturbed_rank, na.rm = T),
               sd_pert = sd(perturbed_rank, na.rm = T),
               mean_nonpert = mean(nonperturbed_rank, na.rm = T),
               sd_nonpert = sd(nonperturbed_rank, na.rm = T))
}) %>%
  pivot_longer(!kinase, names_to = "type", values_to = "scaled_rank") %>%
  mutate(benchmark = map_chr(str_split(type, "_"), 2)) %>%
  mutate(metric = map_chr(str_split(type, "_"), 1))
mean_non <- rank_overview %>% filter(benchmark == "nonpert") %>% filter(metric == "mean") %>% pull(scaled_rank) %>% mean(na.rm = T)

scaled_p <- ggplot(rank_overview %>% filter(metric == "mean"), aes(x = scaled_rank, y = kinase, color = benchmark)) +
 geom_point() +
 ggtitle(paste0(round(mean_non, digits = 2)))

pdf(performance_plot, width = 8, height = 4)
scaled_p
dev.off()
