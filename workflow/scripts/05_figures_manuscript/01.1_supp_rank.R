#'

## Snakemake ---------------------------
if(exists("snakemake")){
  rank_files <- snakemake@input$rank
  prior_id <- snakemake@wildcards$PKN
  method_id <- snakemake@wildcards$methods
  performance_plot <- snakemake@output$plot
}else{
  rank_files <-  "results/03_benchmark/merged/02_mean_rank/phosphositeplus/zscore-phosphositeplus.csv"
  performance_plot <- "results/manuscript_figures/figure_1/supp/comparison_rank_subset.pdf"
  prior_id <- "phosphositeplus"
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

tmp <- map_dfr(kinase_pert, function(kin_id){
  exp_target <- rank_df %>%
    dplyr::filter(targets == kin_id) %>%
    pull(sample)
  class <- rank_df %>%
    dplyr::filter(targets == kin_id) %>%
    pull(class) %>%
    unique()
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
             sd_nonpert = sd(nonperturbed_rank, na.rm = T),
             class = class)
})

rank_overview <- tmp %>%
  dplyr::select(-class) %>%
  pivot_longer(!kinase, names_to = "type", values_to = "scaled_rank") %>%
  mutate(benchmark = map_chr(str_split(type, "_"), 2)) %>%
  mutate(metric = map_chr(str_split(type, "_"), 1)) %>%
  left_join(tmp %>% dplyr::select(kinase,class), by = "kinase")


## sort kinases
kin_order <- rank_overview %>%
  filter(metric == "mean") %>%
  dplyr::select(kinase, scaled_rank, benchmark) %>%
  pivot_wider(names_from = benchmark, values_from = scaled_rank) %>%
  group_by(kinase) %>%
  summarise(diff = nonpert/pert) %>%
  arrange(diff)

rank_overview$kinase <- factor(rank_overview$kinase, levels = kin_order$kinase)

scaled_p <- ggplot(rank_overview %>% filter(metric == "mean"), aes(x = scaled_rank, y = kinase, color = benchmark)) +
  geom_point() +
  scale_color_manual(values = c("darkgrey", "#4292C6")) +  # Muted scientific color palette
  theme_bw() +
  theme(
    panel.spacing.x = unit(0, "lines"),
    text = element_text(family = "Helvetica", size = 10), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  )
  #
rank_overview %>% filter(class == "Tyrosine")
rank_overview %>% filter(class == "Dual-specificity")

tmp$mean_pert %>% mean()
tmp$mean_nonpert %>% mean()

tmp$mean_pert %>% sd()
tmp$mean_nonpert %>% sd()

tmp$mean_pert %>% range()
tmp$mean_nonpert %>% range()

pdf(performance_plot, width = 5, height = 4.8)
scaled_p
dev.off()
