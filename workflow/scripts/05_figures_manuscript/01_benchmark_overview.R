#'

## Snakemake ---------------------------
if(exists("snakemake")){
  bench_files <- snakemake@input$bench
  rank_files <- snakemake@input$rank
  k_phit_c <- snakemake@params$k_phit
  auroc_plot <- snakemake@output$auroc
  rank_plot <- snakemake@output$rank
  p_hit_plot <- snakemake@output$phit
}else{
  bench_files <- list.files("results/03_benchmark/merged/02_benchmark_res/phosphositeplus",
                            pattern = "bench", recursive = TRUE, full.names = T)
  rank_files <-  list.files("results/03_benchmark/merged/02_mean_rank/phosphositeplus",
                            pattern = "csv", recursive = TRUE, full.names = T)
  k_phit_c <- c(5, 10, 20)
  auroc_plot <- "results/manuscript_figures/figure_1/auroc_phosphositeplus.pdf"
  rank_plot <- "results/manuscript_figures/figure_1/scaledrank_phosphositeplus.pdf"
  p_hit_plot <- "results/manuscript_figures/figure_1/phit_phosphositeplus.pdf"
}

## Libraries ---------------------------
library(tidyverse)

## Load scaled rank and pHit ---------------------------
rank_list <- map(rank_files, function(file){
  read_csv(file, col_types = cols())
})

rank_df <- bind_rows(rank_list) %>%
  dplyr::select(method, prior, rank, scaled_rank, sample, targets) %>%
  filter(!is.na(rank))

rank_df <- rank_df %>%
  mutate(method = recode(method,
                         "chisq" = "X\u00B2 test",
                         "fisher" = "Fisher",
                         "KS" = "KS test",
                         "lmRoKAI" = "lm RoKAI",
                         "norm_mean" = "norm mean",
                         "number_of_targets" = "n targets",
                         "ptmsea" = "PTM-SEA",
                         "viper" = "VIPER",
                         "wilcox" = "MWU test",
                         "zscore" = "z-score"))

# number of TP kinases
rank_df %>%
  group_by(prior, method) %>%
  summarise(kinases = length(unique(targets))) %>%
  pull(kinases) %>%
  unique()

pHit <- map_dfr(k_phit_c, function(k){rank_df %>%
  group_by(prior, method) %>%
  summarise(phit = (sum(rank <= k)/n())) %>%
  add_column(k_phit = k) %>%
  ungroup()
})

pHit %>% group_by(method) %>% summarise(phit = mean(phit)) %>% arrange(desc(phit))
pHit %>% filter(k_phit == 5) %>% arrange(desc(phit))
pHit %>% filter(k_phit == 10) %>% arrange(desc(phit))
pHit %>% filter(k_phit == 20) %>% arrange(desc(phit))

## Load AUROC ---------------------------
if (any(str_detect(bench_files, "subset"))){
  net_id <- 6
} else {
  net_id <- 5
}

bench_list <- map(bench_files, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][net_id])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]
bench_df <- bind_rows(bench_list) %>%
  dplyr::filter(metric == "mcauroc") %>%
  dplyr::rename("prior" = net) %>%
  mutate(method = recode(method,
                         "chisq" = "X\u00B2 test",
                         "fisher" = "Fisher",
                         "KS" = "KS test",
                         "lmRoKAI" = "lm RoKAI",
                         "norm_mean" = "norm mean",
                         "number_of_targets" = "n targets",
                         "ptmsea" = "PTM-SEA",
                         "viper" = "VIPER",
                         "wilcox" = "MWU test",
                         "zscore" = "z-score"))

## Plot AUROC ---------------------------
# Sorting methods by mean AUROC
mean_auroc <- bench_df %>%
  group_by(method) %>%
  summarize(mean_auroc = mean(score)) %>%
  arrange(desc(mean_auroc))

mean_auroc

mean_rank <- rank_df %>%
  group_by(method) %>%
  summarize(mean_auroc = mean(scaled_rank)) %>%
  arrange(mean_auroc)

mean_rank

mean_phit <- pHit %>%
  group_by(method) %>%
  summarize(mean_auroc = mean(phit)) %>%
  arrange(desc(mean_auroc))

mean_phit

# Update df to ensure methods are ordered based on mean AUROC
bench_df <- bench_df %>%
  mutate(method = factor(method, levels = mean_auroc$method))

rank_df <- rank_df %>%
  mutate(method = factor(method, levels = mean_rank$method))

pHit <- pHit %>%
  mutate(method = factor(method, levels = mean_phit$method))

pHit <- pHit %>%
  mutate(k_phit = factor(k_phit, levels = rev(k_phit_c)))

## Plot performance ---------------------------
auroc_p <- ggplot(bench_df, aes(x = method, y = score)) +
  geom_boxplot(fill = "#4292C6", linewidth = 0.5, outlier.size = 0.8) + # Boxplot for AUROC
  scale_y_continuous(
    name = "AUROC") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("")

rank_p <- ggplot(rank_df, aes(x = method, y = scaled_rank)) +
  geom_boxplot(fill = "#4292C6", linewidth = 0.5, outlier.size = 0.8) + # Boxplot for AUROC
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("") +
  ylab("scaled rank")

p_hit_p <- ggplot(pHit, aes(x = method, y = phit, fill = k_phit)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.8), width = 0.7) + # Boxplot for AUROC
  theme_bw() +
  scale_fill_manual(values = c("5" = "#9ECAE1", "10" = "#4292C6", "20" = "#08519C"))+  # Use 'Blues' palette
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) + guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7)) +
  labs(fill = "k")+
  xlab("") +
  ylab(bquote(P[Hit](k)))


pdf(auroc_plot, width = 4, height = 2.5)
auroc_p
dev.off()

pdf(rank_plot, width = 4, height = 2.5)
rank_p
dev.off()

pdf(p_hit_plot, width = 8, height = 2.5)
p_hit_p
dev.off()

## Compare metrices ---------------------------
performance <- left_join(bench_df, rank_df, by = c("prior", "method")) %>%
  left_join(pHit, by = c("prior", "method")) %>%
  group_by(method, prior) %>%
  summarise(auroc = mean(score), scaled_rank = mean(scaled_rank), phit = mean(phit))

cor(performance[, 3:ncol(performance)], method = "pearson")
cor.test(performance$scaled_rank, performance$phit, method = "pearson")
cor.test(performance$scaled_rank, performance$auroc, method = "pearson")
cor.test(performance$phit, performance$auroc, method = "pearson")

