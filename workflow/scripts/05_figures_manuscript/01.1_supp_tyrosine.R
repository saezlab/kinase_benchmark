if(exists("snakemake")){
  rank_files <- snakemake@input$rank
  kinase_class_file <- snakemake@input$kinclass
  out_plot <- snakemake@output$out
  phit_plot <- snakemake@output$phit
}else{
  rank_files <- list.files("results/03_benchmark/merged2/02_mean_rank/phosphositeplus",
                           pattern = "csv", recursive = TRUE, full.names = T)
  kinase_class_file <- "resources/kinase_class.csv"
  out_plot <- "results/manuscript_figures/figure_1/supp/rank_class.pdf"
  phit_plot <- "results/manuscript_figures/figure_1/supp/phit_class.pdf"
}


## Libraries ---------------------------
library(tidyverse)

## Benchmark ------------------
## Load rank ---------------------------
rank_list <- map(rank_files, function(file){
  read_csv(file, col_types = cols()) %>%
    add_column(bench = str_split(file, "\\/")[[1]][3])
})

rank_df <- bind_rows(rank_list) %>%
  dplyr::select(method, prior, rank, scaled_rank, sample, targets, bench) %>%
  filter(!is.na(rank)) %>%
  mutate(prior = recode(prior, "phosphositeplus" = "PhosphoSitePlus")) %>%
  group_by(method, prior, targets) %>% 
  summarise(scaled_rank = mean(scaled_rank), rank = mean(rank)) %>% 
  ungroup() %>%
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


kinase_class <- read_csv(kinase_class_file, col_types = cols())

rank_df_class <- rank_df %>%
    left_join(kinase_class %>% select(source, class) %>% distinct(), by = c("targets" = "source"), relationship = "many-to-many") %>% 
    filter(class %in% c("Tyrosine", "Serine/Threonine"))

mean_auroc_tyr <- rank_df_class %>%
  filter(class == "Tyrosine") %>%
  group_by(method) %>%
  summarize(mean_auroc = mean(scaled_rank)) %>%
  arrange(mean_auroc)  # Sorting methods by mean AUROC

rank_df_class %>%
  filter(class == "Tyrosine") %>% pull(targets) %>% unique() %>% length()

mean_auroc_tyr %>%
  filter(!method == "n targets") %>% 
  pull(mean_auroc) %>% 
  mean()

mean_auroc_ser <- rank_df_class %>%
  filter(class == "Serine/Threonine") %>%
  group_by(method) %>%
  summarize(mean_auroc = mean(scaled_rank)) %>%
  arrange(mean_auroc) 

rank_df_class %>%
  filter(class == "Serine/Threonine") %>% pull(targets) %>% unique() %>% length()

mean_auroc_ser %>%
  filter(!method == "n targets") %>% 
  pull(mean_auroc) %>% 
  mean()

full_perf <- left_join(mean_auroc_ser, mean_auroc_tyr, by = "method")
cor.test(full_perf$mean_auroc.x, full_perf$mean_auroc.y)
# Update df to ensure methods are ordered based on mean AUROC
rank_df_tyr <- rank_df_class %>%
  filter(class == "Tyrosine") %>%
  mutate(method = factor(method, levels = mean_auroc_tyr$method))

rank_df_ser <- rank_df_class %>%
  filter(class == "Serine/Threonine") %>%
  mutate(method = factor(method, levels = mean_auroc_ser$method))


ser_p <- ggplot(rank_df_ser, aes(x = method, y = scaled_rank)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.5) +
  scale_fill_manual(values = c("#4292C6", "#AA42C6")) + # Muted scientific color palette
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0, "lines"),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("") +
  ylab("scaled_rank") +
  ggtitle("Serine/Threonine")

tyr_p <- ggplot(rank_df_tyr, aes(x = method, y = scaled_rank)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.5) +
  scale_fill_manual(values = c("#4292C6", "#AA42C6")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +  # Muted scientific color palette
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0, "lines"),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("") +
  ylab("scaled_rank") +
  ggtitle("Tyrosine")

pdf(out_plot, width = 4, height = 3.8)
ser_p
tyr_p
dev.off()

write_csv(rank_df_ser %>% dplyr::select(method, targets, scaled_rank), "results/manuscript_data/suppfig2a.csv")
write_csv(rank_df_tyr %>% dplyr::select(method, targets, scaled_rank), "results/manuscript_data/suppfig2b.csv")

## pHit ---------------------------
pHit <- map_dfr(c(10), function(k){rank_df_class %>%
  group_by(prior, method, class) %>%
  summarise(phit = (sum(rank <= k)/n())) %>%
  add_column(k_phit = k) %>%
  ungroup()
})

mean_phit_ser <- pHit %>%
  filter(class == "Serine/Threonine") %>%
  arrange(desc(phit), desc(method))

mean_phit_tyr <- pHit %>%
  filter(class == "Tyrosine") %>%
  arrange(desc(phit))

pHit_tyr <- pHit %>%
  filter(class == "Tyrosine") %>%
  mutate(method = factor(method, levels = mean_phit_tyr$method))

pHit_ser <- pHit %>%
  filter(class == "Serine/Threonine") %>%
  mutate(method = factor(method, levels = mean_phit_ser$method))

ser_phit_p <- ggplot(pHit_ser, aes(x = method, y = phit)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.8), width = 0.7) + # Boxplot for AUROC
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) + guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7)) +
  labs(fill = "k")+
  xlab("") +
  ylab(bquote(P[Hit](k))) +
  ggtitle("Serine/Threonine")

tyr_phit_p <- ggplot(pHit_tyr, aes(x = method, y = phit)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.8), width = 0.7) + # Boxplot for AUROC
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) + guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7)) +
  labs(fill = "k")+
  xlab("") +
  ylab(bquote(P[Hit](k))) +
  ggtitle("Tyrosine")


pdf(phit_plot, width = 4, height = 3.8)
ser_phit_p
tyr_phit_p
dev.off()
