if(exists("snakemake")){
  rank_files <- snakemake@input$rank
  kinase_class_file <- snakemake@input$kinclass
  performance_pert_ser <- snakemake@output$plt_ser
  performance_pert_tyr <- snakemake@output$plt_tyr
}else{
  rank_files <- list.files("results/03_benchmark/merged2/02_mean_rank",
                           pattern = "csv", recursive = TRUE, full.names = T)
  rank_files <- rank_files[str_detect(rank_files, "/iKiPdb/|/GPS/|/omnipath/|/networkin/|/phosphositeplus/|/ptmsigdb/|/GSknown/|/shuffled2/")]
  kinase_class_file <- "resources/kinase_class.csv"
  performance_pert_ser <- "results/manuscript_figures/figure_3/supp/ser_performance.pdf"
  performance_pert_tyr <- "results/manuscript_figures/figure_3/supp/tyr_performance.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(ggpubr)
library(patchwork)

## Benchmark ------------------
## Load AUROC ---------------------------
rank_list <- map(rank_files, function(file){
  read_csv(file, col_types = cols())
})

rank_df <- bind_rows(rank_list) %>%
  dplyr::select(method, prior, rank, scaled_rank, sample, targets) %>%
  filter(!is.na(rank)) %>%
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
                         "zscore" = "z-score")) %>%
  mutate(prior = recode(prior,
                        "iKiPdb" = "iKiP-DB",
                        "GPS" = "GPS gold",
                        "omnipath" = "OmniPath",
                        "networkin" = "NetworKIN",
                        "phosphositeplus" = "PhosphoSitePlus",
                        "ptmsigdb" = "PTMsigDB",
                        "combined" = "extended combined",
                        "GSknown" = "curated",
                        "shuffled2" = "shuffled")) %>%
    filter(method == "z-score")

kinase_class <- read_csv(kinase_class_file, col_types = cols())

rank_df_class <- rank_df %>%
    left_join(kinase_class %>% dplyr::select(source, class) %>% distinct(), by = c("targets" = "source"), relationship = "many-to-many") %>% 
    filter(class %in% c("Tyrosine", "Serine/Threonine"))

rank_df_class_ser <- rank_df_class %>%
    filter(class == "Serine/Threonine")

rank_df_class_tyr <- rank_df_class %>%
    filter(class == "Tyrosine")

order_ser <- rank_df_class_ser %>% 
    group_by(prior) %>%
    summarise(scaled_rank = mean(scaled_rank)) %>%
    arrange(scaled_rank)

order_ser
order_ser %>% filter(!prior == "shuffled") %>% pull(scaled_rank) %>% mean()

order_tyr <- rank_df_class_tyr %>% 
    group_by(prior) %>%
    summarise(scaled_rank = mean(scaled_rank)) %>%
    arrange(scaled_rank)

order_tyr
order_tyr %>% filter(!prior == "shuffled") %>% pull(scaled_rank) %>% mean()

full_df <- left_join(order_ser, order_tyr, by = "prior")
cor.test(full_df$scaled_rank.x, full_df$scaled_rank.y)

rank_df_class_ser$prior <- factor(rank_df_class_ser$prior, levels = order_ser$prior)
rank_df_class_tyr$prior <- factor(rank_df_class_tyr$prior, levels = order_tyr$prior)

n_kinases_ser <- rank_df_class_ser %>%
  group_by(prior, method) %>%
  summarise(kinases = length(unique(targets))) %>%
  select(prior, kinases)
n_kinases_ser$prior <- factor(n_kinases_ser$prior, levels = order_ser$prior)

n_kinases_ser

n_kinases_tyr <- rank_df_class_tyr %>%
  group_by(prior, method) %>%
  summarise(kinases = length(unique(targets))) %>%
  select(prior, kinases)
n_kinases_tyr$prior <- factor(n_kinases_tyr$prior, levels = order_tyr$prior)

n_kinases_tyr


auroc_p_ser <- ggplot(rank_df_class_ser, aes(x = prior, y = scaled_rank)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.1) +# Muted scientific color palette
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0, "lines"),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 11),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("") +
  ylab("scaled rank")


kin_p_ser <- ggplot(n_kinases_ser, aes(x = prior, y = kinases)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.4)+ # Line connecting the dots
  scale_y_continuous(
    name = "tmp",
    limits = c(0, 40),
    breaks = seq(0, 40, by = 10) #
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(0, "lines"),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_blank(),   # Remove x-axis text labels
    axis.ticks.x = element_blank()
  )+
  ggtitle("Number of Kinases in Evaluation Set")

full_p_ser <- ggarrange(kin_p_ser, auroc_p_ser, ncol = 1, common.legend = T, heights = c(3.1, 9))

pdf(performance_pert_ser, width = 3.9, height = 3.5)
full_p_ser
dev.off()

write_csv(rank_df_class_ser %>% dplyr::select(prior, targets, scaled_rank), "results/manuscript_data/suppfig10a_boxplot.csv")
write_csv(n_kinases_ser, "results/manuscript_data/suppfig10a_barplot.csv")

auroc_p_tyr <- ggplot(rank_df_class_tyr, aes(x = prior, y = scaled_rank)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.1) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +# Muted scientific color palette
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0, "lines"),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 11),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("") +
  ylab("scaled rank")


kin_p_tyr <- ggplot(n_kinases_tyr, aes(x = prior, y = kinases)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.4)+ # Line connecting the dots
  scale_y_continuous(
    name = "tmp",
    limits = c(0, 15),
    breaks = seq(0, 15, by = 5) #
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(0, "lines"),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_blank(),   # Remove x-axis text labels
    axis.ticks.x = element_blank()
  )+
  ggtitle("Number of Kinases in Evaluation Set")

full_p_tyr <- ggarrange(kin_p_tyr, auroc_p_tyr, ncol = 1, common.legend = T, heights = c(3.1, 9))

pdf(performance_pert_tyr, width = 3.9, height = 3.5)
full_p_tyr
dev.off()

write_csv(rank_df_class_tyr %>% dplyr::select(prior, targets, scaled_rank), "results/manuscript_data/suppfig10b_boxplot.csv")
write_csv(n_kinases_tyr, "results/manuscript_data/suppfig10b_barplot.csv")
