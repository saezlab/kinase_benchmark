if(exists("snakemake")){
  rank_files <- snakemake@input$rank
  activating_files <- snakemake@input$act
  tumor_files <- snakemake@input$tumor
  performance_plot <- snakemake@output$plt
}else{
  rank_files <- list.files("results/03_benchmark/merged/02_mean_rank",
                           pattern = "csv", recursive = TRUE, full.names = T)
  rank_files <- rank_files[str_detect(rank_files, "/shuffled2/|/iKiPdb/|/GPS/|/omnipath/|/networkin/|/phosphositeplus/|/ptmsigdb/|/combined/|/GSknown/")]
  activating_strat <- "data/tumor_benchmark/performance_regulon/known_strat_rocs_actsiteBM.Rds"
  activating_kins <- "data/tumor_benchmark/performance_regulon/known_eval_kins_actsiteBM.Rds"
  performance_pert <- "results/manuscript_figures/figure_3/regulon_perturbation.pdf"
  performance_act <- "results/manuscript_figures/figure_3/regulon_activatingsites.pdf"
  performance_tumor <- "results/manuscript_figures/figure_3/regulon_tumor.pdf"
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

rank_df <- bind_rows(rank_list)

# Stratify kinases
kinase_strat <- rank_df %>%
  filter(prior == "GSknown") %>%
  group_by(targets) %>%
  summarise(measured_targets = mean(measured_targets)) %>%
  dplyr::mutate(strat = case_when(
   measured_targets < 11 ~ "small",
   measured_targets > 10 & measured_targets < 26 ~ "medium",
   measured_targets > 25 ~ "large"
  )) %>%
  select(targets, strat) %>%
  dplyr::filter(!is.na(strat))

rank_df <- rank_df %>%
  dplyr::select(method, prior, rank, scaled_rank, sample, targets) %>%
  filter(!is.na(rank)) %>%
  mutate(prior = recode(prior,
                        "iKiPdb" = "iKiP-DB",
                        "GPS" = "GPS gold",
                        "omnipath" = "OmniPath",
                        "networkin" = "NetworKIN",
                        "phosphositeplus" = "PhosphoSitePlus",
                        "ptmsigdb" = "PTMsigDB",
                        "combined" = "extended combined",
                        "GSknown" = "curated combined",
                        "shuffled2" = "Shuffled")) %>%
  left_join(kinase_strat, by = "targets") %>%
  dplyr::filter(!is.na(strat))

order_net <- c("curated combined", "GPS gold", "PhosphoSitePlus", "PTMsigDB", "NetworKIN", "OmniPath", "extended combined", "iKiP-DB", "Shuffled")
rank_df$prior <- factor(rank_df$prior, levels = order_net)
rank_df$strat <- factor(rank_df$strat, levels = c("small", "medium", "large"))

n_kinases <- rank_df %>%
  group_by(prior, method, strat) %>%
  summarise(kinases = length(unique(targets))) %>%
  filter(method == "zscore") %>%
  select(prior, kinases, strat) %>%
  add_column(benchmark = "perturbation-based")
n_kinases$prior <- factor(n_kinases$prior, levels = order_net)
n_kinases$strat <- factor(n_kinases$strat, levels = c("small", "medium", "large"))


auroc_p <- ggplot(rank_df, aes(x = prior, y = scaled_rank, fill = strat)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.1) +
  scale_fill_manual(values = c("#96c3df", "#4292C6", "#265c7f")) +  # Muted scientific color palette
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


kin_p <- ggplot(n_kinases, aes(x = prior, y = kinases, fill = strat)) +
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
  scale_fill_manual(values = c("#96c3df", "#4292C6", "#265c7f")) +
  ggtitle("Number of Kinases in Evaluation Set")

full_p <- ggarrange(kin_p, auroc_p, ncol = 1, common.legend = T, heights = c(3.1, 9))

pdf(performance_pert, width = 3.9, height = 3.5)
full_p
dev.off()


## tumor benchmark
df_act <- readRDS(activating_strat) %>%
  as.data.frame() %>%
  add_column(benchmark = "tumor-based") %>%
  pivot_longer(!benchmark, names_to = "prior", values_to = "auroc") %>%
  mutate(strat = map_chr(str_split(prior, "_"), 2)) %>%
  mutate(strat = case_when(
    strat == "0to10" ~ "small",
    strat == "11to15" ~ "medium",
    strat == "26plus" ~ "large"
  )) %>%
  mutate(prior = map_chr(str_split(prior, "_"), 1)) %>%
  mutate(prior = recode(prior,
                      "ikip" = "iKiP-DB",
                      "gps" = "GPS gold",
                      "omni" = "OmniPath",
                      "nwkin" = "NetworKIN",
                      "psp" = "PhosphoSitePlus",
                      "ptmsig" = "PTMsigDB",
                      "combo" = "extended combined",
                      "known" = "curated combined",
                      "shuffled2" = "shuffled"))


df_act$prior <- factor(df_act$prior, levels = order_net)
df_act$strat <- factor(df_act$strat, levels = c("small", "medium", "large"))

auroc_p <- ggplot(df_act, aes(x = prior, y = auroc, fill = strat)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.1) +
  scale_fill_manual(values = c("#e2b99e", "#C67642", "#58331a")) +  # Muted scientific color palette
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
  ylab("AUROC") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5)

kin_p <- ggplot(n_kinases, aes(x = prior, y = kinases, fill = strat)) +
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
  scale_fill_manual(values = c("#e2b99e", "#C67642", "#58331a")) +
  ggtitle("Number of Kinases in Evaluation Set")

full_p <- ggarrange(kin_p, auroc_p, ncol = 1, common.legend = T, heights = c(3.1, 9))

pdf(performance_act, width = 3.9, height = 3.5)
full_p
dev.off()


## Tumor
auroc_p <- ggplot(df_act, aes(x = prior, y = auroc, fill = strat)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.1) +
  scale_fill_manual(values = c("#d9a4a2", "#b54d4a", "#572523")) +  # Muted scientific color palette
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
  ylab("AUROC") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5)

kin_p <- ggplot(n_kinases, aes(x = prior, y = kinases, fill = strat)) +
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
  scale_fill_manual(values = c("#d9a4a2", "#b54d4a", "#572523")) +
  ggtitle("Number of Kinases in Evaluation Set")

full_p <- ggarrange(kin_p, auroc_p, ncol = 1, common.legend = T, heights = c(3.1, 9))

pdf(performance_tumor, width = 3.9, height = 3.5)
full_p
dev.off()
