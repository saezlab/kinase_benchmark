#'

## Snakemake ---------------------------
if(exists("snakemake")){
  tumor_files <- snakemake@input$bench
  act_files <- snakemake@input$act
  auroc_plot <- snakemake@output$plot
  auroc_plot_act <- snakemake@output$plotact
}else{
  act_files <- "data/results_cptac/overall_performance/actsiteBM/all_kins/actsiteBM_5perThr_psp_roc_table.rds"
  tumor_files <- "data/results_cptac/overall_performance/protBM/all_kins/protBM_5perThr_psp_roc_table.rds"
  auroc_plot <- "results/manuscript_figures/figure_2/auroc_tumor_phosphositeplus.pdf"
  auroc_plot_act <- "results/manuscript_figures/figure_2/auroc_act_phosphositeplus.pdf"
}

## Libraries ---------------------------
library(tidyverse)

## Load AUROC tumor benchmark ---------------------------
psp_roc <- readRDS(tumor_files)
df_tumor <- psp_roc %>%
  as.data.frame() %>%
  add_column(net = "phosphositeplus") %>%
  pivot_longer(!net, values_to = "score", names_to = "method")

df_tumor <- df_tumor %>%
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


# Sorting methods by mean AUROC
mean_auroc <- df_tumor %>%
  group_by(method) %>%
  summarize(mean_auroc = mean(score)) %>%
  arrange(desc(mean_auroc))

# Update df to ensure methods are ordered based on mean AUROC
df_tumor <- df_tumor %>%
  mutate(method = factor(method, levels = mean_auroc$method))

## Plot performance ---------------------------
auroc_p <- ggplot(df_tumor, aes(x = method, y = score)) +
  geom_boxplot(fill = "#b54d4a", linewidth = 0.3, outlier.size = 0.1) + # Boxplot for AUROC
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

pdf(auroc_plot, width = 4, height = 2.5)
auroc_p
dev.off()


## Load AUROC tumor benchmark ---------------------------
psp_roc <- readRDS(act_files)
df_tumor <- psp_roc %>%
  as.data.frame() %>%
  add_column(net = "phosphositeplus") %>%
  pivot_longer(!net, values_to = "score", names_to = "method")


df_tumor <- df_tumor %>%
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


# Sorting methods by mean AUROC
mean_auroc <- df_tumor %>%
  group_by(method) %>%
  summarize(mean_auroc = mean(score)) %>%
  arrange(desc(mean_auroc))

# Update df to ensure methods are ordered based on mean AUROC
df_tumor <- df_tumor %>%
  mutate(method = factor(method, levels = mean_auroc$method))

## Plot performance ---------------------------
auroc_p <- ggplot(df_tumor, aes(x = method, y = score)) +
  geom_boxplot(fill = "#C67642", linewidth = 0.3, outlier.size = 0.1) + # Boxplot for AUROC
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

pdf(auroc_plot_act, width = 4, height = 2.5)
auroc_p
dev.off()


