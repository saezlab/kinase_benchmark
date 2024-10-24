#'

## Snakemake ---------------------------
if(exists("snakemake")){
  tumor_files <- snakemake@input$bench
  auroc_plot <- snakemake@output$plot
}else{
  tumor_files <- "data/tumor_benchmark/activity_scores/roc_data.rds"
  auroc_plot <- "results/manuscript_figures/figure_2/auroc_tumor_phosphositeplus.pdf"
}

## Libraries ---------------------------
library(tidyverse)

## Load AUROC tumor benchmark ---------------------------
roc <- readRDS(tumor_files)
psp_roc <- roc$phosphositeplus
roc_list <- lapply(psp_roc$ROC_results, "[[", "sample_AUROCs")

df_tumor <- do.call(rbind, lapply(names(roc_list), function(method) {
  data.frame(method = method, score = roc_list[[method]])
}))

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
  geom_boxplot(fill = "#b54d4a", linewidth = 0.5, outlier.size = 0.8) + # Boxplot for AUROC
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


