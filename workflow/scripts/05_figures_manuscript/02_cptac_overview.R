#'

## Snakemake ---------------------------
if(exists("snakemake")){
  tumor_files <- snakemake@input$bench
  act_files <- snakemake@input$act
  norm_act_files <- snakemake@input$norm_act
  norm_prot_files <- snakemake@input$prot_act
  auroc_plot <- snakemake@output$plot
  auroc_plot_act <- snakemake@output$plotact
  norm_plot <- snakemake@output$norm_plot
  norm_plot_act <- snakemake@output$normact_plot
}else{
  act_files <- "data/results_cptac/overall_performance/actsiteBM/all_kins/actsiteBM_5perThr_psp_roc_table.rds"
  tumor_files <- "data/results_cptac/overall_performance/protBM/all_kins/protBM_5perThr_psp_roc_table.rds"
  norm_act_files <- list.files("data/results_cptac/normalisation/actsiteBM", full.names = T)
  norm_prot_files <- list.files("data/results_cptac/normalisation/proteinBM", full.names = T)
  auroc_plot <- "results/manuscript_figures/figure_2/auroc_tumor_phosphositeplus.pdf"
  auroc_plot_act <- "results/manuscript_figures/figure_2/auroc_act_phosphositeplus.pdf"
  norm_plot <- "results/manuscript_figures/figure_2/norm_tumor_phosphositeplus.pdf"
  norm_plot_act <- "results/manuscript_figures/figure_2/norm_act_phosphositeplus.pdf"
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

## Normalisation ---------------------------
df_norm_prot <- map_dfr(norm_prot_files, function(file_roc){
    roc <- readRDS(file_roc)
    net_id <- str_extract(file_roc, "(?<=5perThr_).*?(?=_roc_table)")
    map_dfr(colnames(roc), function(method_id){
        data.frame(normalisation = net_id, score = roc[,colnames(roc) == method_id], method = method_id, benchmark = "tumor bench")
    })
}) %>%
   mutate(normalisation = recode(normalisation,
                        "glob" = "global LM",
                        "perprot" = "per prot LM",
                        "persite" = "per site LM",
                        "sub" = "ratio/subtract",
                        "unnorm" = "no norm")) %>%
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
  filter(method %in% c("z-score", "PTM-SEA", "mean", "fgsea", "VIPER", "KSEA", "n targets")) 

order_norm <- df_norm_prot %>%
  group_by(normalisation) %>%
  summarise(score = mean(score)) %>%
  arrange(desc(score))

order_method <- df_norm_prot %>%
  group_by(method) %>%
  summarise(score = mean(score)) %>%
  arrange(desc(score))

df_norm_prot$normalisation <- factor(df_norm_prot$normalisation, levels = order_norm$normalisation)
df_norm_prot$method <- factor(df_norm_prot$method, levels = order_method$method)

norm_prot_plt <- df_norm_prot %>%
  ggplot( aes(x = method, y = score, fill = normalisation)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.1) + # Boxplot for AUROC
  scale_fill_manual(values = c("#6B86C3", "#6D5A79", "#B66577", "#E76B6E", "#ECAD8B")) +
  scale_y_continuous(
    name = "AUROC") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5) +
  theme_bw() +
  theme(
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 11),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("")

pdf(norm_plot, width = 4, height = 3.5)
norm_prot_plt
dev.off()



df_norm_act <- map_dfr(norm_act_files, function(file_roc){
    roc <- readRDS(file_roc)
    net_id <- str_extract(file_roc, "(?<=5perThr_).*?(?=_roc_table)")
    map_dfr(colnames(roc), function(method_id){
        data.frame(normalisation = net_id, score = roc[,colnames(roc) == method_id], method = method_id, benchmark = "activating sites")
    })
}) %>%
   mutate(normalisation = recode(normalisation,
                        "glob" = "global LM",
                        "perprot" = "per prot LM",
                        "persite" = "per site LM",
                        "sub" = "ratio/subtract",
                        "unnorm" = "no norm")) %>%
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
  filter(method %in% c("z-score", "PTM-SEA", "mean", "fgsea", "VIPER", "KSEA", "n targets")) 

order_norm <- df_norm_act %>%
  group_by(normalisation) %>%
  summarise(score = mean(score)) %>%
  arrange(desc(score))

order_method <- df_norm_act %>%
  group_by(method) %>%
  summarise(score = mean(score)) %>%
  arrange(desc(score))

df_norm_act$normalisation <- factor(df_norm_act$normalisation, levels = order_norm$normalisation)
df_norm_act$method <- factor(df_norm_act$method, levels = order_method$method)

norm_act_plt <- df_norm_act %>%
  ggplot( aes(x = method, y = score, fill = normalisation)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.1) + # Boxplot for AUROC
  scale_fill_manual(values = c("#6B86C3", "#6D5A79", "#B66577", "#E76B6E", "#ECAD8B")) +
  scale_y_continuous(
    name = "AUROC") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5) +
  theme_bw() +
  theme(
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 11),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("")

pdf(norm_plot_act, width = 4, height = 3.5)
norm_act_plt
dev.off()
