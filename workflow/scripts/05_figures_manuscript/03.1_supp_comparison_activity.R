#' Process data for benchmark

## Snakemake ---------------------------
if(exists("snakemake")){
  pearson_priors_file <- snakemake@input$pearson_prior
  pearson_methods_file <- snakemake@input$pearson_method

  spearman_priors_file <- snakemake@input$spearman_prior
  spearman_methods_file <- snakemake@input$spearman_method

  jaccard_up_priors_file <- snakemake@input$jaccard_up_prior
  jaccard_down_priors_file <- snakemake@input$jaccard_down_prior

  jaccard_up_methods_file <- snakemake@input$jaccard_up_method
  jaccard_down_methods_file <- snakemake@input$jaccard_down_method

  cor_priors_pearson_p <- snakemake@output$cor_plot_pearson_prior
  cor_priors_spearman_p <- snakemake@output$cor_plot_spearman_prior
  cor_priors_jacc_up_p <- snakemake@output$cor_plot_jacc_up_prior
  cor_priors_jacc_down_p <- snakemake@output$cor_plot_jacc_down_prior

  cor_methods_pearson_p <- snakemake@output$cor_plot_pearson_method
  cor_methods_spearman_p <- snakemake@output$cor_plot_spearman_method
  cor_methods_jacc_up_p <- snakemake@output$cor_plot_jacc_up_method
  cor_methods_jacc_down_p <- snakemake@output$cor_plot_jacc_down_method
}else{
  pearson_priors_file <- "results/04_exploration/merged/correlation/correlation_priors_pearson.rds"
  pearson_methods_file <- "results/04_exploration/merged/correlation/correlation_methods_pearson.rds"

  spearman_priors_file <- "results/04_exploration/merged/correlation/correlation_priors_spearman.rds"
  spearman_methods_file <- "results/04_exploration/merged/correlation/correlation_methods_spearman.rds"

  jaccard_up_priors_file <- "results/04_exploration/merged/jaccard/jaccard_priors_up_10.rds"
  jaccard_down_priors_file <- "results/04_exploration/merged/jaccard/jaccard_priors_down_10.rds"

  jaccard_up_methods_file <- "results/04_exploration/merged/jaccard/jaccard_methods_up_10.rds"
  jaccard_down_methods_file <- "results/04_exploration/merged/jaccard/jaccard_methods_down_10.rds"

  cor_priors_pearson_p <- "results/manuscript_figures/figure_3/supp/corrplot_pearson_priors.pdf"
  cor_priors_spearman_p <- "results/manuscript_figures/figure_3/supp/corrplot_spearman_priors.pdf"
  cor_priors_jacc_up_p <- "results/manuscript_figures/figure_3/supp/corrplot_jacc_up_priors.pdf"
  cor_priors_jacc_down_p <- "results/manuscript_figures/figure_3/supp/corrplot_jacc_down_priors.pdf"

  cor_methods_pearson_p <- "results/manuscript_figures/figure_3/supp/corrplot_pearson_methods.pdf"
  cor_methods_spearman_p <- "results/manuscript_figures/figure_3/supp/corrplot_spearman_methods.pdf"
  cor_methods_jacc_up_p <- "results/manuscript_figures/figure_3/supp/corrplot_jacc_up_methods.pdf"
  cor_methods_jacc_down_p <- "results/manuscript_figures/figure_3/supp/corrplot_jacc_down_methods.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(corrplot)

## Method comparison ---------------------------
## Pearson
method_comparison <- readRDS(pearson_methods_file)
tmp <- data.frame(tmp = colnames(method_comparison[[1]])) %>%
  mutate(method = recode(tmp,
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

Y <- do.call(cbind, method_comparison)
Y_array <- array(Y, dim=c(dim(method_comparison[[1]]), length(method_comparison)))

mean_correlation <- colMeans(aperm(Y_array, c(3, 1, 2)), na.rm = TRUE)
colnames(mean_correlation) <- tmp$method
rownames(mean_correlation) <- tmp$method

pdf(cor_methods_pearson_p, height = 3, width = 3.5)
corrplot(mean_correlation,
         col.lim=c(min(mean_correlation), 1),
         is.corr = FALSE, type = 'lower',
         tl.col = 'black', order = 'hclust', cl.pos = 'r', col = COL1('Blues'))
dev.off()

tmp <- mean_correlation %>%
  as.numeric()
tmp[!tmp == 1] %>%
  quantile(probs = 0.8)

tmp2 <- mean_correlation[colnames(mean_correlation) == "KARP"]
tmp2[!tmp2 == 1] %>%
  range()

## Spearman
method_comparison <- readRDS(spearman_methods_file)
tmp <- data.frame(tmp = colnames(method_comparison[[1]])) %>%
  mutate(method = recode(tmp,
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

Y <- do.call(cbind, method_comparison)
Y_array <- array(Y, dim=c(dim(method_comparison[[1]]), length(method_comparison)))

mean_correlation <- colMeans(aperm(Y_array, c(3, 1, 2)), na.rm = TRUE)
colnames(mean_correlation) <- tmp$method
rownames(mean_correlation) <- tmp$method

pdf(cor_methods_spearman_p, height = 3, width = 3.5)
corrplot(mean_correlation,
         col.lim=c(min(mean_correlation), 1),
         is.corr = FALSE, type = 'lower',
         tl.col = 'black', order = 'hclust', cl.pos = 'r', col = COL1('Reds'))
dev.off()

tmp <- mean_correlation %>%
  as.numeric()
tmp[!tmp == 1] %>%
  quantile(probs = 0.8)

tmp2 <- mean_correlation[colnames(mean_correlation) == "KARP"]
tmp2[!tmp2 == 1] %>%
  range()

## Jaccard
method_comparison <- readRDS(jaccard_up_methods_file)
method_comparison <- map(method_comparison, as.matrix)
tmp <- data.frame(tmp = colnames(method_comparison[[1]])) %>%
  mutate(method = recode(tmp,
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

Y <- do.call(cbind, method_comparison)
Y_array <- array(Y, dim=c(dim(method_comparison[[1]]), length(method_comparison)))

mean_correlation <- colMeans(aperm(Y_array, c(3, 1, 2)), na.rm = TRUE)
colnames(mean_correlation) <- tmp$method
rownames(mean_correlation) <- tmp$method

pdf(cor_methods_jacc_up_p, height = 3, width = 3.5)
corrplot(mean_correlation,
         col.lim=c(min(mean_correlation), 1),
         is.corr = FALSE, type = 'lower',
         tl.col = 'black', order = 'hclust', cl.pos = 'r', col = COL1('Purples'))
dev.off()

tmp <- mean_correlation %>%
  as.numeric()
tmp[!tmp == 1] %>%
  quantile(probs = 0.8)

tmp2 <- mean_correlation[colnames(mean_correlation) == "KARP"]
tmp2[!tmp2 == 1] %>%
  range()

method_comparison <- readRDS(jaccard_down_methods_file)
method_comparison <- map(method_comparison, as.matrix)
tmp <- data.frame(tmp = colnames(method_comparison[[1]])) %>%
  mutate(method = recode(tmp,
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

Y <- do.call(cbind, method_comparison)
Y_array <- array(Y, dim=c(dim(method_comparison[[1]]), length(method_comparison)))

mean_correlation <- colMeans(aperm(Y_array, c(3, 1, 2)), na.rm = TRUE)
colnames(mean_correlation) <- tmp$method
rownames(mean_correlation) <- tmp$method

pdf(cor_methods_jacc_down_p, height = 3, width = 3.5)
corrplot(mean_correlation,
         col.lim=c(min(mean_correlation), 1),
         is.corr = FALSE, type = 'lower',
         tl.col = 'black', order = 'hclust', cl.pos = 'r', col = COL1('Purples'))
dev.off()

tmp <- mean_correlation %>%
  as.numeric()
tmp[!tmp == 1] %>%
  quantile(probs = 0.8)

tmp2 <- mean_correlation[colnames(mean_correlation) == "KARP"]
tmp2[!tmp2 == 1] %>%
  range()

## Priors comparison ---------------------------
## Pearson
method_comparison <- readRDS(pearson_priors_file)
method_comparison <- method_comparison[!names(method_comparison) == "number_of_targets"]
tmp <- data.frame(tmp = colnames(method_comparison[[1]])) %>%
  mutate(method = recode(tmp,
                         "iKiPdb" = "iKiP-DB",
                         "GPS" = "GPS gold",
                         "omnipath" = "OmniPath",
                         "networkin" = "NetworKIN",
                         "phosphositeplus" = "PhosphoSitePlus",
                         "ptmsigdb" = "PTMsigDB",
                         "combined" = "extended combined",
                         "GSknown" = "curated combined"))

Y <- do.call(cbind, method_comparison)
Y_array <- array(Y, dim=c(dim(method_comparison[[1]]), length(method_comparison)))

mean_correlation <- colMeans(aperm(Y_array, c(3, 1, 2)), na.rm = TRUE)
colnames(mean_correlation) <- tmp$method
rownames(mean_correlation) <- tmp$method

pdf(cor_priors_pearson_p, height = 3, width = 3.5)
corrplot(mean_correlation,
         col.lim=c(min(mean_correlation), 1),
         is.corr = FALSE, type = 'lower',
         tl.col = 'black', order = 'hclust', cl.pos = 'r', col = COL1('Blues'))
dev.off()

tmp <- mean_correlation %>%
  as.numeric()
tmp[!tmp == 1] %>%
  quantile(probs = 0.8)


## Spearman
method_comparison <- readRDS(spearman_priors_file)
tmp <- data.frame(tmp = colnames(method_comparison[[1]])) %>%
  mutate(method = recode(tmp,
                         "iKiPdb" = "iKiP-DB",
                         "GPS" = "GPS gold",
                         "omnipath" = "OmniPath",
                         "networkin" = "NetworKIN",
                         "phosphositeplus" = "PhosphoSitePlus",
                         "ptmsigdb" = "PTMsigDB",
                         "combined" = "extended combined",
                         "GSknown" = "curated combined"))

Y <- do.call(cbind, method_comparison)
Y_array <- array(Y, dim=c(dim(method_comparison[[1]]), length(method_comparison)))

mean_correlation <- colMeans(aperm(Y_array, c(3, 1, 2)), na.rm = TRUE)
colnames(mean_correlation) <- tmp$method
rownames(mean_correlation) <- tmp$method

pdf(cor_priors_spearman_p, height = 3, width = 3.5)
corrplot(mean_correlation,
         col.lim=c(min(mean_correlation), 1),
         is.corr = FALSE, type = 'lower',
         tl.col = 'black', order = 'hclust', cl.pos = 'r', col = COL1('Reds'))
dev.off()

tmp <- mean_correlation %>%
  as.numeric()
tmp[!tmp == 1] %>%
  quantile(probs = 0.8)


## Jaccard
method_comparison <- readRDS(jaccard_up_priors_file)
method_comparison <- map(method_comparison, as.matrix)
tmp <- data.frame(tmp = colnames(method_comparison[[1]])) %>%
  mutate(method = recode(tmp,
                         "iKiPdb" = "iKiP-DB",
                         "GPS" = "GPS gold",
                         "omnipath" = "OmniPath",
                         "networkin" = "NetworKIN",
                         "phosphositeplus" = "PhosphoSitePlus",
                         "ptmsigdb" = "PTMsigDB",
                         "combined" = "extended combined",
                         "GSknown" = "curated combined"))

Y <- do.call(cbind, method_comparison)
Y_array <- array(Y, dim=c(dim(method_comparison[[1]]), length(method_comparison)))

mean_correlation <- colMeans(aperm(Y_array, c(3, 1, 2)), na.rm = TRUE)
colnames(mean_correlation) <- tmp$method
rownames(mean_correlation) <- tmp$method

pdf(cor_priors_jacc_up_p, height = 3, width = 3.5)
corrplot(mean_correlation,
         col.lim=c(min(mean_correlation), 1),
         is.corr = FALSE, type = 'lower',
         tl.col = 'black', order = 'hclust', cl.pos = 'r', col = COL1('Purples'))
dev.off()

tmp <- mean_correlation %>%
  as.numeric()
tmp[!tmp == 1] %>%
  quantile(probs = 0.8)

tmp2 <- mean_correlation[colnames(mean_correlation) == "KARP"]
tmp2[!tmp2 == 1] %>%
  range()

method_comparison <- readRDS(jaccard_down_priors_file)
method_comparison <- map(method_comparison, as.matrix)
tmp <- data.frame(tmp = colnames(method_comparison[[1]])) %>%
  mutate(method = recode(tmp,
                         "iKiPdb" = "iKiP-DB",
                         "GPS" = "GPS gold",
                         "omnipath" = "OmniPath",
                         "networkin" = "NetworKIN",
                         "phosphositeplus" = "PhosphoSitePlus",
                         "ptmsigdb" = "PTMsigDB",
                         "combined" = "extended combined",
                         "GSknown" = "curated combined"))

Y <- do.call(cbind, method_comparison)
Y_array <- array(Y, dim=c(dim(method_comparison[[1]]), length(method_comparison)))

mean_correlation <- colMeans(aperm(Y_array, c(3, 1, 2)), na.rm = TRUE)
colnames(mean_correlation) <- tmp$method
rownames(mean_correlation) <- tmp$method

pdf(cor_priors_jacc_down_p, height = 3, width = 3.5)
corrplot(mean_correlation,
         col.lim=c(min(mean_correlation), 1),
         is.corr = FALSE, type = 'lower',
         tl.col = 'black', order = 'hclust', cl.pos = 'r', col = COL1('Purples'))
dev.off()

tmp <- mean_correlation %>%
  as.numeric()
tmp[!tmp == 1] %>%
  quantile(probs = 0.8)

