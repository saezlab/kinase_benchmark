#'

## Snakemake ---------------------------
if(exists("snakemake")){
  bench_files <- snakemake@input$bench_files
  activating_files <-  snakemake@input$act
  tumor_files <-  snakemake@input$tumor

  heat_plot <-  snakemake@output$plt
  heat_plot_tumor <-  snakemake@output$plt_tumor
  heat_plot_act <-  snakemake@output$plt_act

  prior_csv <- snakemake@output$csv_prior
  method_csv <- snakemake@output$csv_method
  prior_csv_tumor <- snakemake@output$csv_prior_tumor
  method_csv_tumor <- snakemake@output$csv_method_tumor
  prior_csv_act <- snakemake@output$csv_prior_act
  method_csv_act <- snakemake@output$csv_method_act
}else{
  bench_files <- list.files("results/03_benchmark/merged/02_benchmark_res",
                            pattern = "bench", recursive = TRUE, full.names = T)
  bench_files <- bench_files[str_detect(bench_files, "/shuffled2/|/iKiPdb/|/GPS/|/omnipath/|/networkin/|/phosphositeplus/|/ptmsigdb/|/GSknown/")]
  activating_files <- list.files("data/results_cptac/overall_performance/actsiteBM/all_kins",
                                 pattern = "table", full.names = T)
  tumor_files <- list.files("data/results_cptac/overall_performance/protBM/all_kins",
                                 pattern = "table", full.names = T)

  heat_plot <- "results/manuscript_figures/figure_3/supp/median_auroc.pdf"
  heat_plot_tumor <- "results/manuscript_figures/figure_3/supp/median_auroc_tumor.pdf"
  heat_plot_act <- "results/manuscript_figures/figure_3/supp/median_auroc_act.pdf"
  
  prior_csv <- "results/manuscript_figures/supp_files/prior_comparison_perturbation.csv"
  method_csv <- "results/manuscript_figures/supp_files/method_comparison_perturbation.csv"
  prior_csv_tumor <- "results/manuscript_figures/supp_files/prior_comparison_tumor.csv"
  method_csv_tumor <- "results/manuscript_figures/supp_files/method_comparison_tumor.csv"
  prior_csv_act <- "results/manuscript_figures/supp_files/prior_comparison_act.csv"
  method_csv_act <- "results/manuscript_figures/supp_files/method_comparison_act.csv"
}

## Libraries ---------------------------
library(tidyverse)
library(corrplot)
library(ComplexHeatmap)
library(circlize)

## AUROC results ---------------------------
bench_list <- map(bench_files, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][5])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]
bench_df <- bind_rows(bench_list) %>%
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
                         "zscore" = "z-score"))  %>%
  mutate(net = recode(net,
                      "iKiPdb" = "iKiP-DB",
                      "GPS" = "GPS gold",
                      "omnipath" = "OmniPath",
                      "networkin" = "NetworKIN",
                      "phosphositeplus" = "PhosphoSitePlus",
                      "ptmsigdb" = "PTMsigDB",
                      "combined" = "extended combined",
                      "GSknown" = "curated",
                      "shuffled2" = "shuffled")) %>%
  add_column(benchmark = "perturbation-based") %>%
  dplyr::select(net, score, method, benchmark)

med_mat <- bench_df %>%
  group_by(net, method) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  pivot_wider(names_from = "method", values_from = "score") %>%
  column_to_rownames("net")
med_mat

col_fun = colorRamp2(
  c(min(med_mat, na.rm = TRUE), mean(unlist(as.vector(med_mat)), na.rm = TRUE), max(med_mat, na.rm = TRUE)),
  c("white", "#d6e2f2", "#4770b2")
)

#col_fun = colorRamp2(c(min(med_mat, na.rm = T), max(med_mat, na.rm = T)), c("white", "deeppink4"))
column_ha = HeatmapAnnotation(mean_method = colMeans(med_mat, na.rm = T), col = list(mean_method = col_fun))
row_ha = rowAnnotation(mean_prior = rowMeans(med_mat, na.rm = T), col = list(mean_prior = col_fun))

pdf(heat_plot, width = 8, height = 4)
Heatmap(med_mat, row_split = 4, column_split = 5,border = TRUE, rect_gp = gpar(col = "white", lwd = 1), top_annotation = column_ha, left_annotation = row_ha, col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", med_mat[i, j]), x, y, gp = gpar(fontsize = 6))})
dev.off()

write_csv(data.frame(med_mat) %>% rownames_to_column("prior"), "results/manuscript_data/suppfig8a.csv")

## tumor bench ---------------------------
df_tumor <- map_dfr(tumor_files, function(file_roc){
    roc <- readRDS(file_roc)
    net_id <- str_extract(file_roc, "(?<=5perThr_).*?(?=_roc_)")
    map_dfr(colnames(roc), function(method_id){
        data.frame(net = net_id, score = roc[,colnames(roc) == method_id], method = method_id, benchmark = "tumor-based")
    })
    
}) %>%
  mutate(net = recode(net,
                      "ikip" = "iKiP-DB",
                      "gps" = "GPS gold",
                      "omni" = "OmniPath",
                      "nwkin" = "NetworKIN",
                      "psp" = "PhosphoSitePlus",
                      "ptmsig" = "PTMsigDB",
                      "combo" = "extended combined",
                      "known" = "curated",
                      "shuffled2" = "shuffled")) %>%
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
  dplyr::select(net, score, method, benchmark)


med_mat <- df_tumor %>%
  group_by(net, method) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  pivot_wider(names_from = "method", values_from = "score") %>%
  column_to_rownames("net")
med_mat

col_fun = colorRamp2(
  c(min(med_mat, na.rm = TRUE), mean(unlist(as.vector(med_mat)), na.rm = TRUE), max(med_mat, na.rm = TRUE)),
  c("white", "#e6c3c2", "#b54d4a")
)

#col_fun = colorRamp2(c(min(med_mat, na.rm = T), max(med_mat, na.rm = T)), c("white", "deeppink4"))
column_ha = HeatmapAnnotation(mean_method = colMeans(med_mat, na.rm = T), col = list(mean_method = col_fun))
row_ha = rowAnnotation(mean_prior = rowMeans(med_mat, na.rm = T), col = list(mean_prior = col_fun))

pdf(heat_plot_tumor, width = 8, height = 4)
Heatmap(med_mat, row_split = 4, column_split = 5,border = TRUE, rect_gp = gpar(col = "white", lwd = 1), top_annotation = column_ha, left_annotation = row_ha, col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", med_mat[i, j]), x, y, gp = gpar(fontsize = 6))})
dev.off()

write_csv(data.frame(med_mat) %>% rownames_to_column("prior"), "results/manuscript_data/suppfig8c.csv")
## activating site
# activating sites benchmark
df_act <- map_dfr(activating_files, function(file_roc){
    roc <- readRDS(file_roc)
    net_id <- str_extract(file_roc, "(?<=5perThr_).*?(?=_roc_)")
    map_dfr(colnames(roc), function(method_id){
        data.frame(net = net_id, score = roc[,colnames(roc) == method_id], method = method_id, benchmark = "activating site")
    })   
}) %>%
  mutate(net = recode(net,
                      "ikip" = "iKiP-DB",
                      "gps" = "GPS gold",
                      "omni" = "OmniPath",
                      "nwkin" = "NetworKIN",
                      "psp" = "PhosphoSitePlus",
                      "ptmsig" = "PTMsigDB",
                      "combo" = "extended combined",
                      "known" = "curated",
                      "shuffled2" = "shuffled")) %>%
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
  dplyr::select(net, score, method, benchmark)

med_mat <- df_act %>%
  group_by(net, method) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  pivot_wider(names_from = "method", values_from = "score") %>%
  column_to_rownames("net")
med_mat

col_fun = colorRamp2(
  c(min(med_mat, na.rm = TRUE), mean(unlist(as.vector(med_mat)), na.rm = TRUE), max(med_mat, na.rm = TRUE)),
  c("white", "#e3bda4", "#C67642")
)

#col_fun = colorRamp2(c(min(med_mat, na.rm = T), max(med_mat, na.rm = T)), c("white", "deeppink4"))
column_ha = HeatmapAnnotation(mean_method = colMeans(med_mat, na.rm = T), col = list(mean_method = col_fun))
row_ha = rowAnnotation(mean_prior = rowMeans(med_mat, na.rm = T), col = list(mean_prior = col_fun))

pdf(heat_plot_act, width = 8, height = 4)
Heatmap(med_mat, row_split = 4, column_split = 5,border = TRUE, rect_gp = gpar(col = "white", lwd = 1), top_annotation = column_ha, left_annotation = row_ha, col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", med_mat[i, j]), x, y, gp = gpar(fontsize = 6))})
dev.off()

write_csv(data.frame(med_mat) %>% rownames_to_column("prior"), "results/manuscript_data/suppfig8b.csv")
## Exploration -----------
comb_df <- rbind(bench_df, df_tumor, df_act) %>%
  group_by(method, net, benchmark) %>%
  summarise(score = mean(score)) %>%
  mutate(id = paste(method, net, sep = "_")) %>%
  ungroup() 

comb_df %>% arrange(desc(score)) %>% group_by(id) %>% summarise(score = mean(score)) %>% arrange(desc(score))
comb_df %>% filter(!net == "shuffled") %>% group_by(method, benchmark) %>% summarise(score = mean(score)) %>% ungroup() %>% arrange(desc(score)) %>% group_by(benchmark) %>% group_split()
comb_df %>% filter(!method == "n_targets") %>% group_by(net, benchmark) %>% summarise(score = mean(score)) %>% ungroup() %>% arrange(desc(score)) %>% group_by(benchmark) %>% group_split()
comb_df %>% filter(net == "curated" & method == "z-score") %>% summarise(score = mean(score)) 

mat_ppsp <- comb_df %>% filter(net == "PhosphoSitePlus") %>% select(method, benchmark, score) %>% pivot_wider(names_from = "benchmark", values_from = "score") %>% column_to_rownames("method")
cor(mat_ppsp)

cor.test(mat_ppsp$`activating site`, mat_ppsp$`tumor-based`, method = "pearson")
cor.test(mat_ppsp$`activating site`, mat_ppsp$`perturbation-based`, method = "pearson")
cor.test(mat_ppsp$`perturbation-based`, mat_ppsp$`tumor-based`, method = "pearson")

mat_comb <- comb_df %>% filter(!net == "shuffled") %>% filter(!method == "n_targets") %>% select(id, benchmark, score) %>% pivot_wider(names_from = "benchmark", values_from = "score") %>% column_to_rownames("id")

cor.test(mat_comb$`activating site`, mat_comb$`tumor-based`, method = "pearson")
cor.test(mat_comb$`activating site`, mat_comb$`perturbation-based`, method = "pearson")
cor.test(mat_comb$`perturbation-based`, mat_comb$`tumor-based`, method = "pearson")
## Correlation ----------- 
comb_mat <- comb_df %>%
  dplyr::select(id, score, benchmark) %>%
  pivot_wider(names_from = "benchmark", values_from = "score") %>%
  column_to_rownames("id")

comb_mat %>% 
  cor(method = "spearman")

comb_mat %>%
  cor(method = "pearson")


## Statistics -----------
multi.wilcox <- function(mat, pVal = T, alternative = "two.sided") {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- wilcox.test(mat[, i], mat[, j], alternative = "two.sided")
      if(pVal){
        p.mat[i, j] <- p.mat[j, i] <- test$p.value
      } else {
        p.mat[i, j] <- p.mat[j, i] <- test$statistic
      }

    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

perform.multi.wilcox <- function(mat){

  # extract p_values
  wilcox.p <- multi.wilcox(mat) %>%
    as.data.frame() %>%
    rownames_to_column("comp1") %>%
    pivot_longer(!comp1, names_to = "comp2", values_to = "p.value")

  # extract t_values
  wilcox.t <- multi.wilcox(mat, pVal = F) %>%
    as.data.frame() %>%
    rownames_to_column("comp1") %>%
    pivot_longer(!comp1, names_to = "comp2", values_to = "t.value")

  # Adjust t-values direction
  # Identify positions were direction needs to be adjusted
  idx <- c()
  for (i in 1:(length(unique(wilcox.p$comp1))-1)){
    new_idx <- length(unique(wilcox.t$comp1))*i + rep(1:i)
    idx <- append(idx, new_idx)
  }
  wilcox.t$t.value[idx] <- wilcox.t$t.value[idx]*-1

  wilcox <- wilcox.p %>%
    add_column(p.adj = p.adjust(wilcox.p$p.value, method = "BH")) %>%
    left_join(wilcox.t, by = c("comp1", "comp2")) %>%
    filter(!comp1 == comp2)

  # remove duplicated comparisons
  idx_rm <- c()
  for (i in 1:(length(unique(wilcox.t$comp1))-1)){
    new_idx <- (length(unique(wilcox.t$comp1))-1)*i + rep(1:i)
    idx_rm <- append(idx_rm, new_idx)
  }

  wilcox <- wilcox %>%
    slice(-idx_rm)

  return(wilcox)
}

## Perturbation
### change format into matrix for AUROC and AUPRC
prior_mat <- bench_df %>%
  filter(method != "n targets") %>%
  group_by(net, method) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  add_column(counter = rep(c(1:(length(unique(bench_df$method))-1)), times = length(unique(bench_df$net)))) %>%
  select(score, net, counter) %>%
  pivot_wider(names_from = net, values_from = score) %>%
  column_to_rownames("counter")

method_mat <- bench_df %>%
  filter(net != "shuffled") %>%
  group_by(net, method) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  add_column(counter = rep(c(1:(length(unique(bench_df$net))-1)), each = length(unique(bench_df$method)))) %>%
  select(score, method, counter) %>%
  pivot_wider(names_from = method, values_from = score) %>%
  column_to_rownames("counter")

### Perform t-test
# AUROC
prior.wilcox <- perform.multi.wilcox(prior_mat) %>%
  mutate(wilcoxon = round(t.value, digits = 1)) %>%
  select(-t.value) %>%
  arrange(desc(wilcoxon))

prior.wilcox$p.adj[prior.wilcox$p.adj == 0] <- 2.2e-16

write_csv(prior.wilcox, prior_csv)

prior.wilcox %>%
  filter(p.adj <= 0.05) %>%
  arrange(desc(wilcoxon))

# method
method.wilcox <- perform.multi.wilcox(method_mat) %>%
  mutate(wilcoxon = round(t.value, digits = 1)) %>%
  select(-t.value)  %>%
  arrange(desc(wilcoxon))

method.wilcox$p.adj[method.wilcox$p.adj == 0] <- 2.2e-16

write_csv(method.wilcox, method_csv)

method.wilcox %>%
  filter(p.adj <= 0.05)  %>%
  arrange(desc(wilcoxon))

## Tumor based
### change format into matrix for AUROC and AUPRC
prior_mat <- df_tumor %>%
  filter(method != "n targets") %>%
  group_by(net, method) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  add_column(counter = rep(c(1:(length(unique(df_tumor$method))-1)), times = length(unique(df_tumor$net)))) %>%
  select(score, net, counter) %>%
  pivot_wider(names_from = net, values_from = score) %>%
  column_to_rownames("counter")

method_mat <- df_tumor %>%
  filter(net != "shuffled") %>%
  group_by(net, method) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  add_column(counter = rep(c(1:(length(unique(df_tumor$net))-1)), each = length(unique(df_tumor$method)))) %>%
  select(score, method, counter) %>%
  pivot_wider(names_from = method, values_from = score) %>%
  column_to_rownames("counter")

### Perform t-test
# AUROC
prior.wilcox <- perform.multi.wilcox(prior_mat) %>%
  mutate(wilcoxon = round(t.value, digits = 1)) %>%
  select(-t.value) %>%
  arrange(desc(wilcoxon))

prior.wilcox$p.adj[prior.wilcox$p.adj == 0] <- 2.2e-16

write_csv(prior.wilcox, prior_csv_tumor)

prior.wilcox %>%
  filter(p.adj <= 0.05) %>%
  arrange(desc(wilcoxon))

# method
method.wilcox <- perform.multi.wilcox(method_mat) %>%
  mutate(wilcoxon = round(t.value, digits = 1)) %>%
  select(-t.value) %>%
  arrange(desc(wilcoxon))

method.wilcox$p.adj[method.wilcox$p.adj == 0] <- 2.2e-16

write_csv(method.wilcox, method_csv_tumor)

method.wilcox %>%
  filter(p.adj <= 0.05) %>%
  arrange(desc(wilcoxon))

## Activating sites based
### change format into matrix for AUROC and AUPRC
prior_mat <- df_act %>%
  filter(method != "n targets") %>%
  group_by(net, method) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  add_column(counter = rep(c(1:(length(unique(df_act$method))-1)), times = length(unique(df_act$net)))) %>%
  select(score, net, counter) %>%
  pivot_wider(names_from = net, values_from = score) %>%
  column_to_rownames("counter")

method_mat <- df_act %>%
  filter(net != "shuffled") %>%
  group_by(net, method) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  add_column(counter = rep(c(1:(length(unique(df_act$net))-1)), each = length(unique(df_act$method)))) %>%
  select(score, method, counter) %>%
  pivot_wider(names_from = method, values_from = score) %>%
  column_to_rownames("counter")

### Perform t-test
# AUROC
prior.wilcox <- perform.multi.wilcox(prior_mat) %>%
  mutate(wilcoxon = round(t.value, digits = 1)) %>%
  select(-t.value) %>%
  arrange(desc(wilcoxon))

prior.wilcox$p.adj[prior.wilcox$p.adj == 0] <- 2.2e-16

write_csv(prior.wilcox, prior_csv_act)

prior.wilcox %>%
  filter(p.adj <= 0.05) %>%
  arrange(desc(wilcoxon))

# method
method.wilcox <- perform.multi.wilcox(method_mat) %>%
  mutate(wilcoxon = round(t.value, digits = 1)) %>%
  select(-t.value) %>%
  arrange(desc(wilcoxon))

method.wilcox$p.adj[method.wilcox$p.adj == 0] <- 2.2e-16

write_csv(method.wilcox, method_csv_act)

method.wilcox %>%
  filter(p.adj <= 0.05) %>%
  arrange(desc(wilcoxon))
