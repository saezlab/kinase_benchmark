#'

## Snakemake ---------------------------
if(exists("snakemake")){
  bench_files <- snakemake@input$bench_files
  meta_file <- snakemake@input$meta
  auroc_plot <- snakemake@output$auroc
  overview_meta <- snakemake@output$meta_over
  rank_files <- snakemake@input$rank
  rank_plot <- snakemake@output$rankPlt
  kin_file <- snakemake@output$rankKin
  heat_plot <- snakemake@output$heat
  medRank_plot <- snakemake@output$medRank
  gene_citations <- snakemake@input$cit
  cor_plot <- snakemake@output$cor
  ppsp_file <- snakemake@input$ppsp
  cor_target_plot <- snakemake@output$target
  prior_csv <- snakemake@output$statPrior
  method_csv <- snakemake@output$statMethod
}else{
  bench_files <- list.files("results/03_benchmark/merged/02_benchmark_res",
                            pattern = "bench", recursive = TRUE, full.names = T)
  meta_file <- "results/01_processed_data/merged/data/benchmark_metadata.csv"
  rank_files <- list.files("results/03_benchmark/merged/02_mean_rank/",
                            recursive = TRUE, full.names = T)
  gene_citations <- "resources/protein_citations.csv"
  overview_meta <- "results/manuscript_figures/figure_3/overview_kin.pdf"
  auroc_plot <- "results/manuscript_figures/figure_3/auroc_res.pdf"
  rank_plot <- "results/manuscript_figures/figure_3/mean_rank.pdf"
  kin_file <- "results/manuscript_figures/figure_3/kinase_GSknown.csv"
  heat_plot <- "results/manuscript_figures/figure_3/median_auroc.pdf"
  medRank_plot <- "results/manuscript_figures/figure_3/median_rank.pdf"
  cor_plot <- "results/manuscript_figures/figure_3/study_bias.pdf"
  ppsp_file <- "results/01_processed_data/cptac/mapped_priors/phosphositeplus.tsv"
  cor_target_plot <- "results/manuscript_figures/figure_3/target_bias.pdf"
  prior_csv <- "results/manuscript_figures/supp_files/prior_comparison.csv"
  method_csv <- "results/manuscript_figures/supp_files/method_comparison.csv"
}

## Libraries ---------------------------
library(tidyverse)
library(corrplot)
library(ComplexHeatmap)
library(circlize)

## Overview kinases ---------------------------
hernandez_meta <- read_csv(meta_file, col_types = cols()) %>%
    dplyr::select(id, sign, target)
hernandez_meta <- separate_rows(hernandez_meta, target, sep = ";")

kin_df <- table(hernandez_meta$target, hernandez_meta$sign) %>%
  as.data.frame() %>%
  dplyr::rename("kinase" = Var1) %>%
  dplyr::rename("perturbation" = Var2) %>%
  mutate(perturbation = recode(perturbation,
                               "1" = "up",
                               "-1" = "down"))
kin_order <- kin_df %>%
  group_by(kinase) %>%
  summarise(total = sum(Freq)) %>%
  arrange(desc(total)) %>%
  pull(kinase)
kin_df$kinase <- factor(kin_df$kinase, levels = rev(kin_order))

kin_p <- ggplot(kin_df, aes(x = kinase, y = Freq, fill = perturbation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  xlab("") +
  ylab("# perturbations") +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 11),
        legend.position = "bottom",
        text = element_text(size = 10))

length(unique(kin_df$kinase))
kin_df %>% group_by(perturbation) %>% summarise(total = sum(Freq))

pdf(overview_meta, width = 2.7, height = 7)
kin_p
dev.off()


## AUROC results ---------------------------
bench_list <- map(bench_files, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][5])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]
bench_df <- bind_rows(bench_list)

order_m <- bench_df %>%
  filter(metric == "mcauroc") %>%
  group_by(net) %>%
  summarise(m_score = mean(score)) %>%
  arrange(desc(m_score)) %>%
  pull(net)
bench_df$net <- factor(bench_df$net, levels = order_m)

order_method <- bench_df %>%
  filter(metric == "mcauroc") %>%
  group_by(method) %>%
  summarise(m_score = mean(score)) %>%
  arrange(desc(m_score)) %>%
  pull(method)
bench_df$method <- factor(bench_df$method, levels = order_method)

custom_palette <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494")

auroc_p <- bench_df %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=method, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.8) +
  theme_bw() +
  scale_fill_manual(values = custom_palette) +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.4, 'cm')) +
  ylim(0.5, 0.83) +
  geom_hline(yintercept = 0.5, linetype="dotted", lwd = 1) +
  ylab("AUROC") +
  xlab("")

pdf(auroc_plot, width = 6.5, height = 3)
auroc_p
dev.off()

## AUROC Median heatmap
med_mat <- bench_df %>%
  group_by(net, method) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  pivot_wider(names_from = "method", values_from = "score") %>%
  column_to_rownames("net")
med_mat

col_fun = colorRamp2(
  c(min(med_mat, na.rm = TRUE), mean(unlist(as.vector(med_mat)), na.rm = TRUE), max(med_mat, na.rm = TRUE)),
  c("white", "#edcac6", "#7F0863")
)

#col_fun = colorRamp2(c(min(med_mat, na.rm = T), max(med_mat, na.rm = T)), c("white", "deeppink4"))
column_ha = HeatmapAnnotation(mean_method = colMeans(med_mat, na.rm = T), col = list(mean_method = col_fun))
row_ha = rowAnnotation(mean_prior = rowMeans(med_mat, na.rm = T), col = list(mean_prior = col_fun))

pdf(heat_plot, width = 8, height = 4)
Heatmap(med_mat, row_split = 4, column_split = 5,border = TRUE, rect_gp = gpar(col = "white", lwd = 1), top_annotation = column_ha, left_annotation = row_ha, col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", med_mat[i, j]), x, y, gp = gpar(fontsize = 6))})
dev.off()

## Mean rank ---------------------------
priors <- map_chr(str_split(bench_files, "/"), 5) %>% unique()
rank_df <- map_dfr(rank_files, function(x) read_csv(x, col_types = cols()))
rank_df <- rank_df %>%
  filter(prior %in% priors)
rank_df$prior <- factor(rank_df$prior, levels = order_m)
rank_df$method <- factor(rank_df$method, levels = order_method)

rank_p <- rank_df %>%
  filter(!is.na(scaled_rank)) %>%
  ggplot(aes(x=method, y=scaled_rank, fill=prior)) +
  geom_boxplot(outlier.size=0.2, lwd=0.8) +
  theme_bw() +
  scale_fill_manual(values = custom_palette) +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.4, 'cm')) +
  geom_hline(yintercept = 0.5, linetype="dotted", lwd = 1) +
  ylab("scaled rank") +
  xlab("")

pdf(rank_plot, width = 6.5, height = 3)
rank_p
dev.off()

## scaled Median rank heatmap
medRank_mat <- rank_df %>%
  filter(!is.na(scaled_rank)) %>%
  group_by(prior, method) %>%
  summarise(scaled_rank = mean(scaled_rank)) %>%
  ungroup() %>%
  pivot_wider(names_from = "method", values_from = "scaled_rank") %>%
  column_to_rownames("prior")
medRank_mat

col_fun = colorRamp2(
  c(min(medRank_mat, na.rm = TRUE), mean(unlist(as.vector(medRank_mat)), na.rm = TRUE), max(medRank_mat, na.rm = TRUE)),
  c("#4770b2", "#d6e2f2", "white")
)

#col_fun = colorRamp2(c(min(medRank_mat), max(medRank_mat)), c("deepskyblue4", "white"))
column_ha = HeatmapAnnotation(mean_method = colMeans(medRank_mat, na.rm = T), col = list(mean_method = col_fun))
row_ha = rowAnnotation(mean_prior = rowMeans(medRank_mat, na.rm = T), col = list(mean_prior = col_fun))

pdf(medRank_plot, width = 8, height = 4)
Heatmap(medRank_mat, row_split = 4, column_split = 5,border = TRUE, rect_gp = gpar(col = "white", lwd = 1), top_annotation = column_ha, left_annotation = row_ha, col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", medRank_mat[i, j]), x, y, gp = gpar(fontsize = 6))})
dev.off()

## Mean rank per kinase ---------------------------
kin_gsknown <- rank_df %>%
  filter(prior == "phosphositeplus") %>%
  filter(!method == "number_of_targets") %>%
  group_by(targets) %>%
  summarise(mean_rank = mean(rank, na.rm = T)) %>%
  filter(!is.na(mean_rank)) %>%
  arrange(mean_rank) %>%
  mutate(mean_rank = round(mean_rank)) %>%
  dplyr::select(targets, mean_rank)

write_csv(kin_gsknown, kin_file)

## Compare to number of citations
overview_citations <- read_csv(gene_citations, col_types = cols())

kin_citations <- kin_gsknown %>%
  left_join(overview_citations %>% rename("targets" = Symbol), by = "targets")

cor_p <- ggplot(kin_citations, aes(x = mean_rank, y = log(citations))) +
  geom_point(size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.6, fullrange = TRUE, color = "steelblue3") +
  theme_bw() +
  theme(text = element_text(size = 11),
        legend.key.size = unit(0.4, 'cm'))  +
  geom_text(aes(label = targets), check_overlap = TRUE, size = 2.2, nudge_x = 0.1, nudge_y = 0.1)+
  ylab("log (# pubmed IDs)") +
  xlab("mean rank") +
  ggtitle(paste0("Pearson correlation: ", round(cor(kin_citations$mean_rank, log(kin_citations$citations)), digits = 2)))

cor.test(kin_citations$mean_rank, log(kin_citations$citations))

pdf(cor_plot, width = 3, height = 2.7)
cor_p
dev.off()

## Compare to number of targets
overview_targets <- read_tsv(ppsp_file, col_types = cols())
overview_targets <- overview_targets %>%
  group_by(source) %>%
  summarise(n = n())

kin_targets <- kin_gsknown %>%
  left_join(overview_targets %>% rename("targets" = source), by = "targets")

cor_targets_p <- ggplot(kin_targets, aes(x = mean_rank, y = log(n))) +
  geom_point(size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.6, fullrange = TRUE, color = "steelblue3") +
  theme_bw() +
  theme(text = element_text(size = 11),
        legend.key.size = unit(0.4, 'cm'))  +
  geom_text(aes(label = targets), check_overlap = TRUE, size = 2.2, nudge_x = 0.1, nudge_y = 0.1)+
  ylab("log (# targets)") +
  xlab("mean rank") +
  ggtitle(paste0("Pearson correlation: ", round(cor(kin_targets$mean_rank, log(kin_targets$n)), digits = 2)))

cor.test(kin_targets$mean_rank, log(kin_targets$n))

pdf(cor_target_plot, width = 3, height = 2.7)
cor_targets_p
dev.off()

## Statistics -----------
### change format into matrix for AUROC and AUPRC
prior_mat <- bench_df %>%
  filter(method != "number_of_targets") %>%
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

### Perform t-test
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

# AUROC
prior.wilcox <- perform.multi.wilcox(prior_mat) %>%
  mutate(wilcoxon = round(t.value, digits = 1)) %>%
  select(-t.value)

prior.wilcox$p.adj[prior.wilcox$p.adj == 0] <- 2.2e-16

write_csv(prior.wilcox, prior_csv)

# method
method.wilcox <- perform.multi.wilcox(method_mat) %>%
  mutate(wilcoxon = round(t.value, digits = 1)) %>%
  select(-t.value)

method.wilcox$p.adj[method.wilcox$p.adj == 0] <- 2.2e-16

write_csv(method.wilcox, method_csv)


