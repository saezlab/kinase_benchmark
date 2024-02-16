#' In this script we perform the statistical analysis of our results

## Snakemake ---------------------------
if(exists("snakemake")){
  bench_files <- snakemake@input$bench_files
  meta_file <- snakemake@input$meta
  meta_hijazi <- snakemake@input$hijazi
  auroc_plot <- snakemake@output$auroc
  overview_meta <- snakemake@output$meta_over
  rank_file <- snakemake@input$rank
  rank_plot <- snakemake@output$rankPlt
  kinase_rank <- snakemake@input$kin
  kin_file <- snakemake@output$rankKin
  prior_size <- snakemake@input$prior
}else{
  bench_files <- list.files("results/hijazi/06_benchmark_res",
                            pattern = "bench", recursive = TRUE, full.names = T)
  bench_files <- bench_files[str_detect(bench_files, "merged")]
  bench_files <- bench_files[!str_detect(bench_files, "number_of_targets")]
  bench_files <- bench_files[!str_detect(bench_files, "GSknown")]
}

## Libraries ---------------------------
library(tidyverse)
library(corrplot)
library(ComplexHeatmap)

## Comparison of benchmark results---------------------------
bench_list <- map(bench_files, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][4])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]

bench_df <- bind_rows(bench_list) %>%
  filter(metric == "mcauroc") %>%
  filter(method != "number_of_targets")

order_net <- bench_df %>%
  group_by(net) %>%
  summarise(AUROC = median(score)) %>%
  arrange(desc(AUROC)) %>%
  pull(net)

order_method <- bench_df %>%
  group_by(method) %>%
  summarise(AUROC = median(score)) %>%
  arrange(desc(AUROC)) %>%
  pull(method)

bench_df$net <- factor(bench_df$net, levels = order_net)
bench_df$method <- factor(bench_df$method, levels = order_method)

### change format into matrix for AUROC and AUPRC
prior_mat <- bench_df %>%
  group_by(net, method) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  add_column(counter = rep(c(1:(length(unique(bench_df$method)))), times = length(unique(bench_df$net)))) %>%
  select(score, net, counter) %>%
  pivot_wider(names_from = net, values_from = score) %>%
  column_to_rownames("counter")

method_mat <- bench_df%>%
  group_by(net, method) %>%
  summarise(score = mean(score)) %>%
  ungroup() %>%
  add_column(counter = rep(c(1:(length(unique(bench_df$net)))), each = length(unique(bench_df$method)))) %>%
  select(score, method, counter) %>%
  pivot_wider(names_from = method, values_from = score) %>%
  column_to_rownames("counter")

comb_mat <-  bench_df %>%
  mutate(comp = paste(net, method, sep = ":")) %>%
  add_column(counter = rep(c(1:1000), times = (length(unique(bench_df$net))*length(unique(bench_df$method))))) %>%
  dplyr::select(score, comp, counter) %>%
  pivot_wider(names_from = comp, values_from = score) %>%
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
  mutate(t.value = round(t.value, digits = 1)) %>%
  mutate(p.adj = round(p.adj, digits = 2))

prior.wilcox <- prior.wilcox %>%
  mutate(comp1_new = case_when(
    t.value > 0 ~ comp1,
    t.value == 0 ~ comp1,
    t.value < 0 ~ comp2))  %>%
  mutate(comp2_new = case_when(
    t.value > 0 ~ comp2,
    t.value == 0 ~ comp2,
    t.value < 0 ~ comp1))  %>%
  mutate(t.value_new = abs(t.value)) %>%
  select(comp1_new, comp2_new, p.adj, t.value_new, p.value) %>%
  rename("comp1" = comp1_new, "comp2" = comp2_new, "t.value" = t.value_new)

prior.wilcox <- prior.wilcox %>% arrange(comp1, comp2)

# method
method.wilcox <- perform.multi.wilcox(method_mat) %>%
  mutate(t.value = round(t.value, digits = 1)) %>%
  mutate(p.adj = round(p.adj, digits = 2))

method.wilcox <- method.wilcox %>%
  mutate(comp1_new = case_when(
    t.value > 0 ~ comp1,
    t.value == 0 ~ comp1,
    t.value < 0 ~ comp2))  %>%
  mutate(comp2_new = case_when(
    t.value > 0 ~ comp2,
    t.value == 0 ~ comp2,
    t.value < 0 ~ comp1))  %>%
  mutate(t.value_new = abs(t.value)) %>%
  select(comp1_new, comp2_new, p.adj, t.value_new, p.value) %>%
  rename("comp1" = comp1_new, "comp2" = comp2_new, "t.value" = t.value_new)

method.wilcox <- method.wilcox %>% arrange(comp1, comp2)

## Comb
comb.wilcox <- perform.multi.wilcox(comb_mat) %>%
  mutate(t.value = round(t.value, digits = 1)) %>%
  mutate(p.adj = round(p.adj, digits = 2))

comb.wilcox <- comb.wilcox %>%
  mutate(comp1_new = case_when(
    t.value > 0 ~ comp1,
    t.value == 0 ~ comp1,
    t.value < 0 ~ comp2))  %>%
  mutate(comp2_new = case_when(
    t.value > 0 ~ comp2,
    t.value == 0 ~ comp2,
    t.value < 0 ~ comp1))  %>%
  mutate(t.value_new = abs(t.value)) %>%
  select(comp1_new, comp2_new, p.adj, t.value_new, p.value) %>%
  rename("comp1" = comp1_new, "comp2" = comp2_new, "t.value" = t.value_new)

comb.wilcox <- comb.wilcox %>% arrange(comp1, comp2)


