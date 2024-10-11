#'

## Snakemake ---------------------------
if(exists("snakemake")){
  bench_files <- snakemake@input$bench
  rank_files <- snakemake@input$rank
  k_phit <- snakemake@params$k_phit
  performance_plot <- snakemake@output$plot
}else{
  bench_files <- list.files("results/03_benchmark/merged/02_benchmark_res_subset/subset",
                           pattern = "bench", recursive = TRUE, full.names = T)
  rank_files <-  list.files("results/03_benchmark/merged/02_mean_rank_subset/subset", 
                            pattern = "csv", recursive = TRUE, full.names = T)                        
  performance_plot <- "results/04_exploration/merged/benchmark/plots/performance_subset.pdf"
  k_phit <- 10
}

## Libraries ---------------------------
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

## Load AUROC ---------------------------
if (any(str_detect(bench_files, "subset"))){
  net_id <- 6
} else {
  net_id <- 5
}

bench_list <- map(bench_files, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][net_id])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]
bench_df <- bind_rows(bench_list)

med_auroc <- bench_df %>%
  group_by(net, method) %>%
  summarise(auroc = mean(score)) %>%
  ungroup() %>%
  rename("prior" = net)
  
med_auroc_wide <- med_auroc %>%
  pivot_wider(names_from = "method", values_from = "auroc") %>%
  column_to_rownames("prior")

## Load scaled rank and pHit ---------------------------
rank_list <- map(rank_files, function(file){
  read_csv(file, col_types = cols()) %>%
    dplyr::select(method, prior, rank, scaled_rank)
})

rank_df <- bind_rows(rank_list) %>% 
  filter(!is.na(rank))

med_rank <- rank_df %>%
  group_by(prior, method) %>%
  summarise(scaled_rank = mean(scaled_rank)) %>%
  ungroup() 

med_rank_wide <- med_rank %>%
  pivot_wider(names_from = "method", values_from = "scaled_rank") %>%
  column_to_rownames("prior")

pHit <- rank_df %>%
  group_by(prior, method) %>%
  summarise(phit = (sum(rank <= k_phit)/n())) %>%
  ungroup() 

pHit_wide <- pHit %>%
  pivot_wider(names_from = "method", values_from = "phit") %>%
  column_to_rownames("prior")

## Compare metrices ---------------------------
performance <- left_join(med_auroc, med_rank, by = c("prior", "method")) %>%
  left_join(pHit, by = c("prior", "method"))

cor(performance[, 3:ncol(performance)], method = "pearson")
cor.test(performance$scaled_rank, performance$phit, method = "pearson")
cor.test(performance$scaled_rank, performance$auroc, method = "pearson")
cor.test(performance$phit, performance$auroc, method = "pearson")


## Plot performance ---------------------------
col_auroc <- colorRamp2(
  c(min(med_auroc_wide, na.rm = TRUE), mean(unlist(as.vector(med_auroc_wide)), na.rm = TRUE), max(med_auroc_wide, na.rm = TRUE)),
  c("white", "#edcac6", "#7F0863")
)
auroc_p <- Heatmap(med_auroc_wide,
        name = "AUROC",                             # Name for the legend
        col = col_auroc,
        cluster_rows = TRUE,                        # Enable clustering for rows
        cluster_columns = TRUE,                     # Enable clustering for columns
        show_row_names = TRUE,                      # Show row names (Methods)
        show_column_names = TRUE,                   # Show column names (Priors)
        row_title = "Methods",                      # Row title
        column_title = "Priors"
)
col_rank <- colorRamp2(
  c(min(med_rank_wide, na.rm = TRUE), mean(unlist(as.vector(med_rank_wide)), na.rm = TRUE), max(med_rank_wide, na.rm = TRUE)),
  c("#4770b2", "#d6e2f2", "white")
)
rank_p <- Heatmap(med_rank_wide,
        name = "scaled_rank",                             # Name for the legend
        col = col_rank,
        cluster_rows = TRUE,                        # Enable clustering for rows
        cluster_columns = TRUE,                     # Enable clustering for columns
        show_row_names = TRUE,                      # Show row names (Methods)
        show_column_names = TRUE,                   # Show column names (Priors)
        row_title = "Methods",                      # Row title
        column_title = "Priors"
)

col_phit <- colorRamp2(
  c(min(pHit_wide, na.rm = TRUE), mean(unlist(as.vector(pHit_wide)), na.rm = TRUE), max(pHit_wide, na.rm = TRUE)),
  c("white", "#edcac6", "#7F0863")
)
pHit_p <- Heatmap(pHit_wide,
        name = paste0("pHit(k=", k_phit,")"),       # Name for the legend
        col =col_phit,
        cluster_rows = TRUE,                        # Enable clustering for rows
        cluster_columns = TRUE,                     # Enable clustering for columns
        show_row_names = TRUE,                      # Show row names (Methods)
        show_column_names = TRUE,                   # Show column names (Priors)
        row_title = "Methods",                      # Row title
        column_title = "Priors"
)

pdf(performance_plot, width = 8, height = 4)
auroc_p
rank_p
pHit_p
dev.off()
