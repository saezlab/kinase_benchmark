#'

## Snakemake ---------------------------
if(exists("snakemake")){
  bench_files <- snakemake@input$bench_files
  bench_files_2 <- snakemake@input$predicted
  output_file <- snakemake@output$merged
  output_file_2 <- snakemake@output$pred

}else{
  bench_files <- list.files("results/03_benchmark/merged/02_benchmark_res_subset/predicted",
                            pattern = "bench", recursive = TRUE, full.names = T)
  output_file <- "results/manuscript_figures/supp_figures/combined.pdf"
  bench_files_full <- list.files("results/03_benchmark/merged/02_benchmark_res",
                            pattern = "bench", recursive = TRUE, full.names = T)
  bench_files_full <- bench_files_full[str_detect(bench_files_full, "phosphositeplus")]
  bench_files_full <- bench_files_full[!str_detect(bench_files_full, "johnson|phosformer")]
  output_file_full <- "results/manuscript_figures/supp_figures/combined_full.pdf"
  bench_files_2 <- list.files("results/03_benchmark/merged/02_benchmark_res_subset/johnson",
                            pattern = "bench", recursive = TRUE, full.names = T)
  output_file_2 <- "results/manuscript_figures/supp_figures/predicted_combined.pdf"
  bench_files_3 <- list.files("results/03_benchmark/merged/02_benchmark_res",
                              pattern = "bench", recursive = TRUE, full.names = T)
  bench_files_3 <- bench_files_3[str_detect(bench_files_3, "phosphositeplus|johnson|phosformer")]
}

## Libraries ---------------------------
library(tidyverse)


## AUROC results ---------------------------
bench_list <- map(bench_files, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][6])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]
bench_df <- bind_rows(bench_list)

order_m <- bench_df %>%
  filter(metric == "mcauroc") %>%
  group_by(net) %>%
  summarise(m_score = mean(score)) %>%
  arrange(desc(m_score)) %>%
  pull(net)
order_m <- order_m[!order_m == "phosphositeplus"]
bench_df$net <- factor(bench_df$net, levels = c("phosphositeplus", order_m))

order_method <- bench_df %>%
  filter(metric == "mcauroc") %>%
  group_by(method) %>%
  summarise(m_score = mean(score)) %>%
  arrange(desc(m_score)) %>%
  pull(method)
bench_df$method <- factor(bench_df$method, levels = c("RoKAI_z", "KSEA_z", "ptmsea", "mean", "fgsea", "number_of_targets"))

custom_palette <- c("#66c2a5", "#fc8d62", "#8da0cb")

auroc_p <- bench_df %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=method, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.8) +
  theme_bw() +
  scale_fill_manual(values = custom_palette) +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.4, 'cm')) +
  geom_hline(yintercept = 0.5, linetype="dotted", lwd = 1) +
  ylab("AUROC") +
  xlab("")

pdf(output_file, width = 6, height = 3)
auroc_p
dev.off()

## full results
bench_list <- map(bench_files_full, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][5])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]
bench_df <- bind_rows(bench_list) %>%
  filter(method %in% c("KSEA_z", "RoKAI_z", "fgsea", "mean", "number_of_targets", "ptmsea"))

order_m <- bench_df %>%
  filter(metric == "mcauroc") %>%
  group_by(net) %>%
  summarise(m_score = mean(score)) %>%
  arrange(desc(m_score)) %>%
  pull(net)
order_m <- order_m[!order_m == "phosphositeplus"]
bench_df$net <- factor(bench_df$net, levels = c("phosphositeplus", rev(order_m)))

order_method <- bench_df %>%
  filter(metric == "mcauroc") %>%
  group_by(method) %>%
  summarise(m_score = mean(score)) %>%
  arrange(desc(m_score)) %>%
  pull(method)
bench_df$method <- factor(bench_df$method, levels = c("RoKAI_z", "KSEA_z", "ptmsea", "mean", "fgsea", "number_of_targets"))

custom_palette <- c("#66c2a5", "#fc8d62", "#8da0cb")

auroc_p <- bench_df %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=method, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.8) +
  theme_bw() +
  scale_fill_manual(values = custom_palette) +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.4, 'cm')) +
  geom_hline(yintercept = 0.5, linetype="dotted", lwd = 1) +
  ylab("AUROC") +
  xlab("")

pdf(output_file_full, width = 6, height = 3)
auroc_p
dev.off()

bench_list <- map(bench_files_2, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][6])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]
bench_df <- bind_rows(bench_list)

order_m <- bench_df %>%
  filter(metric == "mcauroc") %>%
  group_by(net) %>%
  summarise(m_score = mean(score)) %>%
  arrange(desc(m_score)) %>%
  pull(net)
order_m <- order_m[!order_m == "phosphositeplus"]
bench_df$net <- factor(bench_df$net, levels = c("phosphositeplus", rev(order_m)))

order_method <- bench_df %>%
  filter(metric == "mcauroc") %>%
  group_by(method) %>%
  summarise(m_score = mean(score)) %>%
  arrange(desc(m_score)) %>%
  pull(method)
bench_df$method <- factor(bench_df$method, levels = order_method)

custom_palette <- c("#66c2a5", "#fc8d62", "#8da0cb")

auroc_p <- bench_df %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=method, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.8) +
  theme_bw() +
  scale_fill_manual(values = custom_palette) +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.4, 'cm')) +
  geom_hline(yintercept = 0.5, linetype="dotted", lwd = 1) +
  ylab("AUROC") +
  xlab("")

pdf(output_file_2, width = 6, height = 3)
auroc_p
dev.off()
