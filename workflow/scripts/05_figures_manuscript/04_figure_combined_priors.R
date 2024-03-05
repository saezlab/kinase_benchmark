#'

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
  heat_plot <- snakemake@output$heat
  medRank_plot <- snakemake@output$medRank
  gene_citations <- snakemake@input$cit
  cor_plot <- snakemake@output$cor
}else{
  bench_files <- list.files("results/hijazi/06_benchmark_res",
                            pattern = "bench", recursive = TRUE, full.names = T)
  bench_files <- bench_files[str_detect(bench_files, "subset")]
  prior_files <- list.files("results/hernandez/prior", full.names = T)
  prior_files <- prior_files[str_detect(prior_files, "phosphositeplus")]
  output_plot <- "results/manuscript_figures/figure_4/predicted_auroc.pdf"
}

## Libraries ---------------------------
library(tidyverse)


## AUROC results ---------------------------
bench_list <- map(bench_files, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][4])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]
bench_df <- bind_rows(bench_list)

order_m <- c("phosphositeplus", "phosphositeplus_networkin", "phosphositeplus_iKiPdb")
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
  geom_hline(yintercept = 0.5, linetype="dotted", lwd = 1) +
  ylab("AUROC") +
  xlab("")

pdf(output_plot, width = 8, height = 5)
auroc_p + theme(legend.position = "bottom")
dev.off()

prior <- map(prior_files, read_tsv)
names(prior) <- str_remove(map_chr(str_split(prior_files, "/"), 4), ".tsv")
overview_df <- map_dfr(names(prior), function(prior_idx){
  prior_df <- prior[[prior_idx]]
  n_kin <- length(unique(prior_df$source))
  n_edges <- nrow(prior_df)

  sum(table(prior_df$source) >= 5)

  data.frame(prior = prior_idx,
             kinases = n_kin,
             kinase5 = sum(table(prior_df$source) >= 5),
             edges = n_edges)
})

ggplot(overview_df, aes(x = prior, y = edges)) +
  geom_bar(stat = "identity")
