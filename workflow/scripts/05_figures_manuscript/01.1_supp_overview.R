if(exists("snakemake")){
  prior_files <- snakemake@input$prior_files
  kin_pdf <- snakemake@output$kin
  kinase_pdf <- snakemake@output$kin_heat
  edg_pdf <- snakemake@output$edges
  edge_pdf <- snakemake@output$edges_heat
}else{
  meta_file <- "results/01_processed_data/merged/data/benchmark_metadata.csv"
  overview_meta <- "results/manuscript_figures/figure_1/supp/overview_kin.pdf"
    
  bench_files <- list.files("results/03_benchmark/hijazi/02_benchmark_res/phosphositeplus",
                            pattern = "bench", recursive = TRUE, full.names = T)
  bench_files <- c(bench_files, list.files("results/03_benchmark/hijaziDiscoverX/02_benchmark_res/phosphositeplus",
                            pattern = "bench", recursive = TRUE, full.names = T))
  rank_files <- list.files("results/03_benchmark/hijazi/02_mean_rank/phosphositeplus",
                           pattern = "csv", recursive = TRUE, full.names = T)
  rank_files <- c(rank_files, list.files("results/03_benchmark/hijaziDiscoverX/02_mean_rank/phosphositeplus",
                           pattern = "csv", recursive = TRUE, full.names = T))
  out_plot <- "results/manuscript_figures/figure_1/supp/hijazi.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(ggpubr)

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

length(unique(hernandez_meta$id))
length(unique(kin_df$kinase))
kin_df %>% group_by(perturbation) %>% summarise(total = sum(Freq))

pdf(overview_meta, width = 2.7, height = 7.2)
kin_p
dev.off()

## Benchmark ------------------
## Load AUROC ---------------------------
rank_list <- map(rank_files, function(file){
  read_csv(file, col_types = cols()) %>%
    add_column(bench = str_split(file, "\\/")[[1]][3])
})

rank_df <- bind_rows(rank_list) %>%
  dplyr::select(method, prior, rank, scaled_rank, sample, targets, bench) %>%
  filter(!is.na(rank)) %>%
  mutate(prior = recode(prior, "phosphositeplus" = "PhosphoSitePlus"))

n_kinases <- rank_df %>%
  group_by(bench, method) %>%
  summarise(kinases = length(unique(targets))) %>%
  select(method, bench, kinases)

bench_list <- map(bench_files, function(file){
  read.csv(file, col.names = c("rows", "groupby", "group", "source",
                               "method", "metric", "score", "ci")) %>%
    dplyr::select(-rows) %>%
    add_column(net = str_split(file, "/")[[1]][5]) %>%
    add_column(bench = str_split(file, "\\/")[[1]][3])
})

bench_list <- bench_list[!(map_dbl(bench_list, nrow) == 0)]
df_perturb_all <- bind_rows(bench_list)  %>%
  mutate(net = recode(net,"phosphositeplus" = "PhosphoSitePlus"))

df_perturb_all %>% group_by(method) %>% summarise(auroc = mean(score)) %>% arrange(desc(auroc))
df_perturb_all %>% group_by(bench) %>% summarise(auroc = mean(score)) %>% arrange(desc(auroc))

df_perturb <- df_perturb_all %>%
  dplyr::filter(metric == "mcauroc") %>%
  dplyr::select(method, bench, score) 

mean_auroc <- df_perturb %>%
  group_by(method, bench) %>%
  summarize(tmp_auroc = mean(score)) %>%
  ungroup() %>%
  group_by(method) %>%
  summarize(mean_auroc = mean(tmp_auroc)) %>%
  arrange(desc(mean_auroc))  # Sorting methods by mean AUROC

# Update df to ensure methods are ordered based on mean AUROC
df_perturb <- df_perturb %>%
  mutate(method = factor(method, levels = mean_auroc$method))

# Create the boxplot with ggplot2
lines <- mean_auroc$mean_auroc
names(lines) <- mean_auroc$method

auroc_p <- ggplot(df_perturb, aes(x = method, y = score, fill = bench)) +
  geom_boxplot(linewidth = 0.5, outlier.size = 0.8) +
  scale_fill_manual(values = c("#4292C6", "#b54d4a")) +  # Muted scientific color palette
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0, "lines"),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  # Add horizontal lines or dots for mean AUROC
  geom_point(data = data.frame(method = names(lines), mean_auroc = lines), aes(x = method, y = mean_auroc),
             color = "black", size = 2, shape = 3, fill = "white")  +
  xlab("") +
  ylab("AUROC")  +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5)


kin_df <- n_kinases
kin_df$method <- factor(kin_df$method, levels = mean_auroc$method)

kin_df %>% group_by(bench) %>% summarise(targets = mean(kinases))

kin_p <- ggplot(kin_df, aes(x = method, y = kinases, fill = bench)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.4)+ # Line connecting the dots
  scale_y_continuous(
    name = "tmp",
    limits = c(0, 50),
    breaks = seq(0, 50, by = 10) #
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(0, "lines"),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_blank(),   # Remove x-axis text labels
    axis.ticks.x = element_blank()
  )+
  scale_fill_manual(values = c("#4292C6", "#b54d4a")) +
  ggtitle("Number of Kinases in Evaluation Set")

full_p <- ggarrange(kin_p, auroc_p, ncol = 1, common.legend = T, heights = c(3, 9))

pdf(out_plot, width = 7, height = 4)
full_p
dev.off()
