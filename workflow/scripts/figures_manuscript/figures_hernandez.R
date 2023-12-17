#'

## Snakemake ---------------------------
if(exists("snakemake")){
  input_file <- snakemake@input
  output_file <- snakemake@output
}else{
  input_file <- "path_to_input"
  output_file <- "path_to_output"
}

## Libraries ---------------------------
library(tidyverse)

## Load data ---------------------------
hernandez <- read_csv("results/hernandez/processed_data/benchmark_data.csv", col_types = cols()) %>%
  column_to_rownames("ID")

hernandez_meta <- read_csv("results/hernandez/processed_data/benchmark_metadata.csv", col_types = cols())

priors <- list.files("results/hernandez/prior", pattern = ".tsv", full.names = T)
all_targets <- map(priors, function(p_file){
  read_tsv(p_file, col_types = cols()) %>% pull(target)
}) %>%
  unlist %>%
  unique() %>%
  str_remove("\\|auto")

bench_samples <- list.files("results/hernandez/benchmark_files", pattern = "obs_", full.names = T)
bench_samples <- bench_samples[!str_detect(bench_samples, "number_of_targets")]
all_samples <- map(bench_samples, function(b_file){
  read_csv(b_file, col_types = cols()) %>% pull(sample)
}) %>%
  unlist %>%
  unique()

hernandez_filtered <- hernandez[colnames(hernandez) %in% all_samples]
hernandez_meta <- hernandez_meta %>%
  filter(id %in% colnames(hernandez_filtered))

## Figure benchmark ---------------------------
### Coverage ---------------------------
coverage_df <- data.frame(experiment = colnames(hernandez_filtered),
                          measures_pps = colSums(!is.na(hernandez_filtered)))

measured_in_prior <- map_dbl(colnames(hernandez_filtered), function(col_i){
  msk <- !is.na(hernandez_filtered[[col_i]])
  sum(rownames(hernandez_filtered)[msk] %in% all_targets)
})

coverage_df <- coverage_df %>%
  add_column(pps_in_prior = measured_in_prior) %>%
  mutate(not_in_prior = measures_pps - pps_in_prior)

coverage_df <- coverage_df %>%
  pivot_longer(!experiment, names_to = "pps_covered", values_to = "n_pps") %>%
  mutate(pps_covered = recode(pps_covered,
                              "measures_pps" = "total",
                              "pps_in_prior" = "known",
                              "not_in_prior" = "unknown"))

mean_measured_covered <- mean(coverage_df %>% filter(pps_covered == "known") %>% pull(n_pps))
sd_measured_covered <- sd(coverage_df %>% filter(pps_covered == "known") %>% pull(n_pps))

coverage_df$pps_covered <- factor(coverage_df$pps_covered, levels = c("unknown", "known"))

coverage_df %>% filter(is.na(pps_covered)) %>% pull(n_pps) %>% mean()
coverage_df %>%
  mutate(pps_covered = case_when(
    is.na(pps_covered) ~ "all",
    !is.na(pps_covered) ~ pps_covered
  )) %>%
  filter(!pps_covered == "unknown") %>%
  pivot_wider(names_from = "pps_covered", values_from = "n_pps") %>%
  mutate(percentage = known/all) %>%
  pull(percentage) %>%
  mean

overview_p <- ggplot(coverage_df %>%
         filter(!pps_covered == "total"), aes(fill=pps_covered, y=n_pps, x=experiment)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("#D3D3D3", "#6a6a6a")) +
  geom_hline(yintercept = mean_measured_covered) +
  #ggtitle("Phosphoproteome datasets", subtitle = paste(round(mean_measured_covered, digits = 0),
  #                                                      label="±",
  #                                                      round(sd_measured_covered, digits = 1),
  #                                                      sep = " ")) +
  xlab(paste0(length(unique(coverage_df$experiment)), " experiments")) +
  ylab("# of unique phosphorylation sites") +
  guides(fill=guide_legend(title="upstream kinase"))+
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 10),
        legend.position = "bottom",
        text = element_text(size = 10))

pdf("results/manuscript_figures/figure_2/overview_experiment.pdf", width = 2.2, height = 3.0)
overview_p
dev.off()

### kinases ---------------------------
kin_df <- table(hernandez_meta$target, hernandez_meta$sign) %>%
  as.data.frame() %>%
  rename("kinase" = Var1) %>%
  rename("perturbation" = Var2) %>%
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
        legend.title = element_text(size = 9),
        legend.position = "bottom",
        text = element_text(size = 9))

pdf("results/manuscript_figures/hernandez_benchmark/02_overview_kin.pdf", width = 2.7, height = 3.15)
kin_p
dev.off()
