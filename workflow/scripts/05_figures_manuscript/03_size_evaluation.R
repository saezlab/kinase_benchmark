if(exists("snakemake")){
  rank_files <- snakemake@input$rank
  activating_strat <- snakemake@input$act
  tumor_strat <- snakemake@input$tumor
  performance_pert <- snakemake@output$plt_pert
  performance_act <- snakemake@output$plt_act
  performance_tumor <- snakemake@output$plt_tumor
}else{
  rank_files <- list.files("results/03_benchmark/merged/02_mean_rank",
                           pattern = "csv", recursive = TRUE, full.names = T)
  rank_files <- rank_files[str_detect(rank_files, "/iKiPdb/|/GPS/|/omnipath/|/networkin/|/phosphositeplus/|/ptmsigdb/|/GSknown/")]
  activating_strat <- list.files("data/results_cptac/stratified_kinases/all_kins/actsiteBM", full.names = T)
  tumor_strat <- list.files("data/results_cptac/stratified_kinases/all_kins/protBM", full.names = T)
  performance_pert <- "results/manuscript_figures/figure_3/regulon_perturbation.pdf"
  performance_act <- "results/manuscript_figures/figure_3/regulon_activatingsites.pdf"
  performance_tumor <- "results/manuscript_figures/figure_3/regulon_tumor.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(ggpubr)
library(patchwork)

## Benchmark ------------------
## Load AUROC ---------------------------
rank_list <- map(rank_files, function(file){
  read_csv(file, col_types = cols())
})

rank_df <- bind_rows(rank_list)

# Stratify kinases
kinase_strat <- rank_df %>%
  filter(prior == "GSknown") %>%
  group_by(targets) %>%
  summarise(measured_targets = mean(measured_targets)) %>%
  dplyr::mutate(strat = case_when(
   measured_targets < 11 ~ "small",
   measured_targets > 10 & measured_targets < 26 ~ "medium",
   measured_targets > 25 ~ "large"
  )) %>%
  select(targets, strat) %>%
  dplyr::filter(!is.na(strat))

tmp <- rank_df %>%
  filter(prior == "GSknown") %>%
  group_by(sample, targets) %>%
  summarise(measured_targets = mean(measured_targets)) %>%
  dplyr::mutate(strat = case_when(
   measured_targets < 11 ~ "small",
   measured_targets > 10 & measured_targets < 26 ~ "medium",
   measured_targets > 25 ~ "large"
  )) %>%
  select(sample, targets, strat) %>%
  dplyr::filter(!is.na(strat))

table(tmp$strat)

rank_df <- rank_df %>%
  dplyr::select(method, prior, rank, scaled_rank, sample, targets) %>%
  filter(!is.na(rank)) %>%
  mutate(prior = recode(prior,
                        "iKiPdb" = "iKiP-DB",
                        "GPS" = "GPS gold",
                        "omnipath" = "OmniPath",
                        "networkin" = "NetworKIN",
                        "phosphositeplus" = "PhosphoSitePlus",
                        "ptmsigdb" = "PTMsigDB",
                        "combined" = "extended combined",
                        "GSknown" = "curated",
                        "shuffled2" = "shuffled")) %>%
  left_join(kinase_strat, by = "targets") %>%
  dplyr::filter(!is.na(strat))

order_net <- c("curated", "GPS gold", "PhosphoSitePlus", "PTMsigDB", "NetworKIN", "OmniPath", "iKiP-DB", "shuffled")
rank_df$prior <- factor(rank_df$prior, levels = order_net)
rank_df$strat <- factor(rank_df$strat, levels = c("small", "medium", "large"))

n_kinases <- rank_df %>%
  group_by(prior, method, strat) %>%
  summarise(kinases = length(unique(targets))) %>%
  filter(method == "zscore") %>%
  select(prior, kinases, strat) %>%
  add_column(benchmark = "perturbation-based")
n_kinases$prior <- factor(n_kinases$prior, levels = order_net)
n_kinases$strat <- factor(n_kinases$strat, levels = c("small", "medium", "large"))


auroc_p <- ggplot(rank_df, aes(x = prior, y = scaled_rank, fill = strat)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.1) +
  scale_fill_manual(values = c("#96c3df", "#4292C6", "#265c7f")) +  # Muted scientific color palette
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0, "lines"),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 11),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("") +
  ylab("scaled rank")


kin_p <- ggplot(n_kinases, aes(x = prior, y = kinases, fill = strat)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.4)+ # Line connecting the dots
  scale_y_continuous(
    name = "tmp",
    limits = c(0, 15),
    breaks = seq(0, 15, by = 5) #
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
  scale_fill_manual(values = c("#96c3df", "#4292C6", "#265c7f")) +
  ggtitle("Number of Kinases in Evaluation Set")

full_p <- ggarrange(kin_p, auroc_p, ncol = 1, common.legend = T, heights = c(3.1, 9))

pdf(performance_pert, width = 3.9, height = 3.5)
full_p
dev.off()


## tumor benchmark
df_act <-  map_dfr(activating_strat, function(act_file){
  res <- readRDS(act_file)
  roc_list <- lapply(res$ROC_results, "[[", "sample_AUROCs")
  net_id <- str_extract(act_file, "(?<=BM_).*?(?=_)")
  strat_id <- str_extract(act_file, paste0("(?<=", net_id, "_).*?(?=_roc)"))

  data.frame(method = "z-score", prior = net_id, auroc = roc_list$zscore, strat=strat_id )
}) %>%
  mutate(prior = recode(prior,
                      "ikip" = "iKiP-DB",
                      "gps" = "GPS gold",
                      "omni" = "OmniPath",
                      "nwkin" = "NetworKIN",
                      "psp" = "PhosphoSitePlus",
                      "ptmsig" = "PTMsigDB",
                      "combo" = "extended combined",
                      "known" = "curated",
                      "shuffled2" = "shuffled")) %>%
  mutate(strat = recode(strat,
                      "0to10" = "small",
                      "11to25" = "medium",
                      "26plus" = "large")) 
                      
df_act <- df_act %>%
  filter(prior %in% rank_df$prior)

df_act$prior <- factor(df_act$prior, levels = order_net)
df_act$strat <- factor(df_act$strat, levels = c("small", "medium", "large"))

n_kinases <-  map_dfr(activating_strat, function(act_file){
  res <- readRDS(act_file)
  kin <- res$evaluation_kinases %>%
    unlist() %>%
    unique() %>%
    length()

  net_id <- str_extract(act_file, "(?<=BM_).*?(?=_)")
  strat_id <- str_extract(act_file, paste0("(?<=", net_id, "_).*?(?=_roc)"))

  data.frame(method = "z-score", prior = net_id, kinases = kin, strat=strat_id )
}) %>%
  mutate(prior = recode(prior,
                      "ikip" = "iKiP-DB",
                      "gps" = "GPS gold",
                      "omni" = "OmniPath",
                      "nwkin" = "NetworKIN",
                      "psp" = "PhosphoSitePlus",
                      "ptmsig" = "PTMsigDB",
                      "combo" = "extended combined",
                      "known" = "curated",
                      "shuffled2" = "shuffled")) %>%
  mutate(strat = recode(strat,
                      "0to10" = "small",
                      "11to25" = "medium",
                      "26plus" = "large")) 

n_kinases <- n_kinases %>%
  filter(prior %in% df_act$prior)

n_kinases$prior <- factor(n_kinases$prior, levels = order_net)
n_kinases$strat <- factor(n_kinases$strat, levels = c("small", "medium", "large"))

auroc_p <- ggplot(df_act, aes(x = prior, y = auroc, fill = strat)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.1) +
  scale_fill_manual(values = c("#e2b99e", "#C67642", "#58331a")) +  # Muted scientific color palette
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0, "lines"),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 11),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("") +
  ylab("AUROC") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5)

kin_p <- ggplot(n_kinases, aes(x = prior, y = kinases, fill = strat)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.4)+ # Line connecting the dots
  scale_y_continuous(
    name = "tmp",
    limits = c(0, 30),
    breaks = seq(0, 30, by = 10) #
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
  scale_fill_manual(values = c("#e2b99e", "#C67642", "#58331a")) +
  ggtitle("Number of Kinases in Evaluation Set")

full_p <- ggarrange(kin_p, auroc_p, ncol = 1, common.legend = T, heights = c(3.1, 9))

pdf(performance_act, width = 3.9, height = 3.5)
full_p
dev.off()


## Tumor
df_prot <-  map_dfr(tumor_strat, function(act_file){
  res <- readRDS(act_file)
  roc_list <- lapply(res$ROC_results, "[[", "sample_AUROCs")
  net_id <- str_extract(act_file, "(?<=BM_).*?(?=_)")
  strat_id <- str_extract(act_file, paste0("(?<=", net_id, "_).*?(?=_roc)"))

  data.frame(method = "z-score", prior = net_id, auroc = roc_list$zscore, strat=strat_id )
}) %>%
  mutate(prior = recode(prior,
                      "ikip" = "iKiP-DB",
                      "gps" = "GPS gold",
                      "omni" = "OmniPath",
                      "nwkin" = "NetworKIN",
                      "psp" = "PhosphoSitePlus",
                      "ptmsig" = "PTMsigDB",
                      "combo" = "extended combined",
                      "known" = "curated",
                      "shuffled2" = "shuffled")) %>%
  mutate(strat = recode(strat,
                      "0to10" = "small",
                      "11to25" = "medium",
                      "26plus" = "large")) 
                      
df_prot <- df_prot %>%
  filter(prior %in% rank_df$prior)

df_prot$prior <- factor(df_prot$prior, levels = order_net)
df_prot$strat <- factor(df_prot$strat, levels = c("small", "medium", "large"))

n_kinases_prot <-  map_dfr(tumor_strat, function(act_file){
  res <- readRDS(act_file)
  kin <- res$evaluation_kinases %>%
    unlist() %>%
    unique() %>%
    length()

  net_id <- str_extract(act_file, "(?<=BM_).*?(?=_)")
  strat_id <- str_extract(act_file, paste0("(?<=", net_id, "_).*?(?=_roc)"))

  data.frame(method = "z-score", prior = net_id, kinases = kin, strat=strat_id )
}) %>%
  mutate(prior = recode(prior,
                      "ikip" = "iKiP-DB",
                      "gps" = "GPS gold",
                      "omni" = "OmniPath",
                      "nwkin" = "NetworKIN",
                      "psp" = "PhosphoSitePlus",
                      "ptmsig" = "PTMsigDB",
                      "combo" = "extended combined",
                      "known" = "curated",
                      "shuffled2" = "shuffled")) %>%
  mutate(strat = recode(strat,
                      "0to10" = "small",
                      "11to25" = "medium",
                      "26plus" = "large")) 

n_kinases_prot <- n_kinases_prot %>%
  filter(prior %in% df_act$prior)

n_kinases_prot$prior <- factor(n_kinases_prot$prior, levels = order_net)
n_kinases_prot$strat <- factor(n_kinases_prot$strat, levels = c("small", "medium", "large"))

auroc_p <- ggplot(df_prot, aes(x = prior, y = auroc, fill = strat)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.1) +
  scale_fill_manual(values = c("#d9a4a2", "#b54d4a", "#572523")) +  # Muted scientific color palette
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing.x = unit(0, "lines"),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 11),
    text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
    axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
    axis.title.x = element_text(family = "Helvetica", size = 10)
  ) +
  xlab("") +
  ylab("AUROC") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5)

kin_p <- ggplot(n_kinases_prot, aes(x = prior, y = kinases, fill = strat)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.4)+ # Line connecting the dots
  scale_y_continuous(
    name = "tmp",
    limits = c(0, 50),
    breaks = seq(0, 50, by = 20) #
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
  scale_fill_manual(values = c("#d9a4a2", "#b54d4a", "#572523")) +
  ggtitle("Number of Kinases in Evaluation Set")

full_p <- ggarrange(kin_p, auroc_p, ncol = 1, common.legend = T, heights = c(3.1, 9))

pdf(performance_tumor, width = 3.9, height = 3.5)
full_p
dev.off()
