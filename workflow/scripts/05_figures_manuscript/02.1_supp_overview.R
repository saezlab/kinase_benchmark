if(exists("snakemake")){
  meta_file <- snakemake@input$meta
  bench_files <- snakemake@input$bench
  rank_files <- snakemake@input$rank
  overview_meta <- snakemake@output$ove
  out_plot <- snakemake@output$out
}else{
  meta_file <- "data/datasets/cptac/GSsets/protein_5percent.Rds"
  meta_act_file <- "data/datasets/cptac/GSsets/actsite_5percent.Rds"
  overview_meta <- "results/manuscript_figures/figure_2/supp/overview_kin.pdf"
  overview_pairs <- "results/manuscript_figures/figure_2/supp/overview_pairs.pdf"
  kin_class_file <- "resources/kinase_class.csv"
  kin_overview_prot <- "results/manuscript_figures/figure_2/supp/top_kin_prot.pdf"
  kin_overview_act <- "results/manuscript_figures/figure_2/supp/top_kin_act.pdf"
}


## Libraries ---------------------------
library(tidyverse)

## Load data ---------------------------
kinase_class <- read_csv(kin_class_file, col_types = cols())
meta <- readRDS(meta_file)
tumor <- names(meta)
gs_pos <- map_dfr(tumor, function(tumor_id){
    kinases <- meta[[tumor_id]]$GS_pos_pairs
    map_dfr(names(kinases), function(kin_id){
        data.frame(kinase = kin_id, sample = kinases[[kin_id]], tumor = tumor_id, GS = "positive")
    })
})

gs_neg <- map_dfr(tumor, function(tumor_id){
    kinases <- meta[[tumor_id]]$GS_neg_pairs
    map_dfr(names(kinases), function(kin_id){
        data.frame(kinase = kin_id, sample = kinases[[kin_id]], tumor = tumor_id, GS = "negative")
    })
})

full_gs <- rbind(gs_pos, gs_neg) %>%
    filter(kinase %in% kinase_class$source)

meta_act <- readRDS(meta_act_file)
act <- names(meta_act)
gs_act_pos <- map_dfr(act, function(tumor_id){
    kinases <- meta_act[[tumor_id]]$GS_pos_pairs
    map_dfr(names(kinases), function(kin_id){
        data.frame(kinase = kin_id, sample = kinases[[kin_id]], tumor = tumor_id, GS = "positive")
    })
})

gs_act_neg <- map_dfr(act, function(tumor_id){
    kinases <- meta_act[[tumor_id]]$GS_neg_pairs
    map_dfr(names(kinases), function(kin_id){
        data.frame(kinase = kin_id, sample = kinases[[kin_id]], tumor = tumor_id, GS = "negative")
    })
})

full_act_gs <- rbind(gs_act_pos, gs_act_neg) %>%
    filter(kinase %in% kinase_class$source)


## Overview kinases ---------------------------
kin_overview <- full_gs %>%
  group_by(GS) %>%
  summarise(n_pairs = n(), n_kin = length(unique(kinase))) %>%
  ungroup() %>%
  add_column(benchmark = "protein-based")

kin_act_overview <- full_act_gs %>%
  group_by(GS) %>%
  summarise(n_pairs = n(), n_kin = length(unique(kinase))) %>%
  ungroup() %>%
  add_column(benchmark = "activating site-based")

kin_full <- rbind(kin_overview, kin_act_overview)

kin_p <- ggplot(kin_full, aes(x = GS, y = n_kin, fill = benchmark)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  scale_fill_manual(values=c( "#c2933d", "#b54d4a")) +
  xlab("") +
  ylab("# kinases") +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 11),
        legend.position = "bottom",
        text = element_text(size = 10))

pair_p <- ggplot(kin_full, aes(x = GS, y = n_pairs, fill = benchmark)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  scale_fill_manual(values=c( "#c2933d", "#b54d4a")) +
  xlab("") +
  ylab("# kinase-patient pairs") +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 11),
        legend.position = "bottom",
        text = element_text(size = 10))


pdf(overview_meta, width = 2.7, height = 3)
kin_p
dev.off()

pdf(overview_pairs, width = 2.7, height = 3)
pair_p
dev.off()

## Most kinases ---------------------------
kin_df <- full_gs %>%
    group_by(GS, kinase) %>%
    summarise(Freq = n())

kin_order <- kin_df %>%
  group_by(kinase) %>%
  summarise(total = sum(Freq)) %>%
  arrange(desc(total)) %>%
  pull(kinase)
kin_select <- kin_order[1:30]
kin_df <- kin_df %>% filter(kinase %in% kin_select)
kin_df$kinase <- factor(kin_df$kinase, levels = rev(kin_select))

kin_p <- ggplot(kin_df, aes(x = kinase, y = Freq, fill = GS)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values=c("#E9B133", "#4E79A7")) +
  xlab("") +
  ylab("# kinase-patient pairs") +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 11),
        legend.position = "bottom",
        text = element_text(size = 10))


pdf(kin_overview_prot, width = 3, height = 6)
kin_p
dev.off()

kin_df <- full_act_gs %>%
    group_by(GS, kinase) %>%
    summarise(Freq = n())

kin_order <- kin_df %>%
  group_by(kinase) %>%
  summarise(total = sum(Freq)) %>%
  arrange(desc(total)) %>%
  pull(kinase)
kin_select <- kin_order[1:30]
kin_df <- kin_df %>% filter(kinase %in% kin_select)
kin_df$kinase <- factor(kin_df$kinase, levels = rev(kin_select))

kin_p <- ggplot(kin_df, aes(x = kinase, y = Freq, fill = GS)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values=c("#E9B133", "#4E79A7")) +
  xlab("") +
  ylab("# kinase-patient pairs") +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 11),
        legend.position = "bottom",
        text = element_text(size = 10))


pdf(kin_overview_act, width = 3, height = 6)
kin_p
dev.off()

## Jaccard gold standards ---------------------------
kinases_full <- unique(c(full_act_gs$kinase, full_gs$kinase))

overlap <- map_dfr(kinases_full, function(kin_id){
    gs_act <- full_act_gs %>% filter(kinase == kin_id) %>% pull(sample)
    gs_prot <- full_gs %>% filter(kinase == kin_id) %>% pull(sample)

    if (length(gs_prot) == 0 | length(gs_act) == 0){
        ov_coef <- NA
    } else {
        ov_coef <- length(intersect(gs_act, gs_prot))/min(c(length(gs_act), length(gs_prot)))
    }

    data.frame(kinase = kin_id, overlap_coef = ov_coef)
})
mean(overlap$overlap_coef, na.rm = T)

