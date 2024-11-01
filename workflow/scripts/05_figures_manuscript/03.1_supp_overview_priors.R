if(exists("snakemake")){
  prior_files <- snakemake@input$prior_files
  resource_class <- snakemake@input$class_kin
  upset_kin <- snakemake@output$upset
  upset_edge <- snakemake@output$upsetEdge
  jaccard_pdf <- snakemake@output$pps
  kintype_pdf <- snakemake@output$kin_type
  regulonsize_pdf <- snakemake@output$reg
}else{
  prior_files <- list.files("results/00_prior", pattern = "tsv", full.names = T)
  resource_class <- "resources/kinase_class.csv"
  prior_files <- prior_files[c(1,2,3,4,6,7,9,10)]
  jaccard_pdf <- "results/manuscript_figures/figure_3/supp/jaccard.pdf"
  kintype_pdf <- "results/manuscript_figures/figure_3/supp/kinase_type.pdf"
  upset_kin <- "results/manuscript_figures/figure_3/supp/upset_kin.pdf"
  upset_edge <- "results/manuscript_figures/figure_3/supp/upset_edge.pdf"
  regulonsize_pdf <- "results/manuscript_figures/figure_3/supp/regulon_size.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(pheatmap)
library(corrplot)
library(ggpubr)
library(ComplexHeatmap)

## Compare coverage ------------------
prior <- map(prior_files, function(file){read_tsv(file, col_types = cols()) %>%
    dplyr::filter(mor == 1) %>%
    filter(!str_detect(source, "-family")) %>%
    filter(!str_detect(source, "-subfamily"))})
tmp <- data.frame(name = str_remove(str_remove(prior_files, "results/00_prior/"), ".tsv")) %>%
  mutate(PKN = recode(name,
                      "iKiPdb" = "iKiP-DB",
                      "GPS" = "GPS gold",
                      "omnipath" = "OmniPath",
                      "networkin" = "NetworKIN",
                      "phosphositeplus" = "PhosphoSitePlus",
                      "ptmsigdb" = "PTMsigDB",
                      "combined" = "extended combined",
                      "GSknown" = "curated combined"))
names(prior) <- tmp$PKN

## UpSetPlots ------------------
kin_list <- map(prior, function(df){
  df %>% pull(source) %>% unique()
})
kin_comb <- make_comb_mat(kin_list)

pdf(upset_kin, width = 5, height = 3)
UpSet(kin_comb, pt_size = unit(2, "mm"), lwd = 1)
dev.off()

edge_list <- map(prior, function(df){
  df %>% mutate(edge = paste0(source, ":", target)) %>% pull(edge) %>% unique()
})
edge_comb <- make_comb_mat(edge_list)

pdf(upset_edge, width = 5.5, height = 3)
UpSet(edge_comb, pt_size = unit(2, "mm"), lwd = 1)
dev.off()

## Jaccard index ------------------
kinases <- map(prior, 1) %>%
  unlist %>%
  unique

jaccard_df <- map_dfr(kinases, function(kin_idx){
  target_list <- map(prior, function(df){
    df %>% filter(source == kin_idx) %>% pull(target)
  })

  dist <- unlist(lapply(combn(target_list, 2, simplify = FALSE), function(x) {
    length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]])) }))

  res <- cbind(t(combn(names(prior),2)), dist) %>%
    as.data.frame() %>%
    dplyr::rename("prior_1" = V1, "prior_2" = V2, "jaccard" = dist) %>%
    add_column(kinase = kin_idx)
  res2 <- rbind(res, res %>%
                  dplyr::rename("prior2_tmp" = prior_1, "prior1_tmp" = prior_2) %>%
                  dplyr::rename("prior_2" = prior2_tmp, "prior_1" = prior1_tmp)) %>%
    mutate(jaccard = as.numeric(jaccard))

  rbind(res2, data.frame(prior_1 = unique(res2$prior_1),
                         prior_2 = unique(res2$prior_1),
                         jaccard = 1,
                         kinase = unique(res2$kinase))) %>%
    arrange(prior_1,prior_2)
})

jaccard_m <- jaccard_df %>%
  group_by(prior_1, prior_2) %>%
  summarise(mean_jaccard = mean(jaccard, na.rm = T)) %>%
  pivot_wider(names_from = prior_2, values_from = mean_jaccard) %>%
  column_to_rownames("prior_1")

pdf(jaccard_pdf, height = 3, width = 3)
corrplot(as.matrix(jaccard_m), is.corr = F, method = 'color', col = COL1('Purples', 100)[20:100],
         addgrid.col = 'white',
         order = 'hclust',
         tl.col = 'black',
         type = 'lower', cl.pos = 'r'
)
dev.off()

jaccard_m
## Tyrosine versus Serin/Thr
class_kin <- "prior"
kinase_type_df <- read_csv(resource_class, col_types = cols()) %>%
  mutate(resource = recode(resource,
                      "iKiPdb" = "iKiP-DB",
                      "GPS" = "GPS gold",
                      "omnipath" = "OmniPath",
                      "networkin" = "NetworKIN",
                      "phosphositeplus" = "PhosphoSitePlus",
                      "ptmsigdb" = "PTMsigDB",
                      "combined" = "extended combined",
                      "GSknown" = "curated combined")) %>%
  filter(resource %in% names(prior))

if (class_kin == "prior"){
  kin_type <- kinase_type_df %>%
    group_by(resource, class) %>%
    summarise(n = n())
  kin_type$kinase <- factor(kin_type$class, levels = c("Serine/Threonine", "Tyrosine", "Histidine", "Dual-specificity"))
  kin_type$resource <- factor(kin_type$resource, levels = unique(kin_type %>% arrange(desc(n)) %>% pull(resource)))
} else {
  kin_type <- kinase_type_df %>%
    group_by(resource, kinase) %>%
    summarise(n = n())
  kin_type$kinase <- factor(kin_type$kinase, levels = c("Serine/Threonine", "Tyrosine", "Histidine", "Ambiguous"))
  kin_type$resource <- factor(kin_type$resource, levels = unique(kin_type %>% arrange(desc(n)) %>% pull(resource)))
}

cbPalette <- c("#56B4E9", "#009E73",  "#CC79A7", "#E69F00")

kintype_p <- ggplot(kin_type, aes(x = resource, y = n, fill = kinase)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal() +
  scale_fill_manual(values=cbPalette) +
  xlab("") + ylab("number of kinases") + guides(fill=guide_legend(title="kinase class")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 10),
        legend.key.size = unit(0.2, 'cm'))

pdf(kintype_pdf, height = 3, width = 4.2)
kintype_p
dev.off()

## Mean set size
set_size <- map_dfr(names(prior), function(prior_idx){
  prior[[prior_idx]] %>%
    group_by(source) %>%
    summarise(n_targets = n()) %>%
    add_column(prior = prior_idx)
})

set_size %>%
  group_by(prior) %>%
  summarise(mean_n = median(n_targets))

regulonsize_p <- ggplot(set_size, aes(x = prior, y = n_targets)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("regulon size") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 10),
        legend.key.size = unit(0.2, 'cm'))

pdf(regulonsize_pdf, height = 3, width = 2.8)
regulonsize_p
dev.off()

