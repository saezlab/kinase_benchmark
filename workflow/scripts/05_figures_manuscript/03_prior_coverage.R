if(exists("snakemake")){
  prior_files <- snakemake@input$prior_files
  kin_pdf <- snakemake@output$kin
  kinase_pdf <- snakemake@output$kin_heat
  edg_pdf <- snakemake@output$edges
  edge_pdf <- snakemake@output$edges_heat
}else{
  prior_files <- list.files("results/00_prior", pattern = "tsv", full.names = T)
  prior_files <- prior_files[c(1:4,6,7,9,10)]
  kin_pdf <- "results/manuscript_figures/figure_3/coverage_kin.pdf"
  kinase_pdf <- "results/manuscript_figures/figure_3/kinase_overview.pdf"
  edg_pdf <- "results/manuscript_figures/figure_3/coverage_edge.pdf"
  edge_pdf <- "results/manuscript_figures/figure_3/edge_overview.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(pheatmap)
library(corrplot)
library(ggpubr)
library(ComplexHeatmap)
library(patchwork)

## Compare coverage ------------------
prior <- map(prior_files, function(file){read_tsv(file, col_types = cols()) %>%
    dplyr::filter(mor == 1) %>%
    filter(!str_detect(source, "-family")) %>%
    filter(!str_detect(source, "-subfamily"))})
names(prior) <- str_remove(str_remove(prior_files, "results/00_prior/"), ".tsv")

coverage <- map_dfr(names(prior), function(PKN_idx){
  PKN <- prior[[PKN_idx]]
  kinases <- table(PKN$source) %>% as.data.frame()

  data.frame(PKN = PKN_idx,
             value = c(nrow(kinases),
                       nrow(kinases %>% dplyr::filter(Freq >= 5)),
                       length(unique(PKN$target)),
                       nrow(PKN)),
             type = c("all kinases",
                      "kinases with \nat least 5 targets",
                      "none",
                      "none"),
             class = c("kinase", "kinase", "pps", "edges"))
})

coverage <- coverage %>% arrange(desc(value)) %>%
  mutate(PKN = recode(PKN,
                      "iKiPdb" = "iKiP-DB",
                      "GPS" = "GPS gold",
                      "omnipath" = "OmniPath",
                      "networkin" = "NetworKIN",
                      "phosphositeplus" = "PhosphoSitePlus",
                      "ptmsigdb" = "PTMsigDB",
                      "combined" = "extended combined",
                      "GSknown" = "curated combined"))
coverage %>%
  filter(type == "all kinases" & class == "kinase")

coverage %>%
  filter(class == "edges")

PKN_order <- coverage %>%
  dplyr::filter(class == "kinase") %>%
  dplyr::select(PKN, type, value) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  arrange(desc(`all kinases`), desc(`kinases with \nat least 5 targets`)) %>%
  pull(PKN)
coverage$PKN <- factor(coverage$PKN, levels = PKN_order)

text_size <- 10

kin_p <- ggplot(data=coverage %>% filter(class == "kinase"), aes(x=value, y=PKN, fill=type)) +
  geom_bar(stat="identity",color="black", position=position_dodge(), width=0.7) +
  scale_fill_manual(labels = c("all", "> 4 targets"), values=c('#477AA3','#97CAF3'))+
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, family = "Helvetica", size = 10),
        text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
        axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
        legend.position = "bottom",
        legend.key.size = unit(0.2, 'cm')) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_minimal()

edges_p <- ggplot(coverage %>% filter(class == "edges")) +
  aes(x = value, y = PKN) +
  geom_bar(stat="identity", color="black", position=position_dodge(), fill = "#AD477A", width=0.7)+
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(family = "Helvetica", size = 11), # Set label font to Helvetica size 9
        axis.title.y = element_text(family = "Helvetica", size = 10), # Set y-axis label size to 10
        axis.title.x = element_text(family = "Helvetica", size = 10)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_minimal()


## Kinase overview ------------------
resource_df <- map_dfr(names(prior), function(x){
  df <- prior[[x]]
  df %>%
    add_column(resource = x) %>%
    mutate(edge = paste(source, target, sep = "_"))
}) %>%
  mutate(resource = recode(resource,
                           "iKiPdb" = "iKiP-DB",
                           "GPS" = "GPS gold",
                           "omnipath" = "OmniPath",
                           "networkin" = "NetworKIN",
                           "phosphositeplus" = "PhosphoSitePlus",
                           "ptmsigdb" = "PTMsigDB",
                           "combined" = "extended combined",
                           "GSknown" = "curated combined"))

kinase_m <- resource_df %>%
  dplyr::select(source,resource) %>%
  add_column(present = 1) %>%
  distinct() %>%
  pivot_wider(names_from = source, values_from = present, values_fill = 0) %>%
  column_to_rownames("resource")

kinase_tmp <- kinase_m[!str_detect(rownames(kinase_m), "Curated"), ]
sum(colSums(kinase_tmp) > 1)/ncol(kinase_tmp)
kinase_tmp[colSums(kinase_tmp) == 1] %>%
  rowSums()
order_kin <- colSums(kinase_tmp) %>%
  order()
kinase_m <- kinase_m[c(PKN_order), rev(order_kin)]

coverage_df <- pheatmap(kinase_m, cluster_rows = F,
                        cluster_cols = F, show_colnames = F,
                        treeheight_row = 0, color = c("white", "#477AA3"))


pdf(file=kin_pdf, height = 1.91, width = 3.5)
kin_p
dev.off()

pdf(file=edg_pdf, height = 1.91, width = 2.3)
edges_p
dev.off()

pdf(kinase_pdf, height = 1.5, width = 3.7)
print(coverage_df)
dev.off()


## Edge overview ------------------
edge_m <- resource_df %>%
  dplyr::select(edge,resource) %>%
  add_column(present = 1) %>%
  distinct() %>%
  pivot_wider(names_from = edge, values_from = present, values_fill = 0) %>%
  column_to_rownames("resource")

order_edge <- colSums(edge_m) %>%
  order()
edge_m <- edge_m[c(PKN_order), rev(order_edge)]
sum(colSums(edge_m) > 1)/ncol(edge_m)
edge_m[colSums(edge_m) == 1] %>%
  rowSums()

coverage_edge_df <- pheatmap(edge_m, cluster_rows = F,
                             cluster_cols = F, show_colnames = F,
                             treeheight_row = 0, color = c("white", "#AD477A"))

pdf(edge_pdf, height = 1.5, width = 10)
print(coverage_edge_df)
dev.off()
