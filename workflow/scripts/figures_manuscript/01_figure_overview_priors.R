if(exists("snakemake")){
  prior_files <- snakemake@input$prior_files
  coverage_pdf <- snakemake@output$kin
  kinase_pdf <- snakemake@output$kin_heat
  edge_pdf <- snakemake@output$edges
  jaccard_pdf <- snakemake@output$pps
}else{
  prior_files <- list.files("results/prior", pattern = "tsv", full.names = T)
  coverage_pdf <- "results/manuscript_figures/figure_1/coverage_merged.pdf"
  kinase_pdf <- "results/manuscript_figures/figure_1/kinase_overview.pdf"
  edge_pdf <- "results/manuscript_figures/figure_1/edge_overview.pdf"
  jaccard_pdf <- "results/manuscript_figures/figure_1/jaccard.pdf"
}

## Libraries ---------------------------
library(tidyverse)
library(pheatmap)
library(corrplot)
library(ggpubr)

## Compare coverage ------------------
prior <- map(prior_files, function(file){read_tsv(file, col_types = cols()) %>%
    dplyr::filter(mor == 1) %>%
    filter(!str_detect(source, "-family")) %>%
    filter(!str_detect(source, "-subfamily"))})
names(prior) <- str_remove(str_remove(prior_files, "results/prior/"), ".tsv")

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

coverage <- coverage %>% arrange(desc(value))
PKN_order <- coverage %>%
  dplyr::filter(class == "kinase") %>%
  dplyr::select(PKN, type, value) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  arrange(desc(`all kinases`), desc(`kinases with \nat least 5 targets`)) %>%
  pull(PKN)
coverage$PKN <- factor(coverage$PKN, levels = PKN_order)

text_size <- 10

pps_p <- ggplot(coverage %>% filter(class == "pps")) +
  aes(x = PKN, y = value) +
  geom_bar(stat="identity", color="black", position=position_dodge(), fill = "#47AD7A", width=0.7)+
  theme_minimal() + xlab("") + ylab("") +
  ggtitle("unique peptides") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = text_size))

edges_p <- ggplot(coverage %>% filter(class == "edges")) +
  aes(x = PKN, y = value) +
  geom_bar(stat="identity", color="black", position=position_dodge(), fill = "#AD477A", width=0.7)+
  theme_minimal() + xlab("") + ylab("") +
  ggtitle("kinase-peptide links") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = text_size)) +
  scale_y_continuous(expand = c(0, 0))

kin_p <- ggplot(data=coverage %>% filter(class == "kinase"), aes(x=PKN, y=value, fill=type)) +
  geom_bar(stat="identity",color="black", position=position_dodge(), width=0.7) +
  scale_fill_manual(labels = c("all", "> 4 targets"), values=c('#477AA3','#97CAF3')) +
  theme_minimal() + xlab("") + ylab("") + guides(fill=guide_legend(title="kinases")) +
  ggtitle("unique kinases") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = text_size),
        legend.position = "bottom",
        legend.key.size = unit(0.2, 'cm')) +
  scale_y_continuous(expand = c(0, 0))

figure_merged <- ggarrange(kin_p + coord_flip(), edges_p + coord_flip(),
          labels = NULL,
          ncol = 1,
          common.legend = TRUE, legend = "top",
          align = "hv",
          font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

pdf(file=coverage_pdf, height = 5, width = 3)
plot(figure_merged)
dev.off()

## Kinase overview ------------------
resource_df <- map_dfr(names(prior), function(x){
  df <- prior[[x]]
  df %>%
    add_column(resource = x) %>%
    mutate(edge = paste(source, target, sep = "_"))
})

kinase_m <- resource_df %>%
  dplyr::select(source,resource) %>%
  add_column(present = 1) %>%
  distinct() %>%
  pivot_wider(names_from = source, values_from = present, values_fill = 0) %>%
  column_to_rownames("resource")

sum(colSums(kinase_m) > 1)/ncol(kinase_m)
kinase_m[colSums(kinase_m) == 1] %>%
  rowSums()
order_kin <- colSums(kinase_m) %>%
  order()
kinase_m <- kinase_m[c(PKN_order), rev(order_kin)]

coverage_df <- pheatmap(kinase_m, cluster_rows = F,
                        cluster_cols = F, show_colnames = F,
                        treeheight_row = 0, color = c("white", "#477AA3"))

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

pdf(edge_pdf, height = 1.5, width = 7)
print(coverage_edge_df)
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

## Tyrosine versus Serin/Thr
kinase_type_df <- map_dfr(names(prior), function(x){
  df <- prior[[x]]
  tmp <- df %>%
    mutate(kinase = case_when(
      str_detect(position, "T") ~ "Threonine",
      str_detect(position, "S") ~ "Serine",
      str_detect(position, "Y") ~ "Tyrosin",
      str_detect(position, "H") ~ "Histidine"
    ))

  if ("Histidine" %in% tmp$kinase){
    overview_kin_type <- table(tmp$source, tmp$kinase) %>%
      as.data.frame() %>%
      pivot_wider(names_from = "Var2", values_from = "Freq") %>%
      mutate(kinase = case_when(
        (Serine > 0 | Threonine > 0) & Tyrosin == 0 & Histidine == 0 ~ "Serine/Threonine",
        (Serine == 0 & Threonine == 0) & Tyrosin == 0 & Histidine > 0 ~ "Histidine",
        (Serine == 0 & Threonine == 0) & Histidine == 0 & Tyrosin > 0 ~ "Tyrosine",
        ((Serine > 0 | Threonine > 0) & Tyrosin > 0) | (Histidine > 0 & Tyrosin > 0) | ((Serine > 0 | Threonine > 0) & Histidine > 0)~ "Ambiguous"
      ))
    overview_kin_type$n_targets <- rowSums(overview_kin_type[2:5])
    overview_kin_type$resource <- x
  } else{
    overview_kin_type <- table(tmp$source, tmp$kinase) %>%
      as.data.frame() %>%
      pivot_wider(names_from = "Var2", values_from = "Freq") %>%
      mutate(kinase = case_when(
        (Serine > 0 | Threonine > 0) & Tyrosin == 0 ~ "Serine/Threonine",
        (Serine == 0 & Threonine == 0) & Tyrosin > 0 ~ "Tyrosine",
        ((Serine > 0 | Threonine > 0) & Tyrosin > 0)~ "Ambiguous"
      ))
    overview_kin_type$n_targets <- rowSums(overview_kin_type[2:4])
    overview_kin_type$resource <- x
    overview_kin_type <- overview_kin_type %>%
      add_column(Histidine = 0, .before = "Serine")
  }
  overview_kin_type
})

kin_type <- kinase_type_df %>%
  group_by(resource, kinase) %>%
  summarise(n = n())
kin_type$kinase <- factor(kin_type$kinase, levels = c("Serine/Threonine", "Tyrosine", "Histidine", "Ambiguous"))
kin_type$resource <- factor(kin_type$resource, levels = unique(kin_type %>% arrange(desc(n)) %>% pull(resource)))

cbPalette <- c("#56B4E9", "#009E73",  "#CC79A7", "#E69F00")
ggplot(kin_type, aes(x = resource, y = n, fill = kinase)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal() +
  scale_fill_manual(values=cbPalette)
