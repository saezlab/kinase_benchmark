if(exists("snakemake")){
  prior_files <- snakemake@input$prior_files
  coverage_pdf <- snakemake@output$kin
  edges_coverage_pdf <- snakemake@output$edges
  pps_coverage_pdf <- snakemake@output$pps
  kinase_pdf <- snakemake@output$kin_heat
  height <- snakemake@params$plot_height
  width <- snakemake@params$plot_width
}else{
  prior_files <- list.files("results/prior", pattern = "tsv", recursive = T, full.names = T)
  coverage_pdf <- "results/comparison/plots/coverage_kinases.pdf"
  edges_coverage_pdf <- "results/comparison/plots/coverage_edges.pdf"
  pps_coverage_pdf <- "results/comparison/plots/coverage_pps.pdf"
  kinase_pdf <- "results/comparison/plots/coverage_pps.pdf"
  height <- 4
  width <- 3.2
}
height <- as.numeric(height)
width <- as.numeric(width)

## Libraries ---------------------------
library(tidyverse)
library(pheatmap)

## Compare coverage ------------------
prior <- map(prior_files, function(file){read_tsv(file, col_types = cols())})
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
        text = element_text(size = text_size))

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

pdf(file=coverage_pdf, height = height*0.95, width = width*0.6)
plot(kin_p)
dev.off()

pdf(file=edges_coverage_pdf, height = height*0.6, width = width*0.5)
plot(edges_p)
dev.off()

pdf(file=pps_coverage_pdf, height = height*0.54, width = width*0.5)
plot(pps_p)
dev.off()

## Kinase overview ------------------
resource_df <- map_dfr(names(prior), function(x){
  df <- prior[[x]]
  df %>%
    add_column(resource = x)
})

kinase_m <- resource_df %>%
  select(source,resource) %>%
  add_column(present = 1) %>%
  distinct() %>%
  pivot_wider(names_from = source, values_from = present, values_fill = 0) %>%
  column_to_rownames("resource")

order_kin <- colSums(kinase_m) %>%
  order()
kinase_m <- kinase_m[c(PKN_order), order_kin]

coverage_df <- pheatmap(kinase_m, cluster_rows = F,
         cluster_cols = F, show_colnames = F,
         treeheight_row = 0, color = c("white", "#477AA3"))

pdf(kinase_pdf, height = height*0.5, width = width*0.65)
print(coverage_df)
dev.off()

## Mean Jaccard ------------------
kinases <- map(prior, 1) %>%
  unlist %>%
  unique

map_dfr(kinases, function(kin_idx){
  target_list <- map(prior, function(df){
    df %>% filter(source == kin_idx) %>% pull(target)
  })


})
