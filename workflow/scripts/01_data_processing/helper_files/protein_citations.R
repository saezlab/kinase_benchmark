#'

## Snakemake ---------------------------
if(exists("snakemake")){
  gene_citations <- snakemake@input$gene
  info_file <- snakemake@input$info
  output_file <- snakemake@output$out
}else{
  gene_citations <- "data/gene2pubmed"
  info_file <- "data/gene_info"
  output_file <- "resources/protein_citations.csv"
}

## Libraries ---------------------------
library(tidyverse)

## ---------------------------
## Compare to number of citations
overview_citations <- read_table(gene_citations, col_types = cols())

gene_info <- read_table(info_file, col_types = cols())

# Filter for human kinases
citations_human <- overview_citations %>%
  filter(`#tax_id` == "9606") %>%
  group_by(GeneID) %>%
  summarise(citations = n())
rm(overview_citations)

gene_info_human <- gene_info %>%
  filter(`#tax_id` == "9606")
rm(gene_info)

merged_df <- left_join(citations_human, gene_info_human, by = "GeneID")
head(merged_df)

merged_df %>% dplyr::select(citations, Symbol, Synonyms) %>%
  write_csv(output_file)
