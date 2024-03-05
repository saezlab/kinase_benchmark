if(exists("snakemake")){
  ppsp_file <- snakemake@input$ppsp
  output_file <- snakemake@output$tsv
}else{
  ppsp_file <- "data/kinase_libraries/phosformer/phosformer_results_15.csv"
  output_file <- "results/prior/raw/phosformer15.tsv"
}

## Libraries ---------------------------
library(tidyverse)
library(org.Hs.eg.db)

## Construct kinase-substrate interaction network ---------------------------
net <- read_csv(ppsp_file, col_types = cols())

# Load kinase symvols
kins <- read_csv('https://raw.githubusercontent.com/esbgkannan/phosformer/main/data/reference_human_kinases.csv')
kins_mapping <- kins %>%
  dplyr::select(gene, uniprot) %>%
  distinct() %>%
  dplyr::mutate(genesymbol = mapIds(x = org.Hs.eg.db, keys = uniprot, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first"))

# Reformat dataframe
net_long <- left_join(net, kins_mapping, by = "gene", relationship = "many-to-many") %>%
  mutate(source = case_when(
    is.na(genesymbol) ~ gene,
    !is.na(genesymbol) ~ genesymbol
  )) %>%
  mutate(target = paste0(Protein, "_", Aminoacid, Position)) %>%
  dplyr::rename(target_protein = "Protein") %>%
  mutate(position = paste0(Aminoacid, Position)) %>%
  mutate(sequence = str_replace(`11-mer`, pattern = "-", replacement = "_")) %>%
  add_column(mor = 1, .before = "sequence") %>%
  dplyr::select(source, target, target_protein, position, mor, sequence)

net_long <- net_long %>%
  mutate(source = recode(source,
                         "PDPK2" = "PDPK2P",
                         "KS6R" = "RSKR"))

## Save processed phosphositeplus
write_tsv(net_long, output_file)
