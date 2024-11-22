if(exists("snakemake")){
  meta_file <- snakemake@input$meta
  meta_out <- snakemake@output$out
}else{
  meta_file <- "results/01_processed_data/hijaziDiscoverX/data/benchmark_metadata.csv"
  meta_out <- "results/manuscript_figures/supp_files/hijazi_meta.csv"
}


## Libraries ---------------------------
library(tidyverse)

## Format data ---------------------------
meta <- read_csv(meta_file, col_types = cols())

meta_pretty <- meta %>% dplyr::select(cell_line, drug, target, sign) %>%
    mutate(sign = case_when(
        sign == -1 ~ "inhibitor", 
        sign == 1 ~ "activator"
    )) %>%
    rename("drug type" = "sign", "cell line" = "cell_line") 

write_csv(meta_pretty, meta_out)
