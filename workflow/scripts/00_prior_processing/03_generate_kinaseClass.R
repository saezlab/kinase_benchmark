#'

## Snakemake ---------------------------
if(exists("snakemake")){
  prior_files <- snakemake@input$prior
  output_file <- snakemake@output$out
}else{
  prior_files <- list.files("results/00_prior", pattern = "tsv", full.names = T)
  output_file <- "resources/kinase_class.csv"
}

## Libraries ---------------------------
library(tidyverse)

## Load data ---------------------------
prior <- map(prior_files, function(file){read_tsv(file, col_types = cols()) %>%
    dplyr::filter(mor == 1) %>%
    filter(!str_detect(source, "-family")) %>%
    filter(!str_detect(source, "-subfamily"))})
names(prior) <- str_remove(str_remove(prior_files, "results/00_prior/"), ".tsv")

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


write_csv(kinase_type_df, output_file)
