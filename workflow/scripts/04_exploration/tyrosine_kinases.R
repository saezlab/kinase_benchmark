#'

## Snakemake ---------------------------
if(exists("snakemake")){
  bench_files <- snakemake@input$bench
  rank_files <- snakemake@input$rank
  k_phit <- snakemake@params$k_phit
  performance_plot <- snakemake@output$plot
}else{
  cptac_data <- list.files("results/01_processed_data/cptac/data/original", full.names = T)
  hernandez_data <-  "results/01_processed_data/hernandez/data/benchmark_data.csv"                  
  hijazi_data <- "results/01_processed_data/hijazi/data/benchmark_data.csv"
  performance_plot <- "results/04_exploration/merged/benchmark/plots/performance_subset.pdf"
  k_phit <- 10
}

## Libraries ---------------------------
library(tidyverse)

## Load data ---------------------------
hernandez <- read_csv(hernandez_data, col_types = cols())
hernandez_df <- pivot_longer(hernandez, !ID, names_to = "exp", values_to = "logFc") %>%
    filter(!is.na(logFc)) %>%
    mutate(pps = map_chr(str_split(ID,"\\|"),2)) %>%
    mutate(aa = case_when(
        str_detect(pps, "T") ~ "Threonine",
        str_detect(pps, "S") ~ "Serine",
        str_detect(pps, "Y") ~ "Tyrosine",
        str_detect(ID, "H") ~ "Histidine"
    ))

hernandez_df %>% group_by(aa) %>% summarise(perc = (n()/nrow(hernandez_df))*100, total = n())

hijazi <- read_csv(hijazi_data, col_types = cols())
hijazi_df <- pivot_longer(hijazi, !ID, names_to = "exp", values_to = "logFc") %>%
    filter(!is.na(logFc)) %>%
    mutate(pps = map_chr(str_split(ID,"\\|"),3)) %>%
    mutate(aa = case_when(
        str_detect(ID, "T") ~ "Threonine",
        str_detect(ID, "S") ~ "Serine",
        str_detect(ID, "Y") ~ "Tyrosine",
        str_detect(ID, "H") ~ "Histidine"
    ))

hijazi_df %>% group_by(aa) %>% summarise(perc = (n()/nrow(hernandez_df))*100, total = n())

perturb_df <- rbind(hijazi_df, hernandez_df)
perturb_df %>% group_by(aa) %>% summarise(perc = (n()/nrow(hernandez_df))*100, total = n())

## cptac ---------------------------
cptac_df <- map_dfr(cptac_data, function(file){
    df <- read_csv(file, col_types = cols())

    pivot_longer(df, !ID, names_to = "exp", values_to = "logFc") %>%
        filter(!is.na(logFc)) %>%
        mutate(pps = map_chr(str_split(ID,"\\|"),3)) %>%
        mutate(aa = case_when(
            str_detect(pps, "T") ~ "Threonine",
            str_detect(pps, "S") ~ "Serine",
            str_detect(pps, "Y") ~ "Tyrosine",
            str_detect(ID, "H") ~ "Histidine"
        ))
})

cptac_df %>% group_by(aa) %>% summarise(perc = (n()/nrow(cptac_df))*100, total = n())
