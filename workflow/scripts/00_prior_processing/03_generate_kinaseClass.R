#'

## Snakemake ---------------------------
if(exists("snakemake")){
  prior_files <- snakemake@input$prior
  literature_file <- snakemake@input$lit
  output_file <- snakemake@output$out
}else{
  prior_files <- list.files("results/00_prior", pattern = "tsv", full.names = T)
  literature_file <- "data/misc/reference_human_kinases.csv"
  kin_map_file <- "resources/kinase_mapping.tsv"
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

## Add literature information ---------------------------
literature_kinase <- read_csv(literature_file, col_types = cols()) %>% 
  dplyr::rename("Var1" = gene)
dual_spec_kins <- c('CLK1', 'CLK2', 'CLK3', 'CLK4', 'DYRK1A', 'DYRK1B', 'DYRK2',
                    'DYRK3', 'DYRK4', 'MAP2K1', 'MAP2K2', 'MAP2K3', 'MAP2K4',
                    'MAP2K5', 'MAP2K6', 'MAP2K7', 'TESK1', 'TESK2', 'SgK496', 'TTK')

kinase_mapping <- read_tsv(kin_map_file, col_types = cols()) %>%
  dplyr::rename("Var1" = external_gene_name, "uniprot" = uniprotswissprot)

kinase_phospho <- kinase_type_df %>%
  select(Var1, kinase) %>%
  distinct()

# Add Uniprot ID for mapping
kinase_uniprot <- left_join(kinase_phospho, kinase_mapping %>% dplyr::select(Var1, uniprot), by = c("Var1"), relationship = "many-to-many")

# map based on UniProt
kinase_overview <- kinase_uniprot %>%
  left_join(literature_kinase %>% dplyr::select(uniprot, group) %>% rename("group_uniprot" = group), by = c("uniprot"), relationship = "many-to-many")

# map based on Gene Name
kinase_overview <- kinase_overview %>%
  left_join(literature_kinase %>% dplyr::select(Var1, group) %>% rename("group_gene" = group), by = c("Var1"), relationship = "many-to-many")

# Generate final table
kinase_overview <- kinase_overview %>%
  mutate(
    group = coalesce(group_uniprot, group_gene)
  ) %>%
  select(-group_uniprot, -group_gene) 

kinase_overview <- kinase_overview %>%
  mutate(class = case_when( # Assign class based on literature
    is.na(group) ~ "Unknown",
    group %in% c("TK", "RTK") ~ "Tyrosine",
    !group %in% c("TK", "RTK") ~ "Serine/Threonine"
  )) %>%
  mutate(class = case_when( # Add Dual-specificity kinases
    Var1 %in% dual_spec_kins ~ "Dual-specificity",
    TRUE ~ class
  )) %>%
  mutate(class = case_when( # Assign based on databases if kinase class not in list
    class == "Unknown" ~ kinase,
    TRUE ~ class
  )) %>%
  mutate(class = case_when( # Manually check kinases with ambiguous classification
    Var1 %in% c("BLVRA", "EEF2K") ~ "Serine/Threonine",
    Var1 %in% c("JAK2") ~ "Tyrosine",
    Var1 %in% c("NME1") ~ "Histidine",
    TRUE ~ class
  )) %>%
  filter(class != "Ambiguous") %>% # Remove non-protein kinases
  dplyr::select(Var1, class, group) %>%
  distinct()

# Add final classification to the overview
kinase_type_df <- left_join(kinase_type_df, kinase_overview, by = c("Var1"), relationship = "many-to-many") %>%
  rename("source" = Var1)

kinase_type_df <- kinase_type_df[,c("source", "class", "group", "kinase", "n_targets", "resource", "Serine", "Threonine", "Tyrosin", "Histidine")] #rearrange column order

## Save data ---------------------------
write_csv(kinase_type_df, output_file)
