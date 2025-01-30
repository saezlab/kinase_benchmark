#'

## Snakemake ---------------------------
if(exists("snakemake")){
  raw_24362263 <- snakemake@input$id24362263
  raw_17389395 <- snakemake@input$id17389395
  raw_17016520 <- snakemake@input$id17016520
  raw_25147952 <- snakemake@input$id25147952
  raw_24804581 <- snakemake@input$id24804581
  raw_30257219 <- snakemake@input$id30257219
  meta_overview <- snakemake@input$meta
  bench_output <- snakemake@output$mat
  meta_prior_output <- snakemake@output$meta_out
}else{
  raw_24362263 <- "data/datasets/tyrosine/NIHMS551505-supplement-2.xlsx"
  raw_17389395 <- "data/datasets/tyrosine/pnas_0608638104_08638Table3.xls"
  raw_17016520 <- "data/datasets/tyrosine/msb4100094-s3.xls"
  raw_25147952 <- "data/datasets/tyrosine/TableS1.xlsx"
  raw_24804581 <- "data/datasets/tyrosine/cb500116c_si_004.xls"
  raw_30257219 <- "data/datasets/tyrosine/processed.xlsx"
  meta_overview <- "data/datasets/tyrosine/Overview.csv"
  bench_output <- "results/01_processed_data/tyrosine/data/benchmark_data.csv"
  meta_prior_output <- "results/01_processed_data/tyrosine/data/benchmark_metadata.csv"
}

## Libraries ---------------------------
library(tidyverse)
library(readxl)
library(biomaRt)

## Prepare data ---------------------------
# PMID: 24362263, Dasatinib_inhibition_Asmussen, select 20 minute exposure of Dasatinib before washout to observe downregualtion of ABL1
proc_24362263 <- read_excel(raw_24362263, skip = 2) %>%
  dplyr::select(`UniProt Protein Name`, `Phosphorylation site`, Peptide, `Acc .`, `EOE (L/H)`) %>%
  dplyr::rename("PMID_24362263" = `EOE (L/H)`) %>%
  mutate(PMID_24362263 = log2(as.numeric(PMID_24362263))) %>%
  filter(!is.na(PMID_24362263)) %>%
  filter(!is.infinite(PMID_24362263)) %>%
  dplyr::mutate(ID = paste0(`UniProt Protein Name`,"_", `Phosphorylation site`, "|", `UniProt Protein Name`,"|", `Phosphorylation site`)) %>%
  dplyr::select(ID, PMID_24362263)

# PMID: 17389395, Human mammary epithelial cells stimulated with EGF for 1, 4, 8, 16 minutes (#selection based on KSTAR)
proc_17389395 <- read_excel(raw_17389395) %>%
  dplyr::select(Name, `Matched Sequence`,Abreviation, Residue, Number, contains("QSAV")) %>%
  dplyr::filter(!str_detect(Residue, "\\/")) %>% #remove double phosphorylation
  dplyr::filter(!str_detect(`QSAV 0`, "ot found")) %>%
  dplyr::mutate(across(contains("QSAV"), ~ as.numeric(.))) %>%
  dplyr::mutate(across(contains("QSAV"), ~ . / `QSAV 0`)) %>%
  dplyr::mutate(across(contains("QSAV"), ~ log2(.)))

transl_gene <- read.csv("resources/translate_genes_17389395.csv", sep = ";") %>%
  dplyr::select(external_gene_name, Name)

proc_17389395 <- left_join(proc_17389395, transl_gene, by = "Name", relationship = "many-to-many") %>%
  mutate(Number = as.numeric(Number)) %>%
  filter(!is.na(Number))%>%
  dplyr::mutate(ID = paste0(external_gene_name,"_", Residue, Number, "|", external_gene_name,"|", Residue, Number)) %>%
  dplyr::select(ID, contains("QSAV")) %>%
  dplyr::select(-`QSAV 0`, -`QSAV 2`, -`QSAV 32`)

colnames(proc_17389395) <- str_remove(colnames(proc_17389395), " ")


# PMID: 17016520, HER2 overexpressing epithelial cells stimulated with EGF and HRG for 10 minutes
proc_17016520 <- read_excel(raw_17016520) %>%
  filter(`Run#` %in% c("24H_EGF", "24H_HRG")) %>%
  mutate(`10/5`= as.numeric(`10/5`)/ as.numeric(`0/5`)) %>%
  dplyr::select(`Run#`, `Protein Name` , Abbreviation, Residue, Number, Sequence, `10/5`) %>%
  dplyr::filter(!str_detect(Residue, "\\/")) %>% #remove double phosphorylation %>%
  filter(!is.na(`10/5`)) %>%
  filter(!is.infinite(`10/5`)) %>%
  mutate(`10/5` = log(`10/5`)) %>%
  pivot_wider(names_from = `Run#`, values_from = `10/5`)

transl_gene <- read.csv("resources/translate_genes_17016520.csv", sep = ";") %>%
  dplyr::select(external_gene_name, Name)

proc_17016520 <- left_join(proc_17016520, transl_gene, by = c("Protein Name" = "Name"), relationship = "many-to-many") %>%
  mutate(Number = as.numeric(Number)) %>%
  filter(!is.na(Number))%>%
  filter(!is.na(external_gene_name))%>%
  dplyr::mutate(ID = paste0(external_gene_name,"_", Residue, Number, "|", external_gene_name,"|", Residue, Number)) %>%
  dplyr::select(ID, contains("24")) %>%
  group_by(ID) %>%
  summarise(across(everything(), ~ if (all(is.na(.))) NA_real_ else na.omit(.)[1]), .groups = "drop") %>%
  distinct()

# PMID: 25147952, TCR stimulation for 5, 15, 30, 60 seconds
proc_25147952 <- read_excel(raw_25147952)
colnames(proc_25147952)[6:10] <- proc_25147952[1, 6:10]
proc_25147952 <- proc_25147952 %>%
  slice(-1) %>%
  dplyr::select(Position, `Gene names`, Uniprot, `Modified Sequence`, contains("sec")) %>%
  dplyr::select(-`0 sec`)%>%
  mutate(Position = as.numeric(Position)) %>%
  filter(!is.na(Position)) %>%
  mutate(aa = str_extract(`Modified Sequence`, ".(?=\\(ph\\))")) %>%
  mutate(gene = map_chr(str_split(`Gene names`, ";"), 1)) %>%
  mutate(ID = paste0(gene,"_", aa, Position, "|", gene,"|", aa, Position) ) %>%
  dplyr::select(ID, contains("sec")) %>%
  mutate(across(matches("\\d+ sec$"), as.numeric))

colnames(proc_25147952) <- str_remove(colnames(proc_25147952), " ")

# PMID: 24804581, TCR stimulation for 5, 15, 30, 60 seconds
proc_24804581 <- read_excel(raw_24804581, sheet = 2) %>%
  dplyr::filter(NumberOfPhospho == "n1") %>%
  dplyr::select(`Protein Group Accessions`, Sequence, Genes, proteinSites, contains("Log2")) %>%
  dplyr::filter(!is.na(proteinSites)) %>%
  separate_rows(`Protein Group Accessions`, sep = ";") %>%
  mutate(Accession = `Protein Group Accessions`) %>%
  dplyr::mutate(proteinSites = map_chr(str_split(proteinSites, ";"), 1)) %>%
  mutate(aa = str_extract(Sequence, ".(?=\\*)"))


# Define the BioMart dataset
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # Replace "hsapiens_gene_ensembl" for other species if needed

# Query BioMart to get gene names corresponding to protein accessions
results <- getBM(
  attributes = c("uniprotswissprot", "external_gene_name"),  # UniProt ID and gene name
  filters = "uniprotswissprot",
  values = unique(proc_24804581$Accession),
  mart = mart
)

proc_24804581 <- left_join(proc_24804581, results, by = c("Accession" = "uniprotswissprot"), relationship = "many-to-many") %>%
  filter(!is.na(external_gene_name)) %>%
  mutate(proteinSites = as.numeric(proteinSites)) %>%
  filter(!is.na(proteinSites)) %>%
  mutate(ID = paste0(external_gene_name,"_", aa, proteinSites, "|", external_gene_name,"|", aa, proteinSites)) %>%
  dplyr::select(ID, contains("CT")) %>%
  group_by(ID) %>%
  summarise(across(everything(), ~ if (all(is.na(.))) NA_real_ else na.omit(.)[1]), .groups = "drop")

colnames(proc_24804581) <- str_remove(colnames(proc_24804581), "/CT")
## Prepare data ---------------------------
# PMID: 30257219, EGF treatment, select early time points (2,4,8 minutes) as strongest EGFR activation is expected
proc_30257219 <- read_excel(raw_30257219, sheet = 2) %>%
  dplyr::select(gene.name, modified.sites, log2.fold.change.2min, log2.fold.change.4min, log2.fold.change.8min) %>%
  filter(!str_detect(modified.sites, ",")) %>%
  filter(!gene.name == "Putative uncharacterized protein") %>%
  mutate(ID = paste0(gene.name, "_", modified.sites, "|", gene.name, "|", modified.sites)) %>%
  dplyr::select(ID, contains("log2"))
proc_30257219 <- proc_30257219[!duplicated(proc_30257219$ID),]

colnames(proc_30257219) <- str_remove(colnames(proc_30257219), "log2.fold.change.")
colnames(proc_30257219)[2:ncol(proc_30257219)] <- paste("EGF", colnames(proc_30257219)[2:ncol(proc_30257219)], sep = "_")

## Merge datasets ---------------------------
full_df <- full_join(proc_24362263, proc_17389395, by = "ID") %>%
  full_join(proc_17016520, by = "ID") %>%
  full_join(proc_25147952, by = "ID") %>%
  full_join(proc_24804581, by = "ID") %>%
  full_join(proc_30257219, by = "ID")

full_df <- full_df[!duplicated(full_df$ID),]

## Meta data ---------------------------
meta_df <- read_csv(meta_overview, col_types = cols())
meta_proc <- meta_df %>%
  dplyr::rename("id" = Name, "sign" = Direction, "target" = Kinase) %>%
  mutate(target = str_remove(target, " Down")) %>%
  mutate(target = str_remove(target, " Up")) %>%
  dplyr::select(id, sign, target) %>%
  separate_rows(target, sep = ", ")

meta_proc <- meta_proc %>%
  filter(id %in% colnames(full_df))


# Save data
write_csv(full_df, bench_output)
write_csv(meta_proc, meta_prior_output)

