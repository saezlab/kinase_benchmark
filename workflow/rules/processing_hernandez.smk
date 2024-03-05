# This rule is to prepare the datasets for kinase activity estimation of hernandez
# ------------------------------------ DATA PROCESSING ------------------------------------
rule format_input:
    input:
        phospho = "data/datasets/hernandez/annotations.xlsx",
        meta = "data/datasets/hernandez/benchmark_data.xlsx"
    output:
        mat = "results/01_processed_data/hernandez/data/benchmark_data.csv",
        meta_out = "results/01_processed_data/hernandez/data/benchmark_metadata.csv"
    params:
        msk_exp = ["705_225", "1298_272", "1288_272", "1291_272", "387_117", "1289_272", "1290_272", "1308_272", "699_225"]
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/hernandez/00_prepareData.R"

# ------------------------------ INPUT PREPARATION ------------------------------
# ------------------------------ Prior knowledge preparation ------------------------------
rule map_priors:
    input:
        ppsp = "results/00_prior/{prior}.tsv",
        file_dataset = "results/01_processed_data/hernandez/data/benchmark_data.csv"
    output:
        tsv = "results/01_processed_data/hernandez/mapped_priors/{prior}.tsv"
    wildcard_constraints:
        prior = '[a-zA-Z0-9]+'
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/hernandez/01_prior_mapping.R"

rule map_merged_priors:
    input:
        ppsp = "results/00_prior/merged/{known}_{predicted}.tsv",
        file_dataset = "results/01_processed_data/hernandez/data/benchmark_data.csv"
    output:
        tsv = "results/01_processed_data/hernandez/mapped_priors/{known}_{predicted}.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/hernandez/01_prior_mapping.R"


# ------------------------------- PTM-SEA input preparation -------------------------------
rule ptmsea_datasets:
    input:
        file_dataset = "results/01_processed_data/hernandez/data/benchmark_data.csv"
    output:
        gct = "results/01_processed_data/hernandez/datasets/benchmark.gct"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/hernandez/02.1_prepare_datasets_ptm-sea.R"

rule ptmsea_prior:
    input:
        file_PKN = "results/01_processed_data/hernandez/mapped_priors/{PKN}.tsv"
    output:
        gmt = "results/01_processed_data/hernandez/mapped_priors/ptm-sea/{PKN}.gmt"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/hernandez/02.2_prepare_prior_ptm-sea.R"
