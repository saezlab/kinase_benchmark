# This rule is to prepare the datasets for kinase activity estimation of tyrosine

# ------------------------------------ DATA PROCESSING ------------------------------------
rule format_input:
    input:
        id24362263 = "data/datasets/tyrosine/NIHMS551505-supplement-2.xlsx",
        id17389395 = "data/datasets/tyrosine/pnas_0608638104_08638Table3.xls",
        id17016520 = "data/datasets/tyrosine/msb4100094-s3.xls",
        id25147952 = "data/datasets/tyrosine/TableS1.xlsx",
        id24804581 = "data/datasets/tyrosine/cb500116c_si_004.xls",
        id30257219 = "data/datasets/tyrosine/processed.xlsx",
        meta = "data/datasets/tyrosine/Overview.csv"
    output:
        mat = "results/01_processed_data/tyrosine/data/benchmark_data.csv",
        meta_out = "results/01_processed_data/tyrosine/data/benchmark_metadata.csv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/tyrosine/00_prepareData.R"

# ------------------------------ INPUT PREPARATION ------------------------------
# ------------------------------ Prior knowledge preparation ------------------------------
rule map_priors:
    input:
        ppsp = "results/00_prior/{prior}.tsv",
        file_dataset = "results/01_processed_data/tyrosine/data/benchmark_data.csv"
    output:
        tsv = temp("results/01_processed_data/tyrosine/mapped_priors/{prior}.tsv")
    wildcard_constraints:
        prior = '(?!shuffled)[a-zA-Z0-9]+'
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/tyrosine/01_prior_mapping.R"

rule map_merged_priors:
    input:
        ppsp = "results/00_prior/merged/{known}_{predicted}.tsv",
        file_dataset = "results/01_processed_data/tyrosine/data/benchmark_data.csv"
    output:
        tsv = temp("results/01_processed_data/tyrosine/mapped_priors/{known}_{predicted}.tsv")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/tyrosine/01_prior_mapping.R"


# ------------------------------- PTM-SEA input preparation -------------------------------
rule ptmsea_datasets:
    input:
        file_dataset = "results/01_processed_data/tyrosine/data/benchmark_data.csv"
    output:
        gct = "results/01_processed_data/tyrosine/datasets/benchmark.gct"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/02.1_prepare_datasets_ptm-sea.R"

rule ptmsea_prior:
    input:
        file_PKN = "results/01_processed_data/tyrosine/mapped_priors/{PKN}.tsv"
    output:
        gmt = "results/01_processed_data/tyrosine/mapped_priors/ptm-sea/{PKN}.gmt"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/02.2_prepare_prior_ptm-sea.R"
