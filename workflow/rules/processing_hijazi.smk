# This rule is to prepare the datasets for kinase activity estimation of hijazi

# ------------------------------------ DATA PROCESSING ------------------------------------
rule format_input:
    input:
        HL60 = "data/datasets/hijazi/HL60_fc.tsv",
        MCF7 = "data/datasets/hijazi/MCF7_fc.tsv",
        HL60p = "data/datasets/hijazi/HL60_pval.tsv",
        MCF7p = "data/datasets/hijazi/MCF7_pval.tsv",
        meta = "data/datasets/hijazi/inhibitor_selectivity.tsv",
        targets = "data/datasets/hijazi/targets_hijazi.tsv"
    output:
        mat = "results/01_processed_data/hijazi/data/benchmark_data.csv",
        meta_out = "results/01_processed_data/hijazi/data/benchmark_metadata.csv",
        meta_discover = "results/01_processed_data/hijazi/data/benchmark_metadataDiscoverX.csv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/hijazi/00_prepareData.R"

# ------------------------------ INPUT PREPARATION ------------------------------
# ------------------------------ Prior knowledge preparation ------------------------------
rule map_priors:
    input:
        ppsp = "results/00_prior/{prior}.tsv",
        file_dataset = "results/01_processed_data/hijazi/data/benchmark_data.csv"
    output:
        tsv = temp("results/01_processed_data/hijazi/mapped_priors/{prior}.tsv")
    wildcard_constraints:
        prior = '(?!shuffled)[a-zA-Z0-9]+'
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/hijazi/01_prior_mapping.R"

rule map_merged_priors:
    input:
        ppsp = "results/00_prior/merged/{known}_{predicted}.tsv",
        file_dataset = "results/01_processed_data/hijazi/data/benchmark_data.csv"
    output:
        tsv = temp("results/01_processed_data/hijazi/mapped_priors/{known}_{predicted}.tsv")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/hijazi/01_prior_mapping.R"


# ------------------------------- PTM-SEA input preparation -------------------------------
rule ptmsea_datasets:
    input:
        file_dataset = "results/01_processed_data/hijazi/data/benchmark_data.csv"
    output:
        gct = "results/01_processed_data/hijazi/datasets/benchmark.gct"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/02.1_prepare_datasets_ptm-sea.R"

rule ptmsea_prior:
    input:
        file_PKN = "results/01_processed_data/hijazi/mapped_priors/{PKN}.tsv"
    output:
        gmt = "results/01_processed_data/hijazi/mapped_priors/ptm-sea/{PKN}.gmt"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/02.2_prepare_prior_ptm-sea.R"
