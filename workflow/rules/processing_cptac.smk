# This rule is to prepare the datasets for kinase activity estimation of cptac
# ------------------------------ DATA FORMATTING ------------------------------
rule fromat_data:
    input:
        phospho = "data/datasets/cptac/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds"
    output:
        out = "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/cptac/00_format_data.R"


# ------------------------------ INPUT PREPARATION ------------------------------
rule map_priors:
    input:
        ppsp = "results/00_prior/{prior}.tsv",
        file_dataset = expand("results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv", dataset = config["cptac"]["datasets"], normalisation = config["cptac"]["normalisation"])
    output:
        tsv = "results/01_processed_data/cptac/mapped_priors/{prior}.tsv"
    wildcard_constraints:
        prior = '(?!shuffled)[a-zA-Z0-9]+'
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/cptac/01_prior_mapping.R"

rule map_merged_priors:
    input:
        ppsp = "results/00_prior/merged/{known}_{predicted}.tsv",
        file_dataset = expand("results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv", dataset = config["cptac"]["datasets"], normalisation = config["cptac"]["normalisation"])
    output:
        tsv = "results/01_processed_data/cptac/mapped_priors/{known}_{predicted}.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/cptac/01_prior_mapping.R"

rule map_kinase_ids:
    input:
        prior_files = expand("results/01_processed_data/cptac/mapped_priors/{PKN}.tsv", PKN = config["cptac"]["cptac_PKNs"])
    output:
        tsv = "resources/kinase_mapping.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/cptac/02_convert_kinase_ids.R"
        
rule prepare_cptacjohnson:
    input:
        ppsp = "data/kinase_libraries/johnson_library/kinase_benchmarking_pancancer_site_table_st_percentiles.tsv",
        tyr = "data/kinase_libraries/johnson_library/kinase_benchmarking_pancancer_site_table_tyr_percentiles.tsv",
        known = "results/01_processed_data/cptac/mapped_priors/phosphositeplus.tsv"
    params:
        perc = lambda w: w.perc
    output:
        tsv = "results/01_processed_data/cptac/mapped_priors/phosphositeplus_cptacjohnson{perc}.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/04_prepare_johnson_cptac.R"

# ------------------------------- PTM-SEA input preparation -------------------------------
rule ptmsea_datasets:
    input:
        file_dataset = "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv"
    output:
        gct = temp("results/01_processed_data/cptac/datasets/{dataset}_{normalisation}.gct")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/02.1_prepare_datasets_ptm-sea.R"

rule ptmsea_prior:
    input:
        file_PKN = "results/01_processed_data/cptac/mapped_priors/{PKN}.tsv"
    output:
        gmt = temp("results/01_processed_data/cptac/mapped_priors/ptm-sea/{PKN}.gmt")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/02.2_prepare_prior_ptm-sea.R"
