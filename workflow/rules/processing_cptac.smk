# This rule is to prepare the datasets for kinase activity estimation of cptac
# ------------------------------ DATA FORMATTING ------------------------------
rule fromat_data:
    input:
        phospho = "data/datasets/cptac_original/{dataset}_original_medcent_30plus.tsv"
    output:
        out = "data/datasets/cptac_phospho/final/{dataset}_norm2prot_original_lm_log2_medCentRatio.rds"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/cptac/00_format_unnormalized_data.R"

# ------------------------------ INPUT PREPARATION ------------------------------
rule map_priors:
    input:
        ppsp = "results/00_prior/{prior}.tsv",
        file_dataset = expand("data/datasets/cptac_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds", dataset = config["cptac"]["datasets"], normalisation = config["cptac"]["normalisation"])
    output:
        tsv = "results/01_processed_data/cptac/mapped_priors/{prior}.tsv"
    wildcard_constraints:
        prior = '[a-zA-Z0-9]+'
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/cptac/01_prior_mapping.R"

rule map_merged_priors:
    input:
        ppsp = "esults/00_prior/merged/{known}_{predicted}.tsv",
        file_dataset = expand("data/datasets/cptac_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds", dataset = config["cptac"]["datasets"], normalisation = config["cptac"]["normalisation"])
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

# ------------------------------- PTM-SEA input preparation -------------------------------
rule ptmsea_datasets:
    input:
        file_dataset = "data/datasets/cptac_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds"
    output:
        gct = "results/01_processed_data/cptac/datasets/{dataset}_{normalisation}.gct"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/cptac/02.1_prepare_datasets_ptm-sea.R"

rule ptmsea_prior:
    input:
        file_PKN = "results/01_processed_data/cptac/mapped_priors/{PKN}.tsv"
    output:
        gmt = "results/01_processed_data/cptac/mapped_priors/ptm-sea/{PKN}.gmt"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/cptac/02.2_prepare_prior_ptm-sea.R"
