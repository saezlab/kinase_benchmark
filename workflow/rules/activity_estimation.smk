# This rule is to run the activity estimation

# ------------------------------ CPTAC ------------------------------
rule cptac_activity_estimation:
    input:
        file_dataset ="data/CPTAC_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds",
        file_PKN = "results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts = expand("workflow/scripts/methods/run_{method}.R", method = ["INKA", "KARP", "lm_rokai", "zscore", "erics_methods"]),
        script_support = "workflow/scripts/methods/support_functions.R"
    output:
        rds = "results/02_activity_scores/cptac/scores/{normalisation}/{normalisation}_{dataset}-{PKN}.rds"
    params:
        rm_auto = "T",
        minsize = "5"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_activity_estimation.R"

rule cptac_activity_estimation_ptmsea:
    input:
        file_dataset = "results/01_processed_data/cptac/datasets/{dataset}_{normalisation}.gct",
        file_PKN = "results/01_processed_data/cptac/mapped_priors/ptm-sea/{PKN}.gmt"
    output:
        rds = "results/02_activity_scores/cptac/ptmsea/log/{normalisation}_{dataset}-{PKN}.log",
        gct = "results/02_activity_scores/cptac/ptmsea/{normalisation}_{dataset}-{PKN}-scores.gct"
    params:
        output_folder = "results/02_activity_scores/cptac/ptmsea"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01.1_activity_estimation_ptmsea_cptac.R"

rule cptac_combine_scores:
    input:
        file_ptmsea = "results/02_activity_scores/cptac/ptmsea/{normalisation}_{dataset}-{PKN}-scores.gct",
        file_scores = "results/02_activity_scores/cptac/scores/{normalisation}/{normalisation}_{dataset}-{PKN}.rds"
    output:
        rds = "results/02_activity_scores/cptac/final_scores/{normalisation}/{normalisation}_{dataset}-{PKN}.rds"
    params:
        rm_methods = ["corr_wmean", "corr_wsum", "norm_wsum", "wmean", "INKA_kinase_centric", "INKA_substrate_centric", "INKA"]
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/02_combine_scores.R"

# ---------------------------------- hernandez ----------------------------------
rule hernandez_activity_estimation:
    input:
        file_dataset = "results/01_processed_data/hernandez/data/benchmark_data.csv",
        file_PKN = "results/01_processed_data/hernandez/mapped_priors/{PKN}.tsv",
        scripts = expand("workflow/scripts/methods/run_{method}.R", method = ["INKA", "KARP", "lm_rokai", "zscore", "erics_methods"]),
        script_support = "workflow/scripts/methods/support_functions.R"
    output:
        rds = "results/02_activity_scores/hernandez/scores/{PKN}.rds"
    params:
        rm_auto = "T",
        minsize = "3"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_activity_estimation.R"

rule hernandez_activity_estimation_ptmsea:
    input:
        file_dataset = "results/01_processed_data/hernandez/datasets/benchmark.gct",
        file_PKN = "results/01_processed_data/hernandez/mapped_priors/ptm-sea/{PKN}.gmt"
    output:
        rds = "results/02_activity_scores/hernandez/ptmsea/log/{PKN}.log",
        gct = "results/02_activity_scores/hernandez/ptmsea/{PKN}-scores.gct"
    params:
        output_folder = "results/02_activity_scores/hernandez/ptmsea",
        minsize = "3"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01.1_activity_estimation_ptmsea.R"

rule hernandez_combine_scores:
    input:
        file_ptmsea = "results/02_activity_scores/hernandez/ptmsea/{PKN}-scores.gct",
        file_scores = "results/02_activity_scores/hernandez/scores/{PKN}.rds"
    output:
        rds = "results/02_activity_scores/hernandez/final_scores/{PKN}.rds"
    params:
        rm_methods = ["corr_wmean", "corr_wsum", "norm_wsum", "wmean", "INKA_kinase_centric", "INKA_substrate_centric", "INKA"]
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/02_combine_scores.R"


# ---------------------------------- hijazi ----------------------------------
rule hijazi_activity_estimation:
    input:
        file_dataset = "results/01_processed_data/hijazi/data/benchmark_data.csv",
        file_PKN = "results/01_processed_data/hijazi/mapped_priors/{PKN}.tsv",
        scripts = expand("workflow/scripts/methods/run_{method}.R", method = ["INKA", "KARP", "lm_rokai", "zscore", "erics_methods"]),
        script_support = "workflow/scripts/methods/support_functions.R"
    output:
        rds = "results/02_activity_scores/hijazi/scores/{PKN}.rds"
    params:
        rm_auto = "T",
        minsize = "3"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_activity_estimation.R"

rule hijazi_activity_estimation_ptmsea:
    input:
        file_dataset = "results/01_processed_data/hijazi/datasets/benchmark.gct",
        file_PKN = "results/01_processed_data/hijazi/mapped_priors/ptm-sea/{PKN}.gmt"
    output:
        rds = "results/02_activity_scores/hijazi/ptmsea/log/{PKN}.log",
        gct = "results/02_activity_scores/hijazi/ptmsea/{PKN}-scores.gct"
    params:
        output_folder = "results/02_activity_scores/hijazi/ptmsea",
        minsize = "3"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01.1_activity_estimation_ptmsea.R"

rule hijazi_combine_scores:
    input:
        file_ptmsea = "results/02_activity_scores/hijazi/ptmsea/{PKN}-scores.gct",
        file_scores = "results/02_activity_scores/hijazi/scores/{PKN}.rds"
    output:
        rds = "results/02_activity_scores/hijazi/final_scores/{PKN}.rds"
    params:
        rm_methods = ["corr_wmean", "corr_wsum", "norm_wsum", "wmean", "INKA_kinase_centric", "INKA_substrate_centric", "INKA"]
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/02_combine_scores.R"
