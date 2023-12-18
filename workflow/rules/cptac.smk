# ------------------------------ DATA FORMATTING ------------------------------
rule prepare_omnipath:
    input:
        phospho = "data/CPTAC_original/{dataset}_original_medcent_30plus.tsv"
    output:
        out = "data/CPTAC_phospho/final/{dataset}_norm2prot_original_lm_log2_medCentRatio.rds"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/cptac/format_unnormalized_data.R"

# ------------------------------ INPUT PREPARATION ------------------------------
rule map_priors:
    input:
        ppsp = "results/prior/{prior}.tsv",
        file_dataset = expand("data/CPTAC_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds", dataset = config["cptac"]["datasets"], normalisation = config["cptac"]["normalisation"])
    output:
        tsv = "results/cptac/prior/{prior}.tsv"
    wildcard_constraints:
        prior = '[a-zA-Z]+'
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/cptac/01_prior_mapping.R"

rule map_merged_priors:
    input:
        ppsp = "results/prior/merged/{known}_{predicted}.tsv",
        file_dataset = expand("data/CPTAC_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds", dataset = config["cptac"]["datasets"], normalisation = config["cptac"]["normalisation"])
    output:
        tsv = "results/cptac/prior/{known}_{predicted}.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/cptac/01_prior_mapping.R"

rule map_kinase_ids:
    input:
        prior_files = expand("results/cptac/prior/{PKN}.tsv", PKN = config["cptac"]["cptac_PKNs"])
    output:
        tsv = "resources/kinase_mapping.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/cptac/02_convert_kinase_ids.R"

# ------------------------------- PTM-SEA input preparation -------------------------------
rule ptmsea_datasets:
    input:
        file_dataset = "data/CPTAC_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds"
    output:
        gct = "results/cptac/datasets/{dataset}_{normalisation}.gct"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/cptac/02.1_prepare_datasets_ptm-sea.R"

rule ptmsea_prior:
    input:
        file_PKN = "results/cptac/prior/{PKN}.tsv"
    output:
        gmt = "results/cptac/prior/ptm-sea/{PKN}.gmt"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/cptac/02.2_prepare_prior_ptm-sea.R"


# ---------------------------------- ACTIVITY ESTIMATION ----------------------------------
# ------------------------------ CPTAC ------------------------------
rule activity_estimation:
    input:
        file_dataset ="data/CPTAC_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds",
        file_PKN = "results/cptac/prior/{PKN}.tsv",
        scripts = expand("workflow/scripts/methods/run_{method}.R", method = ["INKA", "KARP", "lm_rokai", "zscore", "erics_methods"]),
        script_support = "workflow/scripts/methods/support_functions.R"
    output:
        rds = "results/cptac/activity_scores/{PKN}/{normalisation}/{normalisation}_{dataset}-{PKN}.rds"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_kinase_activity_estimation/cptac/03_activity_estimation.R"

rule activity_estimation_ptmsea:
    input:
        file_dataset = "results/cptac/datasets/{dataset}_{normalisation}.gct",
        file_PKN = "results/cptac/prior/ptm-sea/{PKN}.gmt"
    output:
        rds = "results/cptac/activity_scores_ptmsea/log/{normalisation}_{dataset}-{PKN}.log",
        gct = "results/cptac/activity_scores_ptmsea/{normalisation}_{dataset}-{PKN}-scores.gct"
    params:
        output_folder = "results/cptac/activity_scores_ptmsea"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_kinase_activity_estimation/cptac/03.1_activity_estimation_ptmsea.R"

rule combine_scores:
    input:
        file_ptmsea = "results/cptac/activity_scores_ptmsea/{normalisation}_{dataset}-{PKN}-scores.gct",
        file_scores = "results/cptac/activity_scores/{PKN}/{normalisation}/{normalisation}_{dataset}-{PKN}.rds"
    output:
        rds = "results/cptac/final_scores/{PKN}/{normalisation}/{normalisation}_{dataset}-{PKN}.rds"
    params:
        rm_methods = ["corr_wmean", "corr_wsum", "norm_wsum", "wmean", "INKA_kinase_centric", "INKA_substrate_centric", "INKA"]
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_kinase_activity_estimation/cptac/04_combine_scores.R"


# -------------------------------------- Comparison ---------------------------------------
rule prior_comparison:
    input:
        prior_files = expand("results/prior/{PKN}.tsv", PKN = config["cptac"]["cptac_PKNs"])
    output:
        kin = "results/comparison/plots/coverage_kinases.pdf",
        edges = "results/comparison/plots/coverage_edges.pdf",
        pps = "results/comparison/plots/coverage_pps.pdf"
    params:
        plot_width = "13.7",
        plot_height = "6"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04.1_compare_prior.R"

rule activity_comparison:
    input:
        act_files = expand("results/final_scores/{PKN}/{dataset}_{PKN}.rds", PKN = config["cptac"]["cptac_PKNs"], dataset = config["cptac"]["datasets"])
    output:
        plotSpearman = "results/comparison/plots/spearman_heatmap.pdf",
        plotPearson = "results/comparison/plots/pearson_heatmap.pdf",
        plotJaccardUp = "results/comparison/plots/jaccard_up_heatmap.pdf",
        plotJaccardDown = "results/comparison/plots/jaccard_down_heatmap.pdf",
        dist_csv = "results/comparison/mean_distance.csv"
    params:
        jaccard_i = "10"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04.2_compare_activities.R"
