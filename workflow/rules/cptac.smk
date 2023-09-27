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
rule prepare_omnipath:
    input:
        file_dataset = expand("data/CPTAC_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds", dataset = config["activity_estimation"]["datasets"], normalisation = config["activity_estimation"]["normalisation"])
    output:
        tsv = "results/cptac/prior/omnipath.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/cptac/01.1_prepare_omnipath.R"

rule prepare_phosphositeplus:
    input:
        ppsp = "data/prior/phosphositeplus",
        file_dataset = expand("data/CPTAC_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds", dataset = config["activity_estimation"]["datasets"], normalisation = config["activity_estimation"]["normalisation"])
    output:
        tsv = "results/cptac/prior/phosphositeplus.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/cptac/01.2_prepare_phosphositeplus.R"

rule prepare_ptmsigdb:
    input:
        ptmsig_file = "data/prior/ptm.sig.db.all.uniprot.human.v1.9.0.gmt",
        file_dataset = expand("data/CPTAC_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds", dataset = config["activity_estimation"]["datasets"], normalisation = config["activity_estimation"]["normalisation"])
    output:
        tsv = "results/cptac/prior/ptmsigdb.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/cptac/01.3_prepare_ptmsigdb.R"

rule prepare_ikipdb:
    input:
        ptmsig_file = "data/prior/iKiP-DB-Table.tsv",
        file_dataset = expand("data/CPTAC_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds", dataset = config["activity_estimation"]["datasets"], normalisation = config["activity_estimation"]["normalisation"])
    output:
        tsv = "results/cptac/prior/iKiPdb.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/cptac/01.4_prepare_ikipdb.R"

rule prepare_GPS:
    input:
        GPS_file = "data/prior/mmc4.xlsx",
        file_dataset = expand("data/CPTAC_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds", dataset = config["activity_estimation"]["datasets"], normalisation = config["activity_estimation"]["normalisation"])
    output:
        tsv = "results/cptac/prior/GPS.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/cptac/01.5_prepare_GPS.R"

rule prepare_NetworKIN:
    input:
        networkin_file = "data/prior/networkin_human_predictions_3.1.tsv",
        file_dataset = expand("data/CPTAC_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds", dataset = config["activity_estimation"]["datasets"], normalisation = config["activity_estimation"]["normalisation"])
    output:
        tsv = "results/cptac/prior/networkin.tsv"
    params:
    	score = 5
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/cptac/01.6_prepare_NetworKIN.R"

rule prepare_jhonson:
    input:
        ppsp = "data/prior/simplified_jhonson_with_psite.csv",
        file_dataset = expand("data/CPTAC_phospho/final/{dataset}_norm2prot_{normalisation}_lm_log2_medCentRatio.rds", dataset = config["activity_estimation"]["datasets"], normalisation = config["activity_estimation"]["normalisation"])
    output:
        tsv = "results/cptac/prior/jhonson.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/cptac/01.9_prepare_jhonson.R"

rule merge_GPS_PPSP:
    input:
        gps = "results/cptac/prior/GPS.tsv",
        ppsp = "results/cptac/prior/phosphositeplus.tsv"
    output:
        tsv = "results/cptac/prior/GPSppsp.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/cptac/01.7_merge_GPS_PPSP.R"

rule merge_known_predicted:
    input:
        known_file = "results/cptac/prior/{known_targets}.tsv",
        predicted_file = "results/cptac/prior/{predicted_targets}.tsv"
    output:
        tsv = "results/cptac/prior/{known_targets}_{predicted_targets}.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/cptac/01.8_merge_known_predicted.R"

rule map_kinase_ids:
    input:
        prior_files = expand("results/cptac/prior/{PKN}.tsv", PKN = config["activity_estimation"]["PKNs"])
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
        prior_files = expand("results/prior/{PKN}.tsv", PKN = config["activity_estimation"]["PKNs"])
    output:
        kin = "results/comparison/plots/coverage_kinases.pdf",
        edges = "results/comparison/plots/coverage_edges.pdf",
        pps = "results/comparison/plots/coverage_pps.pdf"
    params:
        plot_width = "13",
        plot_height = "6"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04.1_compare_prior.R"

rule activity_comparison:
    input:
        act_files = expand("results/final_scores/{PKN}/{dataset}_{PKN}.rds", PKN = config["activity_estimation"]["PKNs"], dataset = config["activity_estimation"]["datasets"])
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
