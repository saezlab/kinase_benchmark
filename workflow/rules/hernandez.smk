# ------------------------------------ DATA PROCESSING ------------------------------------
rule format_input:
    input:
        phospho = "data/hernandez/annotations.xlsx",
        meta = "data/hernandez/benchmark_data.xlsx"
    output:
        mat = "results/hernandez/processed_data/benchmark_data.csv",
        meta_out = "results/hernandez/processed_data/benchmark_metadata.csv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/hernandez/prepareData.R"

# ------------------------------ INPUT PREPARATION ------------------------------
# ------------------------------ Prior knowledge preparation ------------------------------
rule prepare_omnipath:
    input:
        file_dataset = "results/hernandez/processed_data/benchmark_data.csv"
    output:
        tsv = "results/hernandez/prior/omnipath.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hernandez/01.1_prepare_omnipath.R"

rule prepare_phosphositeplus:
    input:
        ppsp = "data/prior/phosphositeplus",
        file_dataset = "results/hernandez/processed_data/benchmark_data.csv"
    output:
        tsv = "results/hernandez/prior/phosphositeplus.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hernandez/01.2_prepare_phosphositeplus.R"

rule prepare_ptmsigdb:
    input:
        ptmsig_file = "data/prior/ptm.sig.db.all.uniprot.human.v1.9.0.gmt",
        file_dataset = "results/hernandez/processed_data/benchmark_data.csv"
    output:
        tsv = "results/hernandez/prior/ptmsigdb.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hernandez/01.3_prepare_ptmsigdb.R"

rule prepare_ikipdb:
    input:
        ptmsig_file = "data/prior/iKiP-DB-Table.tsv",
        file_dataset = "results/hernandez/processed_data/benchmark_data.csv"
    output:
        tsv = "results/hernandez/prior/iKiPdb.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hernandez/01.4_prepare_ikipdb.R"

rule prepare_GPS:
    input:
        GPS_file = "data/prior/mmc4.xlsx",
        file_dataset = "results/hernandez/processed_data/benchmark_data.csv"
    output:
        tsv = "results/hernandez/prior/GPS.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hernandez/01.5_prepare_GPS.R"

rule prepare_NetworKIN:
    input:
        networkin_file = "data/prior/networkin_human_predictions_3.1.tsv",
        file_dataset = "results/hernandez/processed_data/benchmark_data.csv"
    output:
        tsv = "results/hernandez/prior/networkin.tsv"
    params:
    	score = 5
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hernandez/01.6_prepare_NetworKIN.R"

rule prepare_jhonson:
    input:
        ppsp = "data/prior/simplified_jhonson_with_psite.csv",
        file_dataset = "results/hernandez/processed_data/benchmark_data.csv"
    output:
        tsv = "results/hernandez/prior/jhonson.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hernandez/01.9_prepare_jhonson.R"

rule merge_GPS_PPSP:
    input:
        gps = "results/hernandez/prior/GPS.tsv",
        ppsp = "results/hernandez/prior/phosphositeplus.tsv"
    output:
        tsv = "results/hernandez/prior/GPSppsp.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hernandez/01.7_merge_GPS_PPSP.R"

rule merge_known_predicted:
    input:
        known_file = "results/hernandez/prior/{known_targets}.tsv",
        predicted_file = "results/hernandez/prior/{predicted_targets}.tsv"
    output:
        tsv = "results/hernandez/prior/{known_targets}_{predicted_targets}.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hernandez/01.8_merge_known_predicted.R"


# ------------------------------- PTM-SEA input preparation -------------------------------
rule ptmsea_datasets:
    input:
        file_dataset = "results/hernandez/processed_data/benchmark_data.csv"
    output:
        gct = "results/hernandez/datasets/benchmark.gct"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hernandez/02.1_prepare_datasets_ptm-sea.R"

rule ptmsea_prior:
    input:
        file_PKN = "results/hernandez/prior/{PKN}.tsv"
    output:
        gmt = "results/hernandez/prior/ptm-sea/{PKN}.gmt"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hernandez/02.2_prepare_prior_ptm-sea.R"


# ---------------------------------- ACTIVITY ESTIMATION ----------------------------------
rule activity_estimation:
    input:
        file_dataset = "results/hernandez/processed_data/benchmark_data.csv",
        file_PKN = "results/hernandez/prior/{PKN}.tsv",
        scripts = expand("workflow/scripts/methods/run_{method}.R", method = ["INKA", "KARP", "lm_rokai", "zscore", "erics_methods"]),
        script_support = "workflow/scripts/methods/support_functions.R"
    output:
        rds = "results/hernandez/activity_scores/{PKN}.rds"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_kinase_activity_estimation/hernandez/03_activity_estimation_hernandez.R"


rule activity_estimation_ptmsea:
    input:
        file_dataset = "results/hernandez/datasets/benchmark.gct",
        file_PKN = "results/hernandez/prior/ptm-sea/{PKN}.gmt"
    output:
        rds = "results/hernandez/activity_scores_ptmsea/log/{PKN}.log",
        gct = "results/hernandez/activity_scores_ptmsea/{PKN}-scores.gct"
    params:
        output_folder = "results/hernandez/activity_scores_ptmsea"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_kinase_activity_estimation/hernandez/03.1_activity_estimation_ptmsea.R"

rule combine_scores:
    input:
        file_ptmsea = "results/hernandez/activity_scores_ptmsea/{PKN}-scores.gct",
        file_scores = "results/hernandez/activity_scores/{PKN}.rds"
    output:
        rds = "results/hernandez/final_scores/{PKN}.rds"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_kinase_activity_estimation/hernandez/04_combine_scores.R"


# -------------------------------------- BENCHMARK ---------------------------------------
rule prepare_benchmark:
    input:
        rds = "results/hernandez/final_scores/{PKN}.rds",
        meta = "results/hernandez/processed_data/benchmark_metadata.csv",
    output:
        output = "results/hernandez/benchmark_files/{hernandez_methods}-{PKN}.csv",
        meta_out = "results/hernandez/benchmark_files/obs_{hernandez_methods}-{PKN}.csv"

    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/01_prepare_bench_input.R"

rule run_benchmark:
    input:
        scores = "results/hernandez/benchmark_files/{hernandez_methods}-{PKN}.csv",
        meta = "results/hernandez/benchmark_files/obs_{hernandez_methods}-{PKN}.csv"
    output:
        output = "results/hernandez/benchmark_res/{PKN}/bench_{hernandez_methods}-{PKN}.csv"
    conda:
        "../envs/benchmark.yml"
    script:
        "../scripts/04_benchmark/hernandez/02_decouple_bench.py"
