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

rule prior_overview:
    input:
        prior_files = expand("results/hernandez/prior/{PKN}.tsv", PKN = config["hernandez"]["hernandez_PKNs"])
    output:
        csv = "results/hernandez/overview_priors/coverage.csv",
        kin = "results/hernandez/overview_priors/coverage_kinases.pdf",
        edges = "results/hernandez/overview_priors/coverage_edges.pdf",
        pps = "results/hernandez/overview_priors/coverage_pps.pdf"
    params:
        height = "6",
        width = "13"
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
rule scale_scores:
    input:
        rds = "results/hernandez/final_scores/{PKN}.rds"
    output:
        output = "results/hernandez/final_scores/scaled/{PKN}.rds"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/00_scale_scores.R"

rule prepare_benchmark:
    input:
        rds = "results/hernandez/final_scores/scaled/{PKN}.rds",
        meta = "results/hernandez/processed_data/benchmark_metadata.csv"
    output:
        output = "results/hernandez/benchmark_files/{hernandez_methods}-{PKN}.csv",
        meta_out = "results/hernandez/benchmark_files/obs_{hernandez_methods}-{PKN}.csv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/01_prepare_bench_input.R"

rule run_mean_rank:
    input:
        rds = expand("results/hernandez/final_scores/scaled/{PKN}.rds", PKN = config["hernandez"]["hernandez_PKNs"]),
        meta =  "results/hernandez/processed_data/benchmark_metadata.csv",
        overview = "results/hernandez/overview_priors/coverage.csv"
    output:
        output = "results/hernandez/benchmark_mean_rank/mean_rank.csv",
        per_exp = "results/hernandez/benchmark_mean_rank/performance_per_exp.csv",
        per_kin = "results/hernandez/benchmark_mean_rank/performance_per_kin.csv",
        pdf = "results/hernandez/benchmark_mean_rank/mean_rank.pdf",
        boxplot = "results/hernandez/benchmark_mean_rank/bp_rank.pdf"
    params:
        mth = config["hernandez"]["hernandez_PKNs"]
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/02_decouple_bench.py"

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

rule compare_performance:
    input:
        bench = expand("results/hernandez/benchmark_res/{PKN}/bench_{hernandez_methods}-{PKN}.csv", hernandez_methods = config["hernandez"]["hernandez_methods"], PKN = config["hernandez"]["hernandez_PKNs"])
    output:
        auroc = "results/hernandez/benchmark_res/plots/AUROC.pdf",
        auprc = "results/hernandez/benchmark_res/plots/AUPRC.pdf"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/03_compare_performance.R"

rule overview_bench:
    input:
        bench = expand("results/hernandez/benchmark_files/obs_{hernandez_methods}-{PKN}.csv", hernandez_methods = config["hernandez"]["hernandez_methods"], PKN = config["hernandez"]["hernandez_PKNs"])
    output:
        ov = "results/hernandez/benchmark_res/overview/overview_bench.res",
        ov_prior = "results/hernandez/benchmark_res/overview/overview_bench_prior.res"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/04_overview_bench.R"
