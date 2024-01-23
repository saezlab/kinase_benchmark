# ------------------------------------ DATA PROCESSING ------------------------------------
rule format_input:
    input:
        HL60 = "data/hijazi/HL60_fc.tsv",
        MCF7 = "data/hijazi/MCF7_fc.tsv",
        HL60p = "data/hijazi/HL60_pval.tsv",
        MCF7p = "data/hijazi/MCF7_pval.tsv",
        meta = "data/hijazi/inhibitor_selectivity.tsv",
        targets = "data/hijazi/targets_hijazi.tsv"
    output:
        mat = "results/hijazi/01_processed_data/benchmark_data.csv",
        meta_out = "results/hijazi/01_processed_data/benchmark_metadataPrior.csv",
        meta_discover = "results/hijazi/01_processed_data/benchmark_metadataDiscoverX.csv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/01_data_processing/hijazi/prepareData.R"

# ------------------------------ INPUT PREPARATION ------------------------------
# ------------------------------ Prior knowledge preparation ------------------------------
rule map_priors:
    input:
        ppsp = "results/prior/{prior}.tsv",
        file_dataset = "results/hijazi/01_processed_data/benchmark_data.csv"
    output:
        tsv = "results/hijazi/02_prior/{prior}.tsv"
    wildcard_constraints:
        prior = '[a-zA-Z]+'
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hijazi/01_prior_mapping.R"

rule map_merged_priors:
    input:
        ppsp = "results/prior/merged/{known}_{predicted}.tsv",
        file_dataset = "results/hijazi/01_processed_data/benchmark_data.csv"
    output:
        tsv = "results/hijazi/02_prior/{known}_{predicted}.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hijazi/01_prior_mapping.R"


# ------------------------------- PTM-SEA input preparation -------------------------------
rule ptmsea_datasets:
    input:
        file_dataset = "results/hijazi/01_processed_data/benchmark_data.csv"
    output:
        gct = "results/hijazi/02_datasets/benchmark.gct"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hijazi/02.1_prepare_datasets_ptm-sea.R"

rule ptmsea_prior:
    input:
        file_PKN = "results/hijazi/02_prior/{PKN}.tsv"
    output:
        gmt = "results/hijazi/02_prior/ptm-sea/{PKN}.gmt"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/hijazi/02.2_prepare_prior_ptm-sea.R"

# ---------------------------------- ACTIVITY ESTIMATION ----------------------------------
rule activity_estimation:
    input:
        file_dataset = "results/hijazi/01_processed_data/benchmark_data.csv",
        file_PKN = "results/hijazi/02_prior/{PKN}.tsv",
        scripts = expand("workflow/scripts/methods/run_{method}.R", method = ["INKA", "KARP", "lm_rokai", "zscore", "erics_methods"]),
        script_support = "workflow/scripts/methods/support_functions.R"
    output:
        rds = "results/hijazi/03_activity_scores/{PKN}.rds"
    params:
        rm_auto = "T",
        minsize = "3"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_kinase_activity_estimation/hijazi/03_activity_estimation_hijazi.R"


rule activity_estimation_ptmsea:
    input:
        file_dataset = "results/hijazi/02_datasets/benchmark.gct",
        file_PKN = "results/hijazi/02_prior/ptm-sea/{PKN}.gmt"
    output:
        rds = "results/hijazi/03_activity_scores_ptmsea/log/{PKN}.log",
        gct = "results/hijazi/03_activity_scores_ptmsea/{PKN}-scores.gct"
    params:
        output_folder = "results/hijazi/03_activity_scores_ptmsea",
        minsize = "3"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_kinase_activity_estimation/hijazi/03.1_activity_estimation_ptmsea.R"

rule combine_scores:
    input:
        file_ptmsea = "results/hijazi/03_activity_scores_ptmsea/{PKN}-scores.gct",
        file_scores = "results/hijazi/03_activity_scores/{PKN}.rds"
    output:
        rds = "results/hijazi/04_final_scores/raw/{PKN}.rds"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_kinase_activity_estimation/hijazi/04_combine_scores.R"


# -------------------------------------- BENCHMARK ---------------------------------------
rule scale_scores:
    input:
        rds = "results/hijazi/04_final_scores/raw/{PKN}.rds"
    output:
        output = "results/hijazi/04_final_scores/scaled/{PKN}.rds"
    params:
        scale = "sd"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/00_scale_scores.R"


rule prepare_benchmark:
    input:
        rds = "results/hijazi/04_final_scores/scaled/{PKN}.rds",
        meta = "results/hijazi/01_processed_data/benchmark_metadataPrior.csv"
    output:
        output = "results/hijazi/05_benchmark_files/scaled/{hernandez_methods}-{PKN}.csv",
        meta_out = "results/hijazi/05_benchmark_files/scaled/obs_{hernandez_methods}-{PKN}.csv"
    params:
        rm_exp = "F"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/01_prepare_bench_input.R"

rule merge_benchmark:
    input:
        rds = "results/hijazi/04_final_scores/scaled/{PKN}.rds",
        hernandez = "results/hernandez/final_scores/scaled/GSknown.rds",
        meta = "results/hijazi/01_processed_data/benchmark_metadataPrior.csv",
        metah = "results/hernandez/processed_data/benchmark_metadata.csv"
    output:
        output = "results/hijazi/05_benchmark_files/merged/{hernandez_methods}-{PKN}.csv",
        meta_out = "results/hijazi/05_benchmark_files/merged/obs_{hernandez_methods}-{PKN}.csv"
    params:
        rm_exp = "F"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hijazi/01_prepare_bench_input_merged.R"

rule run_benchmark:
    input:
        scores = "results/hijazi/05_benchmark_files/{overlap}/{hernandez_methods}-{PKN}.csv",
        meta = "results/hijazi/05_benchmark_files/{overlap}/obs_{hernandez_methods}-{PKN}.csv"
    output:
        output = "results/hijazi/06_benchmark_res/{PKN}/{overlap}/bench_{hernandez_methods}-{PKN}.csv"
    conda:
        "../envs/benchmark.yml"
    script:
        "../scripts/04_benchmark/hernandez/02_decouple_bench.py"

rule compare_performance:
    input:
        bench = expand("results/hijazi/06_benchmark_res/{PKN}/{{overlap}}/bench_{hernandez_methods}-{PKN}.csv", hernandez_methods = config["hernandez"]["hernandez_methods"], PKN = config["hernandez"]["hernandez_PKNs"])
    output:
        auroc = "results/hijazi/06_benchmark_res/plots/AUROC_{overlap}.pdf",
        auprc = "results/hijazi/06_benchmark_res/plots/AUPRC_{overlap}.pdf"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/03_compare_performance.R"

# -------------------------------------- RANK ---------------------------------------
rule run_mean_rank:
    input:
        rds = expand("results/hijazi/05_benchmark_files/{{overlap}}/{hernandez_methods}-{PKN}.csv", hernandez_methods = config["hernandez"]["hernandez_methods"], PKN = config["hernandez"]["hernandez_PKNs"]),
        meta =  "results/hernandez/processed_data/benchmark_metadata.csv",
        hijazi = "results/hijazi/01_processed_data/benchmark_metadataPrior.csv",
        overview = "results/hernandez/overview_priors/coverage.csv"
    output:
        output = "results/hijazi/06_mean_rank/mean_rank_{overlap}.csv",
        cov_kin = "results/hijazi/06_mean_rank/overview/covered_kinases_{overlap}.csv",
        per_exp = "results/hijazi/06_mean_rank/performance_per_exp_{overlap}.csv",
        per_kin = "results/hijazi/06_mean_rank/performance_per_kin_{overlap}.csv",
        full = "results/hijazi/06_mean_rank/full_rank_{overlap}.csv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hijazi/02_mean_rank_merged.R"
