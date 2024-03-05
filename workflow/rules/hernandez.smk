
# -------------------------------------- BENCHMARK ---------------------------------------
rule scale_scores:
    input:
        rds = "results/hernandez/final_scores/raw/{PKN}.rds"
    output:
        output = "results/hernandez/final_scores/scaled/{PKN}.rds"
    params:
        scale = "max"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/00_scale_scores.R"

rule get_subset:
    input:
        rds_raw = expand("results/hernandez/final_scores/raw/{PKN}.rds", PKN = config["hernandez"]["hernandez_PKNs_subset"]),
        rds_scaled = expand("results/hernandez/final_scores/scaled/{PKN}.rds", PKN = config["hernandez"]["hernandez_PKNs_subset"])
    wildcard_constraints:
        PKN = '[.*phosphositeplus.*]'
    output:
        output = "results/hernandez/final_scores/subset/{PKN}.rds",
        output_scaled = "results/hernandez/final_scores/scaled_subset/{PKN}.rds"
    params:
        methods = config["hernandez"]["hernandez_methods"]
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/001_get_shared_subset.R"

rule prepare_benchmark:
    input:
        rds = "results/hernandez/final_scores/{overlap}/{PKN}.rds",
        meta = "results/hernandez/processed_data/benchmark_metadata.csv"
    output:
        output = "results/hernandez/benchmark_files/{overlap}/{hernandez_methods}-{PKN}.csv",
        meta_out = "results/hernandez/benchmark_files/{overlap}/obs_{hernandez_methods}-{PKN}.csv"
    params:
        rm_exp = "F"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/01_prepare_bench_input.R"

rule run_mean_rank:
    input:
        rds = expand("results/hernandez/benchmark_files/{{overlap}}/{hernandez_methods}-{PKN}.csv", hernandez_methods = config["hernandez"]["hernandez_methods"], PKN = config["hernandez"]["hernandez_PKNs"]),
        meta =  "results/hernandez/processed_data/benchmark_metadata.csv",
        overview = "results/hernandez/overview_priors/coverage.csv"
    output:
        output = "results/hernandez/benchmark_mean_rank/mean_rank_{overlap}.csv",
        cov_kin = "results/hernandez/benchmark_res/overview/covered_kinases_{overlap}.csv",
        per_exp = "results/hernandez/benchmark_mean_rank/performance_per_exp_{overlap}.csv",
        per_kin = "results/hernandez/benchmark_mean_rank/performance_per_kin_{overlap}.csv",
        pdf = "results/hernandez/benchmark_mean_rank/mean_rank_{overlap}.pdf",
        boxplot = "results/hernandez/benchmark_mean_rank/bp_rank_{overlap}.pdf"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/02.2_get_mean_rank.R"

rule run_benchmark:
    input:
        scores = "results/hernandez/benchmark_files/{overlap}/{hernandez_methods}-{PKN}.csv",
        meta = "results/hernandez/benchmark_files/{overlap}/obs_{hernandez_methods}-{PKN}.csv"
    output:
        output = "results/hernandez/benchmark_res/{PKN}/{overlap}/bench_{hernandez_methods}-{PKN}.csv"
    conda:
        "../envs/benchmark.yml"
    script:
        "../scripts/04_benchmark/hernandez/02_decouple_bench.py"

rule compare_performance:
    input:
        bench = expand("results/hernandez/benchmark_res/{PKN}/{{overlap}}/bench_{hernandez_methods}-{PKN}.csv", hernandez_methods = config["hernandez"]["hernandez_methods"], PKN = config["hernandez"]["hernandez_PKNs"])
    output:
        auroc = "results/hernandez/benchmark_res/plots/AUROC_{overlap}.pdf",
        auprc = "results/hernandez/benchmark_res/plots/AUPRC_{overlap}.pdf"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/03_compare_performance.R"

rule overview_bench:
    input:
        bench = expand("results/hernandez/benchmark_files/{{overlap}}/obs_{hernandez_methods}-{PKN}.csv", hernandez_methods = config["hernandez"]["hernandez_methods"], PKN = config["hernandez"]["hernandez_PKNs"])
    output:
        ov = "results/hernandez/benchmark_res/overview/overview_bench_{overlap}.res",
        ov_prior = "results/hernandez/benchmark_res/overview/overview_bench_prior_{overlap}.res"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_benchmark/hernandez/04_overview_bench.R"



## Hijazi
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
        rds = expand("results/hijazi/05_benchmark_files/{{overlap}}/{hernandez_methods}-{PKN}.csv", hernandez_methods = config["hernandez"]["hernandez_methods"], PKN = config["figures"]["PKN_figure2"]),
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

