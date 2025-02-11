
# -------------------------------------- BENCHMARK ---------------------------------------
rule merge_scores:
    input:
        rds =  "results/02_activity_scores/hijazi/scores/{PKN}.rds",
        hernandez = "results/02_activity_scores/hernandez/scores/{PKN}.rds"
    output:
        output = "results/02_activity_scores/merged/scores/{PKN}.rds"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_benchmark/00_merge_scores.R"
        
rule merge_scores2:
    input:
        rds =  "results/02_activity_scores/hijazi/scores/{PKN}.rds",
        hernandez = "results/02_activity_scores/hernandez/scores/{PKN}.rds",
        tyr = "results/02_activity_scores/tyrosine/scores/{PKN}.rds"
    output:
        output = "results/02_activity_scores/merged2/scores/{PKN}.rds"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_benchmark/00_merge_scores_new.R"

rule copy_scores:
    input:
        rds =  "results/02_activity_scores/hijazi/scores/{PKN}.rds",
    output:
        output = temp("results/02_activity_scores/hijaziDiscoverX/scores/{PKN}.rds")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_benchmark/00_copy_scores.R"

rule merge_meta:
    input:
        meta = "results/01_processed_data/hijazi/data/benchmark_metadata.csv",
        metah = "results/01_processed_data/hernandez/data/benchmark_metadata.csv"
    output:
        meta_out = "results/01_processed_data/merged/data/benchmark_metadata.csv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_benchmark/00_merge_meta.R"

rule merge_meta2:
    input:
        meta = "results/01_processed_data/hijazi/data/benchmark_metadata.csv",
        metah = "results/01_processed_data/hernandez/data/benchmark_metadata.csv",
        metatyr = "results/01_processed_data/tyrosine/data/benchmark_metadata.csv"
    output:
        meta_out = "results/01_processed_data/merged2/data/benchmark_metadata.csv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_benchmark/00_merge_meta_new.R"

rule prepare_benchmark:
    input:
        rds = "results/02_activity_scores/{dataset}/scores/{PKN}.rds",
        meta = "results/01_processed_data/{dataset}/data/benchmark_metadata.csv"
    output:
        output = temp("results/03_benchmark/{dataset}/01_input_bench/{hernandez_methods}-{PKN}.csv"),
        meta_out = temp("results/03_benchmark/{dataset}/01_input_bench/obs_{hernandez_methods}-{PKN}.csv")
    params:
        rm_exp = "F",
        scale = "T"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_benchmark/01_prepare_input.R"

rule run_benchmark:
    input:
        scores = "results/03_benchmark/{dataset}/01_input_bench/{hernandez_methods}-{PKN}.csv",
        meta = "results/03_benchmark/{dataset}/01_input_bench/obs_{hernandez_methods}-{PKN}.csv"
    output:
        output = "results/03_benchmark/{dataset}/02_benchmark_res/{PKN}/bench_{hernandez_methods}-{PKN}.csv"
    conda:
        "../envs/benchmark.yml"
    script:
        "../scripts/03_benchmark/02_run_bench.py"

rule mean_rank:
    input:
        scores = "results/03_benchmark/{dataset}/01_input_bench/{hernandez_methods}-{PKN}.csv",
        meta = "results/03_benchmark/{dataset}/01_input_bench/obs_{hernandez_methods}-{PKN}.csv",
        target = "results/02_activity_scores/{dataset}/scores/{PKN}.rds",
        kinclass = "resources/kinase_class.csv"
    output:
        output = "results/03_benchmark/{dataset}/02_mean_rank/{PKN}/{hernandez_methods}-{PKN}.csv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_benchmark/02_mean_rank.R"

rule compare_rank:
    input:
        bench = expand("results/03_benchmark/{{dataset}}/02_mean_rank/{PKN}/{hernandez_methods}-{PKN}.csv", hernandez_methods = config["perturbation"]["methods"], PKN = config["perturbation"]["PKNs"])
    output:
        ov =  "results/03_benchmark/{dataset}/03_benchmark_comp/overview.csv",
        rank = "results/03_benchmark/{dataset}/03_benchmark_comp/combined_ranks.csv",
        kin = "results/03_benchmark/{dataset}/03_benchmark_comp/kinase_ranks.csv",
        exp = "results/03_benchmark/{dataset}/03_benchmark_comp/experiment_ranks.csv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_benchmark/03_compare_rank.R"


# -------------------------------------- SUBSET ---------------------------------------
def get_methods(wildcards):
    if (wildcards['subset'] == "johnson"):
        methods = config["perturbation"]["methods_johnson"]
    elif (wildcards['subset'] == "subset"):
        methods = config["perturbation"]["methods"]
    else:
        methods = config["perturbation"]["methods_subset"]
    return expand("results/03_benchmark/{dataset}/01_input_bench/{methods}-{PKNs}.csv", methods = methods, PKNs = config["perturbation"][wildcards['subset']], dataset = wildcards['dataset'])

rule generate_subset:
    input:
        scores = get_methods
    params:
    	  subset = lambda w: w.subset
    output:
        output = temp("results/03_benchmark/{dataset}/01_input_bench_subset/{subset}/filter_subset.csv")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_benchmark/001_generate_subset.R"

rule prepare_subset:
    input:
        scores = "results/03_benchmark/{dataset}/01_input_bench/{methods}-{PKNs}.csv",
        filter = "results/03_benchmark/{dataset}/01_input_bench_subset/{subset}/filter_subset.csv",
        meta = "results/03_benchmark/{dataset}/01_input_bench/obs_{methods}-{PKNs}.csv"
    params:
    	  subset = lambda w: w.subset
    output:
        output = temp("results/03_benchmark/{dataset}/01_input_bench_subset/{subset}/{methods}-{PKNs}.csv"),
        meta_out = temp("results/03_benchmark/{dataset}/01_input_bench_subset/{subset}/obs_{methods}-{PKNs}.csv")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_benchmark/001_input_subset.R"

rule run_benchmark_subset:
    input:
        scores = "results/03_benchmark/{dataset}/01_input_bench_subset/{subset}/{hernandez_methods}-{PKN}.csv",
        meta = "results/03_benchmark/{dataset}/01_input_bench_subset/{subset}/obs_{hernandez_methods}-{PKN}.csv"
    params:
    	  subset = lambda w: w.subset
    output:
        output = "results/03_benchmark/{dataset}/02_benchmark_res_subset/{subset}/{PKN}/bench_{hernandez_methods}-{PKN}.csv"
    conda:
        "../envs/benchmark.yml"
    script:
        "../scripts/03_benchmark/02_run_bench.py"

rule mean_rank_subset:
    input:
        scores = "results/03_benchmark/{dataset}/01_input_bench_subset/{subset}/{hernandez_methods}-{PKN}.csv",
        meta = "results/03_benchmark/{dataset}/01_input_bench_subset/{subset}/obs_{hernandez_methods}-{PKN}.csv",
        target = "results/02_activity_scores/{dataset}/scores/{PKN}.rds",
        kinclass = "resources/kinase_class.csv"
    output:
        output = "results/03_benchmark/{dataset}/02_mean_rank_subset/{subset}/{PKN}/{hernandez_methods}-{PKN}.csv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/03_benchmark/02_mean_rank.R"
