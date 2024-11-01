rule compare_performance:
    input:
        bench = expand("results/03_benchmark/{{dataset}}/02_benchmark_res/{PKN}/bench_{hernandez_methods}-{PKN}.csv", hernandez_methods = config["perturbation"]["methods"], PKN = config["perturbation"]["PKNs"]),
        rank = expand("results/03_benchmark/{{dataset}}/02_mean_rank/{PKN}/{hernandez_methods}-{PKN}.csv", hernandez_methods = config["perturbation"]["methods"], PKN = config["perturbation"]["PKNs"])
    params:
        k_phit = 10
    output:
        plot = "results/04_exploration/{dataset}/benchmark/plots/performance.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/04_exploration/compare_performance.R"

rule compare_performance_subset:
    input:
        bench = expand("results/03_benchmark/{{dataset}}/02_benchmark_res_subset/{{subset}}/{PKN}/bench_{hernandez_methods}-{PKN}.csv", hernandez_methods = config["perturbation"]["methods"], PKN = config["perturbation"]["subset"]),
        rank = expand("results/03_benchmark/{{dataset}}/02_mean_rank_subset/{{subset}}/{PKN}/{hernandez_methods}-{PKN}.csv", hernandez_methods = config["perturbation"]["methods"], PKN = config["perturbation"]["subset"])
    params:
        k_phit = 10
    output:
        plot = "results/04_exploration/{dataset}/benchmark/plots/performance_{subset}.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/04_exploration/compare_performance.R"

rule compare_performance_tumort:
    input:
        bench = expand("results/03_benchmark/merged/02_benchmark_res/{PKN}/bench_{hernandez_methods}-{PKN}.csv", hernandez_methods = config["perturbation"]["methods"], PKN = ["phosphositeplus"]),
        tumor = "data/misc/known_roc.Rds"
    output:
        plot = "results/04_exploration/all/performance.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/04_exploration/compare_tumor_perturbation.R"

rule compare_performance_tumor_subset:
    input:
        bench = expand("results/03_benchmark/merged/02_benchmark_res_subset/{{subset}}/{PKN}/bench_{hernandez_methods}-{PKN}.csv", hernandez_methods = config["perturbation"]["methods"], PKN = ["phosphositeplus"]),
        tumor = "data/misc/known_roc_filt.Rds"
    output:
        plot = "results/04_exploration/all/performance_{subset}.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/04_exploration/compare_tumor_perturbation.R"


rule test_performance:
    input:
        scores = "results/03_benchmark/{dataset}/02_mean_rank/{PKN}/{method}-{PKN}.csv"
    output:
        out = "results/04_exploration/{dataset}/mean_rank/plots/{method}-{PKN}_targets.pdf"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_exploration/compare_class_regulon.R"

rule calculate_correlation:
    input:
        act_files = expand("results/02_activity_scores/{{dataset}}/scores/{PKN}.rds", PKN = config["figures"]["coverage"]),
    params:
        corr_type = lambda w: w.corr_type
    output:
        corrMethods = "results/04_exploration/{dataset}/correlation/correlation_methods_{corr_type}.rds",
        corrPriors = "results/04_exploration/{dataset}/correlation/correlation_priors_{corr_type}.rds",
        plotMethods = "results/04_exploration/{dataset}/correlation/plots/correlation_methods_{corr_type}.pdf",
        plotPriors = "results/04_exploration/{dataset}/correlation/plots/correlation_priors_{corr_type}.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/04_exploration/calculate_correlation.R"

rule calculate_jaccard:
    input:
        act_files = expand("results/02_activity_scores/{{dataset}}/scores/{PKN}.rds", PKN = config["figures"]["coverage"]), 
    params:
        jaccard_i = lambda w: w.jaccard_i,
        jacc_type = lambda w: w.jacc_type
    output:
        jaccMethods = "results/04_exploration/{dataset}/jaccard/jaccard_methods_{jacc_type}_{jaccard_i}.rds",
        jaccPriors = "results/04_exploration/{dataset}/jaccard/jaccard_priors_{jacc_type}_{jaccard_i}.rds",
        plotMethods = "results/04_exploration/{dataset}/jaccard/plots/jaccard_methods_{jacc_type}_{jaccard_i}.pdf",
        plotPriors = "results/04_exploration/{dataset}/jaccard/plots/jaccard_priors_{jacc_type}_{jaccard_i}.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/04_exploration/calculate_jaccard.R"

rule explore_rank:
    input:
        rank = "results/03_benchmark/{dataset}/02_mean_rank/{PKN}/{hernandez_methods}-{PKN}.csv",
    output:
        plot = "results/04_exploration/{dataset}/rank/comparison_rank_{hernandez_methods}_{PKN}.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/04_exploration/explore_rank.R"
