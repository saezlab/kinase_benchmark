# ------------------------------------ FIGURE 1 ------------------------------------
rule overview_prior:
    input:
        prior_files = expand("results/prior/{PKN}.tsv", PKN = config["figures"]["PKN_figure1"])
    output:
        kin = "results/manuscript_figures/figure_1/coverage_merged.pdf",
        edges = "results/manuscript_figures/figure_1/edge_overview.pdf",
        pps =  "results/manuscript_figures/figure_1/jaccard.pdf",
        kin_heat = "results/manuscript_figures/figure_1/kinase_overview.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/figures_manuscript/01_figure_overview_priors.R"

# ------------------------------------ FIGURE 2 ------------------------------------
rule act_comparison:
    input:
        bench = "results/hernandez/processed_data/benchmark_data.csv",
        benchMeta = "results/hernandez/processed_data/benchmark_metadata.csv",
        hijazi = "results/hijazi/01_processed_data/benchmark_data.csv",
        hijMeta = "results/hijazi/01_processed_data/benchmark_metadataPrior.csv",
        prior_files = expand("results/hernandez/prior/{PKN}.tsv", PKN = config["figures"]["PKN_figure1"]),
        hijPrior = expand("results/hijazi/02_prior/{PKN}.tsv", PKN = config["figures"]["PKN_figure1"]),
        act =  expand("results/hernandez/final_scores/scaled/{PKN}.rds", PKN = config["figures"]["PKN_figure1"]),
        hijAct = expand("results/hijazi/04_final_scores/scaled/{PKN}.rds", PKN = config["figures"]["PKN_figure1"])
    output:
        overview = "results/manuscript_figures/figure_2/overview_experiment.pdf",
        corrMeth = "results/manuscript_figures/figure_2/corrplot_methods.pdf",
        corrPrior = "results/manuscript_figures/figure_2/corrplot_priors.pdf"
    params:
        methods = config["hernandez"]["hernandez_methods"]
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/figures_manuscript/02_figure_comparison_activity.R"

# ------------------------------------ FIGURE 3 ------------------------------------
rule benchmark_figure:
    input:
        bench_files = expand("results/hijazi/06_benchmark_res/{PKN}/merged/bench_{methods}-{PKN}.csv", PKN = config["figures"]["PKN_figure2"], methods = config["hernandez"]["hernandez_methods"]),
        meta = "results/hernandez/processed_data/benchmark_metadata.csv",
        rank = "results/hijazi/06_mean_rank/full_rank_merged.csv",
        hijazi = "results/hijazi/01_processed_data/benchmark_metadataPrior.csv",
        kin = "results/hijazi/06_mean_rank/performance_per_kin_merged.csv",
        prior = "results/hernandez/overview_priors/coverage.csv"
    output:
        auroc = "results/manuscript_figures/figure_3/auroc_res.pdf",
        meta_over = "results/manuscript_figures/figure_3/overview_kin.pdf",
        rankPlt = "results/manuscript_figures/figure_3/mean_rank.pdf",
        rankKin = "results/manuscript_figures/figure_3/kinase_GSknown.csv",
        heat = "results/manuscript_figures/figure_3/median_auroc.pdf",
        medRank = "results/manuscript_figures/figure_3/median_rank.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/figures_manuscript/03_figure_benchmark.R"
