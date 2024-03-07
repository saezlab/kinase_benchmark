# ------------------------------------ FIGURE 1 ------------------------------------
rule overview_prior:
    input:
        prior_files = expand("results/00_prior/{PKN}.tsv", PKN = config["figures"]["PKN_figure1"])
    output:
        kin = "results/manuscript_figures/figure_1/coverage_merged.pdf",
        edges = "results/manuscript_figures/figure_1/edge_overview.pdf",
        pps =  "results/manuscript_figures/figure_1/jaccard.pdf",
        kin_heat = "results/manuscript_figures/figure_1/kinase_overview.pdf",
        kin_type = "results/manuscript_figures/figure_1/kinase_type.pdf",
        reg = "results/manuscript_figures/figure_1/regulon_size.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/01_figure_overview_priors.R"

# ------------------------------------ FIGURE 2 ------------------------------------
rule act_comparison:
    input:
        bench = "results/01_processed_data/hernandez/data/benchmark_data.csv",
        benchMeta = "results/01_processed_data/hernandez/data/benchmark_metadata.csv",
        hijazi = "results/01_processed_data/hijazi/data/benchmark_data.csv",
        hijMeta = "results/01_processed_data/hijazi/data/benchmark_metadata.csv",
        prior_files = expand("results/01_processed_data/hernandez/mapped_priors/{PKN}.tsv", PKN = config["figures"]["PKN_figure1"]),
        hijPrior = expand("results/01_processed_data/hijazi/mapped_priors/{PKN}.tsv", PKN = config["figures"]["PKN_figure1"]),
        act =  expand("results/02_activity_scores/merged/final_scores/{PKN}.rds", PKN = config["figures"]["PKN_figure1"])
    output:
        overview = "results/manuscript_figures/figure_2/overview_experiment.pdf",
        corrMeth = "results/manuscript_figures/figure_2/corrplot_methods.pdf",
        corrPrior = "results/manuscript_figures/figure_2/corrplot_priors.pdf"
    params:
        methods = config["hernandez"]["hernandez_methods"]
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/02_figure_comparison_activity.R"

# ------------------------------------ FIGURE 3 ------------------------------------
rule benchmark_figure:
    input:
        bench_files = expand("results/03_benchmark/merged/02_benchmark_res/{PKN}/bench_{methods}-{PKN}.csv", PKN = config["figures"]["PKN_figure2"], methods = config["hernandez"]["hernandez_methods"]),
        meta = "results/01_processed_data/merged/data/benchmark_metadata.csv",
        rank = expand("results/03_benchmark/merged/02_mean_rank/{PKN}/{methods}-{PKN}.csv", PKN = config["figures"]["PKN_figure2"], methods = config["hernandez"]["hernandez_methods"]),
        cit = "resources/protein_citations.csv"
    output:
        auroc = "results/manuscript_figures/figure_3/auroc_res.pdf",
        meta_over = "results/manuscript_figures/figure_3/overview_kin.pdf",
        rankPlt = "results/manuscript_figures/figure_3/mean_rank.pdf",
        rankKin = "results/manuscript_figures/figure_3/kinase_GSknown.csv",
        heat = "results/manuscript_figures/figure_3/median_auroc.pdf",
        medRank = "results/manuscript_figures/figure_3/median_rank.pdf",
        cor = "results/manuscript_figures/figure_3/study_bias.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/03_figure_benchmark.R"

rule citation_info:
    input:
        gene = "data/gene2pubmed",
        info = "data/gene_info"
    output:
        out = "resources/protein_citations.csv"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/01_data_processing/protein_citations.R"
