# ------------------------------------ FIGURE 1 ------------------------------------
rule overview_bench:
    input:
        bench_files = expand("results/03_benchmark/merged/02_benchmark_res/phosphositeplus/bench_{methods}-phosphositeplus.csv", methods = config["perturbation"]["methods"]),
        rank = expand("results/03_benchmark/merged/02_mean_rank/phosphositeplus/{methods}-phosphositeplus.csv", methods = config["perturbation"]["methods"])
    params:
        k_phit_c = [5, 10, 20]
    output:
        auroc = "results/manuscript_figures/figure_1/auroc_phosphositeplus.pdf",
        rank = "results/manuscript_figures/figure_1/scaledrank_phosphositeplus.pdf",
        phit =  "results/manuscript_figures/figure_1/phit_phosphositeplus.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/01_benchmark_overview.R"

rule supp_figure1_overview:
    input:
        meta = "results/01_processed_data/merged/data/benchmark_metadata.csv",
        bench = expand("results/03_benchmark/{dataset}/02_benchmark_res/phosphositeplus/bench_{methods}-phosphositeplus.csv", dataset = ["hijazi", "hijaziDiscoverX"], methods = config["perturbation"]["methods"]),
        rank = expand("results/03_benchmark/{dataset}/02_mean_rank/phosphositeplus/{methods}-phosphositeplus.csv", dataset = ["hijazi", "hijaziDiscoverX"], methods = config["perturbation"]["methods"])
    output:
        ove = "results/manuscript_figures/figure_1/supp/overview_kin.pdf",
        out = "results/manuscript_figures/figure_1/supp/hijazi.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/01.1_supp_overview.R"

rule supp_figure1_rank:
    input:
        rank = expand("results/03_benchmark/merged/02_mean_rank/{PKN}/{methods}-{PKN}.csv", PKN = ["phosphositeplus"], methods = ["zscore"])
    output:
        plot = "results/manuscript_figures/figure_1/supp/comparison_rank_subset.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/01.1_supp_rank.R"
        

# ------------------------------------ FIGURE 2 ------------------------------------
rule overview_tumor:
    input:
        bench = "data/tumor_benchmark/activity_scores/roc_data.rds",
        act = "data/tumor_benchmark/activity_scores/psp_roc_actsiteBM.rds"
    output:
        plot = "results/manuscript_figures/figure_2/auroc_tumor_phosphositeplus.pdf",
        plotact = "results/manuscript_figures/figure_2/auroc_act_phosphositeplus.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/01_cptac_overview.R"
      

# ------------------------------------ FIGURE 3 ------------------------------------
rule prior_coverage:
    input:
        prior_files = expand("results/00_prior/{PKN}.tsv", PKN = config["figures"]["coverage"])
    output:
        kin="results/manuscript_figures/figure_3/coverage_kin.pdf",
        kin_heat="results/manuscript_figures/figure_3/kinase_overview.pdf",
        edges="results/manuscript_figures/figure_3/coverage_edge.pdf",
        edges_heat="results/manuscript_figures/figure_3/edge_overview.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/03_prior_coverage.R"
        
rule performance_priors:
    input:
        bench = expand("results/03_benchmark/merged/02_benchmark_res/{PKN}/bench_zscore-{PKN}.csv", PKN = config["figures"]["evaluation"]),
        rank = expand("results/03_benchmark/merged/02_mean_rank/{PKN}/zscore-{PKN}.csv", PKN = config["figures"]["evaluation"]),
        act = expand("data/tumor_benchmark/activity_scores/{PKN_cptac}_roc_actsiteBM.rds", PKN_cptac = congif["figures"]["evaluation_cptac"])
        tumor = "data/tumor_benchmark/activity_scores/roc_data.rds"
    output:
        plt="results/manuscript_figures/figure_3/zscore.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/03_prior_evaluation.R"
        


        
rule benchmark_figure:
    input:
        bench_files = expand("results/03_benchmark/merged/02_benchmark_res/{PKN}/bench_{methods}-{PKN}.csv", PKN = config["figures"]["PKN_figure2"], methods = config["perturbation"]["methods"]),
        meta = "results/01_processed_data/merged/data/benchmark_metadata.csv",
        rank = expand("results/03_benchmark/merged/02_mean_rank/{PKN}/{methods}-{PKN}.csv", PKN = config["figures"]["PKN_figure2"], methods = config["perturbation"]["methods"]),
        cit = "resources/protein_citations.csv",
        ppsp = "results/01_processed_data/cptac/mapped_priors/phosphositeplus.tsv"
    output:
        auroc = "results/manuscript_figures/figure_3/auroc_res.pdf",
        meta_over = "results/manuscript_figures/figure_3/overview_kin.pdf",
        rankPlt = "results/manuscript_figures/figure_3/mean_rank.pdf",
        rankKin = "results/manuscript_figures/figure_3/kinase_GSknown.csv",
        heat = "results/manuscript_figures/figure_3/median_auroc.pdf",
        medRank = "results/manuscript_figures/figure_3/median_rank.pdf",
        cor = "results/manuscript_figures/figure_3/study_bias.pdf",
        target = "results/manuscript_figures/figure_3/target_bias.pdf",
        statPrior = "results/manuscript_figures/supp_files/prior_comparison.csv",
        statMethod = "results/manuscript_figures/supp_files/method_comparison.csv"
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
        methods = config["perturbation"]["methods"]
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/02_figure_comparison_activity.R"

# ------------------------------------ SUPP. FIGURE ------------------------------------
rule supp_figure:
    input:
        bench_files = expand("results/03_benchmark/merged/02_benchmark_res_subset/predicted/{PKN}/bench_{methods}-{PKN}.csv", PKN = config["perturbation"]["predicted"], methods = config["perturbation"]["methods"]),
        #predicted = expand("results/03_benchmark/merged/02_benchmark_res_subset/johnson/{PKN}/bench_{methods}-{PKN}.csv", PKN = config["perturbation"]["johnson"], methods = config["perturbation"]["methods_johnson"])
    output:
        merged = "results/manuscript_figures/supp_figures/combined.pdf",
        pred = "results/manuscript_figures/supp_figures/predicted_combined.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/04_figure_combined_priors.R"
