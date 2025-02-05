# ------------------------------------ FIGURE 1 ------------------------------------
rule overview_bench:
    input:
        bench_files = expand("results/03_benchmark/merged2/02_benchmark_res/phosphositeplus/bench_{methods}-phosphositeplus.csv", methods = config["perturbation"]["methods"]),
        rank = expand("results/03_benchmark/merged2/02_mean_rank/phosphositeplus/{methods}-phosphositeplus.csv", methods = config["perturbation"]["methods"])
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
        meta = "results/01_processed_data/merged2/data/benchmark_metadata.csv",
        bench = expand("results/03_benchmark/{dataset}/02_benchmark_res/phosphositeplus/bench_{methods}-phosphositeplus.csv", dataset = ["hijazi", "hijaziDiscoverX"], methods = config["perturbation"]["methods"]),
        rank = expand("results/03_benchmark/{dataset}/02_mean_rank/phosphositeplus/{methods}-phosphositeplus.csv", dataset = ["hijazi", "hijaziDiscoverX"], methods = config["perturbation"]["methods"]),
        ppspHer = "results/01_processed_data/hernandez/mapped_priors/phosphositeplus.tsv",
        ppspHij = "results/01_processed_data/hijazi/mapped_priors/phosphositeplus.tsv",
        ppspTyr = "results/01_processed_data/tyrosine/mapped_priors/phosphositeplus.tsv",
        her = "results/01_processed_data/hernandez/data/benchmark_data.csv",
        hij = "results/01_processed_data/hijazi/data/benchmark_data.csv",
        tyr = "results/01_processed_data/tyrosine/data/benchmark_data.csv"
    output:
        ove = "results/manuscript_figures/figure_1/supp/overview_kin.pdf",
        oveFil = "results/manuscript_figures/figure_1/supp/overview_kin_filtered.pdf",
        out = "results/manuscript_figures/figure_1/supp/hijazi.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/01.1_supp_overview.R"

rule supp_figure1_rank:
    input:
        rank = expand("results/03_benchmark/merged2/02_mean_rank/{PKN}/{methods}-{PKN}.csv", PKN = ["phosphositeplus"], methods = ["zscore"])
    output:
        plot = "results/manuscript_figures/figure_1/supp/comparison_rank_subset.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/01.1_supp_rank.R"

rule supp_figure1_tyr:
    input:
        rank = expand("results/03_benchmark/merged2/02_mean_rank/{PKN}/{methods}-{PKN}.csv", PKN = ["phosphositeplus"], methods = config["perturbation"]["methods"]),
        kinclass = "resources/kinase_class.csv"
    output:
        out = "results/manuscript_figures/figure_1/supp/rank_class.pdf",
        phit = "results/manuscript_figures/figure_1/supp/phit_class.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/01.1_supp_tyrosine.R"

rule supp_table_hijazi:
    input:
        meta = "results/01_processed_data/hijaziDiscoverX/data/benchmark_metadata.csv"
    output:
        out = "results/manuscript_figures/supp_files/hijazi_meta.csv"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/01.1_supp_hijazi.R"
        

# ------------------------------------ FIGURE 2 ------------------------------------
rule overview_tumor:
    input:
        act="data/results_cptac/overall_performance/actsiteBM/all_kins/actsiteBM_5perThr_psp_roc_table.rds",
        bench="data/results_cptac/overall_performance/protBM/all_kins/protBM_5perThr_psp_roc_table.rds",
        norm_act=expand("data/results_cptac/normalisation/actsiteBM/actsiteBM_5perThr_{norm}_roc_table.rds", norm = config["figures"]["normalisation"]),
        prot_act=expand("data/results_cptac/normalisation/proteinBM/protBM_5perThr_{norm}_roc_table.rds", norm = config["figures"]["normalisation_prot"])
    output:
        plot = "results/manuscript_figures/figure_2/auroc_tumor_phosphositeplus.pdf",
        plotact = "results/manuscript_figures/figure_2/auroc_act_phosphositeplus.pdf",
        norm_plot = "results/manuscript_figures/figure_2/norm_tumor_phosphositeplus.pdf",
        normact_plot = "results/manuscript_figures/figure_2/norm_act_phosphositeplus.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/02_cptac_overview.R"
      

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
        bench = expand("results/03_benchmark/merged2/02_benchmark_res/{PKN}/bench_zscore-{PKN}.csv", PKN = config["figures"]["evaluation"]),
        rank = expand("results/03_benchmark/merged2/02_mean_rank/{PKN}/zscore-{PKN}.csv", PKN = config["figures"]["evaluation"]),
        act = expand("data/results_cptac/overall_performance/actsiteBM/all_kins/actsiteBM_5perThr_{PKN_cptac}_roc_table.rds", PKN_cptac = config["figures"]["evaluation_cptac"]),
        act_kin = expand("data/results_cptac/overall_performance/actsiteBM/all_kins/actsiteBM_5perThr_{PKN_cptac}_roc_kins.rds", PKN_cptac = config["figures"]["evaluation_cptac"]),
        tumor = expand("data/results_cptac/overall_performance/protBM/all_kins/protBM_5perThr_{PKN_cptac}_roc_table.rds", PKN_cptac = config["figures"]["evaluation_cptac"]),
        tumor_kin = expand("data/results_cptac/overall_performance/protBM/all_kins/protBM_5perThr_{PKN_cptac}_roc_kins.rds", PKN_cptac = config["figures"]["evaluation_cptac"]),
        kinclass = "resources/kinase_class.csv"
    output:
        plt="results/manuscript_figures/figure_3/zscore.pdf",
        csv="results/manuscript_figures/supp_files/overview_benchmark.csv"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/03_prior_evaluation.R"

rule performance_strat:
    input:
        rank = expand("results/03_benchmark/merged2/02_mean_rank/{PKN}/zscore-{PKN}.csv", PKN = config["figures"]["evaluation"]),
        act = expand("data/results_cptac/stratified_kinases/all_kins/actsiteBM/actsiteBM_{PKN_cptac}_{strat}_roc_list.Rds", PKN_cptac = config["figures"]["evaluation_cptac"], strat = ["0to10", "11to25", "26plus"]),
        tumor = expand("data/results_cptac/stratified_kinases/all_kins/protBM/protBM_{PKN_cptac}_{strat}_roc_list.Rds", PKN_cptac = config["figures"]["evaluation_cptac"], strat = ["0to10", "11to25", "26plus"])
    output:
        plt_pert="results/manuscript_figures/figure_3/regulon_perturbation.pdf",
        plt_act="results/manuscript_figures/figure_3/regulon_activatingsite.pdf",
        plt_tumor="results/manuscript_figures/figure_3/regulon_tumor.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/03_size_evaluation.R"


rule performance_priors_subset_supp:
    input:
        bench = expand("results/03_benchmark/merged2/02_benchmark_res_subset/subset/{PKN}/bench_zscore-{PKN}.csv", PKN = config["figures"]["coverage"]),
        rank = expand("results/03_benchmark/merged2/02_mean_rank_subset/subset/{PKN}/zscore-{PKN}.csv", PKN = config["figures"]["coverage"]),
        act = expand("data/results_cptac/overall_performance/actsiteBM/same_kins/actsiteBM_5perThr_{PKN_cptac}_roc_filt_table.rds", PKN_cptac = config["figures"]["cptac"]),
        act_kin = expand("data/results_cptac/overall_performance/actsiteBM/same_kins/actsiteBM_5perThr_{PKN_cptac}_roc_filt_kins.rds", PKN_cptac = config["figures"]["cptac"]),
        tumor = expand("data/results_cptac/overall_performance/protBM/same_kins/protBM_5perThr_{PKN_cptac}_roc_filt_table.rds", PKN_cptac = config["figures"]["cptac"]),
        tumor_kin = expand("data/results_cptac/overall_performance/protBM/same_kins/protBM_5perThr_{PKN_cptac}_roc_filt_kins.rds", PKN_cptac = config["figures"]["cptac"])
    output:
        plt="results/manuscript_figures/figure_3/supp/zscore_supp.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/03_prior_evaluation.R"
        
rule prior_coverage_supp:
    input:
        prior_files = expand("results/00_prior/{PKN}.tsv", PKN = config["figures"]["coverage"]),
        class_kin = "resources/kinase_class.csv"
    output:
        upset="results/manuscript_figures/figure_3/supp/upset_kin.pdf",
        upsetEdge="results/manuscript_figures/figure_3/supp/upset_edge.pdf",
        pps="results/manuscript_figures/figure_3/supp/jaccard.pdf",
        kin_type="results/manuscript_figures/figure_3/supp/kinase_type.pdf",
        reg="results/manuscript_figures/figure_3/supp/regulon_size.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/03.1_supp_overview_priors.R"       

rule correlation_analysis_supp:
    input:
        pearson_prior="results/04_exploration/merged2/correlation/correlation_priors_pearson.rds",
        pearson_method="results/04_exploration/merged2/correlation/correlation_methods_pearson.rds",
        spearman_prior="results/04_exploration/merged2/correlation/correlation_priors_spearman.rds",
        spearman_method="results/04_exploration/merged2/correlation/correlation_methods_spearman.rds",
        jaccard_up_prior="results/04_exploration/merged2/jaccard/jaccard_priors_up_10.rds",
        jaccard_down_prior="results/04_exploration/merged2/jaccard/jaccard_priors_down_10.rds",
        jaccard_up_method="results/04_exploration/merged2/jaccard/jaccard_methods_up_10.rds",
        jaccard_down_method="results/04_exploration/merged2/jaccard/jaccard_methods_down_10.rds"
    output:
        cor_plot_pearson_prior="results/manuscript_figures/figure_3/supp/corrplot_pearson_priors.pdf",
        cor_plot_spearman_prior="results/manuscript_figures/figure_3/supp/corrplot_spearman_priors.pdf",
        cor_plot_jacc_up_prior="results/manuscript_figures/figure_3/supp/corrplot_jacc_up_priors.pdf",
        cor_plot_jacc_down_prior="results/manuscript_figures/figure_3/supp/corrplot_jacc_down_priors.pdf",
        cor_plot_pearson_method="results/manuscript_figures/figure_3/supp/corrplot_pearson_methods.pdf",
        cor_plot_spearman_method="results/manuscript_figures/figure_3/supp/corrplot_spearman_methods.pdf",
        cor_plot_jacc_up_method="results/manuscript_figures/figure_3/supp/corrplot_jacc_up_methods.pdf",
        cor_plot_jacc_down_method="results/manuscript_figures/figure_3/supp/corrplot_jacc_down_methods.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/03.1_supp_comparison_activity.R"       
  
rule performance_priors_all_supp:
    input:
        bench_files = expand("results/03_benchmark/merged2/02_benchmark_res/{PKN}/bench_{methods}-{PKN}.csv", PKN = config["figures"]["evaluation"],  methods = config["perturbation"]["methods"]),
        act = expand("data/results_cptac/overall_performance/actsiteBM/all_kins/actsiteBM_5perThr_{PKN_cptac}_roc_table.rds", PKN_cptac = config["figures"]["evaluation_cptac"]),
        tumor = expand("data/results_cptac/overall_performance/protBM/all_kins/protBM_5perThr_{PKN_cptac}_roc_table.rds", PKN_cptac = config["figures"]["evaluation_cptac"])
    output:
        plt = "results/manuscript_figures/figure_3/supp/median_auroc.pdf",
        plt_tumor = "results/manuscript_figures/figure_3/supp/median_auroc_tumor.pdf",
        plt_act = "results/manuscript_figures/figure_3/supp/median_auroc_act.pdf",
        csv_prior = "results/manuscript_figures/supp_files/prior_comparison_perturbation.csv",
        csv_method = "results/manuscript_figures/supp_files/method_comparison_perturbation.csv",
        csv_prior_tumor = "results/manuscript_figures/supp_files/prior_comparison_tumor.csv",
        csv_method_tumor = "results/manuscript_figures/supp_files/method_comparison_tumor.csv",
        csv_prior_act = "results/manuscript_figures/supp_files/prior_comparison_act.csv",
        csv_method_act = "results/manuscript_figures/supp_files/method_comparison_act.csv"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/03.1_supp_combo_evaluation.R"

rule performance_priors_class:
    input:
        rank = expand("results/03_benchmark/merged2/02_mean_rank/{PKN}/zscore-{PKN}.csv", PKN = config["figures"]["evaluation"]),
        kinclass = "resources/kinase_class.csv"
    output:
        plt_ser = "results/manuscript_figures/figure_3/supp/ser_performance.pdf",
        plt_tyr = "results/manuscript_figures/figure_3/supp/tyr_performance.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/03.1_supp_class_evaluation.R"

rule supp_table_comparison:
    input:
        prior = expand("results/manuscript_figures/supp_files/prior_comparison_{bench}.csv", bench = ["act", "perturbation", "tumor"]),
        method = expand("results/manuscript_figures/supp_files/method_comparison_{bench}.csv", bench = ["act", "perturbation", "tumor"])
    output:
        prior_out = "results/manuscript_figures/supp_files/prior_comparison.csv",
        method_out = "results/manuscript_figures/supp_files/method_comparison.csv"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/03.1_supp_combo_table.R"

# ------------------------------------ FIGURE 4 ------------------------------------
rule performance_combinations:
    input:
        bench = expand("results/03_benchmark/merged2/02_benchmark_res/{PKN}/bench_{methods}-{PKN}.csv", PKN = config["figures"]["combination"], methods = config["figures"]["combination_methods"]),
        rank = expand("results/03_benchmark/merged2/02_mean_rank/{PKN}/{methods}-{PKN}.csv", PKN = config["figures"]["combination"], methods = config["figures"]["combination_methods"]),
        act_roc = expand("data/results_cptac/performance_combinations/GSknown/all_kins/actsiteBM_5perThr_combofig_{PKN_cptac}_roc_table.rds", PKN_cptac = config["figures"]["combination_cptac"]),
        act_kin = expand("data/results_cptac/performance_combinations/GSknown/all_kins/actsiteBM_5perThr_combofig_{PKN_cptac}_roc_kins.rds", PKN_cptac = config["figures"]["combination_cptac"]),
        tumor_roc = expand("data/results_cptac/performance_combinations/GSknown/all_kins/protBM_5perThr_combofig_{PKN_cptac}_roc_table.rds", PKN_cptac = config["figures"]["combination_cptac"]),
        tumor_kin = expand("data/results_cptac/performance_combinations/GSknown/all_kins/protBM_5perThr_combofig_{PKN_cptac}_roc_kins.rds", PKN_cptac = config["figures"]["combination_cptac"]),
        kinclass = "resources/kinase_class.csv"
    output:
        plt="results/manuscript_figures/figure_4/combinations_zscore_perturbation.pdf",
        plt_act="results/manuscript_figures/figure_4/combinations_zscore_actsite.pdf",
        plt_prot="results/manuscript_figures/figure_4/combinations_zscore_protein.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/04_combination_evaluation.R"


rule performance_combinations_subset:
    input:
        bench = expand("results/03_benchmark/merged2/02_benchmark_res_subset/GSknownSub/{PKN}/bench_{methods}-{PKN}.csv", PKN = config["figures"]["combination"], methods = config["figures"]["combination_methods"]),
        rank = expand("results/03_benchmark/merged2/02_mean_rank_subset/GSknownSub/{PKN}/{methods}-{PKN}.csv", PKN = config["figures"]["combination"], methods = config["figures"]["combination_methods"]),
        act_roc = expand("data/results_cptac/performance_combinations/GSknown/same_kins/actsiteBM_5perThr_combofig_{PKN_cptac}_roc_filt_table.rds", PKN_cptac = config["figures"]["combination_cptac"]),
        act_kin = "data/results_cptac/performance_combinations/GSknown/same_kins/actsiteBM_5perThr_combofig_known_roc_kins_used4all.rds",
        tumor_roc = expand("data/results_cptac/performance_combinations/GSknown/same_kins/protBM_5perThr_combofig_{PKN_cptac}_roc_filt_table.rds", PKN_cptac = config["figures"]["combination_cptac"]),
        tumor_kin = "data/results_cptac/performance_combinations/GSknown/same_kins/protBM_5perThr_combofig_known_roc_kins_used4all.rds"
    output:
        plt="results/manuscript_figures/figure_4/combinations_subset_zscore_perturbation.pdf",
        plt_act="results/manuscript_figures/figure_4/combinations_subset_zscore_actsite.pdf",
        plt_prot="results/manuscript_figures/figure_4/combinations_subset_zscore_protein.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/04_combination_subset_evaluation.R"

rule performance_combinations_supp:
    input:
        bench = expand("results/03_benchmark/merged2/02_benchmark_res/{PKN}/bench_{methods}-{PKN}.csv", PKN = config["perturbation"]["johnson1Sub"], methods = config["figures"]["combination_methods"]),
        rank = expand("results/03_benchmark/merged2/02_mean_rank/{PKN}/{methods}-{PKN}.csv", PKN = config["perturbation"]["johnson1Sub"], methods = config["figures"]["combination_methods"]),
        act_roc = expand("data/results_cptac/performance_combinations/GSknown/all_kins/actsiteBM_5perThr_combofig_{PKN_cptac}_roc_table.rds", PKN_cptac = config["figures"]["combination_cptac_2"]),
        act_kin = expand("data/results_cptac/performance_combinations/GSknown/all_kins/actsiteBM_5perThr_combofig_{PKN_cptac}_roc_kins.rds", PKN_cptac = config["figures"]["combination_cptac_2"]),
        tumor_roc = expand("data/results_cptac/performance_combinations/GSknown/all_kins/protBM_5perThr_combofig_{PKN_cptac}_roc_table.rds", PKN_cptac = config["figures"]["combination_cptac_2"]),
        tumor_kin = expand("data/results_cptac/performance_combinations/GSknown/all_kins/protBM_5perThr_combofig_{PKN_cptac}_roc_kins.rds", PKN_cptac = config["figures"]["combination_cptac_2"])
    output:
        plt="results/manuscript_figures/figure_4/supp/combinations_zscore_perturbation.pdf",
        plt_act="results/manuscript_figures/figure_4/supp/combinations_zscore_actsite.pdf",
        plt_prot="results/manuscript_figures/figure_4/supp/combinations_zscore_protein.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/04_combination_evaluation.R"


rule performance_combinations_subset_supp:
    input:
        bench = expand("results/03_benchmark/merged2/02_benchmark_res_subset/GSknownSub/{PKN}/bench_{methods}-{PKN}.csv", PKN = config["perturbation"]["johnson1Sub"], methods = config["figures"]["combination_methods"]),
        rank = expand("results/03_benchmark/merged2/02_mean_rank_subset/GSknownSub/{PKN}/{methods}-{PKN}.csv", PKN = config["perturbation"]["johnson1Sub"], methods = config["figures"]["combination_methods"]),
        act_roc = expand("data/results_cptac/performance_combinations/GSknown/same_kins/actsiteBM_5perThr_combofig_{PKN_cptac}_roc_filt_table.rds", PKN_cptac = config["figures"]["combination_cptac_2"]),
        act_kin = "data/results_cptac/performance_combinations/GSknown/same_kins/actsiteBM_5perThr_combofig_known_roc_kins_used4all.rds",
        tumor_roc = expand("data/results_cptac/performance_combinations/GSknown/same_kins/protBM_5perThr_combofig_{PKN_cptac}_roc_filt_table.rds", PKN_cptac = config["figures"]["combination_cptac_2"]),
        tumor_kin = "data/results_cptac/performance_combinations/GSknown/same_kins/protBM_5perThr_combofig_known_roc_kins_used4all.rds"
    output:
        plt="results/manuscript_figures/figure_4/supp/combinations_subset_zscore_perturbation.pdf",
        plt_act="results/manuscript_figures/figure_4/supp/combinations_subset_zscore_actsite.pdf",
        plt_prot="results/manuscript_figures/figure_4/supp/combinations_subset_zscore_protein.pdf"
    conda:
        "../envs/figures.yml"
    script:
        "../scripts/05_figures_manuscript/04_combination_subset_evaluation.R"