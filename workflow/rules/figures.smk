# ------------------------------------ FIGURE 1 ------------------------------------
rule overview_prior:
    input:
        prior_files = expand("results/prior/{PKN}.tsv", PKN = config["activity_estimation"]["PKNs"])
    output:
        kin = "results/manuscript_figures/figure_1/coverage_kinases.pdf",
        edges = "results/manuscript_figures/figure_1/coverage_edges.pdf",
        pps = "results/manuscript_figures/figure_1/coverage_pps.pdf",
        kin_heat = "results/manuscript_figures/figure_1/kinase_heatmap.pdf"
    params:
        plot_width = "6",
        plot_height = "4"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/figures_manuscript/overview_priors.R"
