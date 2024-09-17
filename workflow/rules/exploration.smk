rule test_performance:
    input:
        scores = "results/03_benchmark/{dataset}/02_mean_rank/{PKN}/{method}-{PKN}.csv"
    output:
        out = "results/03_benchmark/{dataset}/02_mean_rank/plots/{method}-{PKN}_targets.pdf"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/04_exploration/test_performance.R"
