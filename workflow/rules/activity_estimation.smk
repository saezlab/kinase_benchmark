# Define common parameters
rm_auto = True
minsize = 5
threads = 3

rule run_chisq:
    input:
        file_dataset="results/01_processed_data/{dataset}/data/benchmark_data.csv",
        file_PKN="results/01_processed_data/{dataset}/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_chisq.R",
        script_support="workflow/scripts/methods/support_functions.R"
    output:
        rds="results/02_activity_scores/{dataset}/chisq/{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize,
        background=20000
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_chisq.R"
        
rule run_fgsea:
    input:
        file_dataset="results/01_processed_data/{dataset}/data/benchmark_data.csv",
        file_PKN="results/01_processed_data/{dataset}/mapped_priors/{PKN}.tsv"
    output:
        rds="results/02_activity_scores/{dataset}/fgsea/{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_fgsea.R"

rule run_fisher:
    input:
        file_dataset="results/01_processed_data/{dataset}/data/benchmark_data.csv",
        file_PKN="results/01_processed_data/{dataset}/mapped_priors/{PKN}.tsv"
    output:
        rds="results/02_activity_scores/{dataset}/fisher/{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize,
        background=20000
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_fisher.R"

rule run_INKA:
    input:
        file_dataset="results/01_processed_data/{dataset}/data/benchmark_data.csv",
        file_PKN="results/01_processed_data/{dataset}/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_INKA.R",
        script_support="workflow/scripts/methods/support_functions.R"
    output:
        rds="results/02_activity_scores/{dataset}/INKA/{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_INKA.R"

rule run_KARP:
    input:
        file_dataset="results/01_processed_data/{dataset}/data/benchmark_data.csv",
        file_PKN="results/01_processed_data/{dataset}/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_KARP.R",
        script_support="workflow/scripts/methods/support_functions.R"
    output:
        rds="results/02_activity_scores/{dataset}/KARP/{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_KARP.R"

rule run_KSEA:
    input:
        file_dataset="results/01_processed_data/{dataset}/data/benchmark_data.csv",
        file_PKN="results/01_processed_data/{dataset}/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_zscore.R",
        script_support="workflow/scripts/methods/support_functions.R"
    output:
        rds="results/02_activity_scores/{dataset}/KSEA/{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_KSEA.R"

rule run_lmRoKAI:
    input:
        file_dataset="results/01_processed_data/{dataset}/data/benchmark_data.csv",
        file_PKN="results/01_processed_data/{dataset}/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_lm_rokai.R",
        script_support="workflow/scripts/methods/support_functions.R"
    output:
        rds="results/02_activity_scores/{dataset}/lmRoKAI/{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_lmRoKAI.R"

rule run_mean:
    input:
        file_dataset="results/01_processed_data/{dataset}/data/benchmark_data.csv",
        file_PKN="results/01_processed_data/{dataset}/mapped_priors/{PKN}.tsv"
    output:
        rds="results/02_activity_scores/{dataset}/mean/{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_mean.R"

rule run_misc:
    input:
        file_dataset="results/01_processed_data/{dataset}/data/benchmark_data.csv",
        file_PKN="results/01_processed_data/{dataset}/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_erics_methods.R"
    output:
        rds="results/02_activity_scores/{dataset}/misc/{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_misc.R"

rule run_mlm:
    input:
        file_dataset="results/01_processed_data/{dataset}/data/benchmark_data.csv",
        file_PKN="results/01_processed_data/{dataset}/mapped_priors/{PKN}.tsv"
    output:
        rds="results/02_activity_scores/{dataset}/mlm/{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_mlm.R"

rule run_ptmsea:
    input:
        file_dataset="results/01_processed_data/{dataset}/datasets/benchmark.gct",
        file_PKN="results/01_processed_data/{dataset}/mapped_priors/ptm-sea/{PKN}.gmt"
    output:
        log=temp("results/02_activity_scores/{dataset}/ptmsea/res/log/{PKN}.log"),
        rds="results/02_activity_scores/{dataset}/ptmsea/{PKN}.csv"
    params:
        output_folder=temp("results/02_activity_scores/{dataset}/ptmsea/res"),
        minsize=minsize
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_ptmsea.R"

rule run_ulm:
    input:
        file_dataset="results/01_processed_data/{dataset}/data/benchmark_data.csv",
        file_PKN="results/01_processed_data/{dataset}/mapped_priors/{PKN}.tsv"
    output:
        rds="results/02_activity_scores/{dataset}/ulm/{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_ulm.R"

rule run_viper:
    input:
        file_dataset="results/01_processed_data/{dataset}/data/benchmark_data.csv",
        file_PKN="results/01_processed_data/{dataset}/mapped_priors/{PKN}.tsv"
    output:
        rds="results/02_activity_scores/{dataset}/viper/{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_viper.R"

rule run_zscore:
    input:
        file_dataset="results/01_processed_data/{dataset}/data/benchmark_data.csv",
        file_PKN="results/01_processed_data/{dataset}/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_zscore.R",
        script_support="workflow/scripts/methods/support_functions.R"
    output:
        rds="results/02_activity_scores/{dataset}/zscore/{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_zscore.R"

rule combine_scores:
    input:
        file_scores=expand("results/02_activity_scores/{{dataset}}/{method}/{{PKN}}.csv", method = config["perturbation"]["methods"]),
    output:
        rds="results/02_activity_scores/{dataset}/scores/{PKN}.rds"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/02_combine_scores.R"

