# Define common parameters
rm_auto = True
minsize = 5
threads = 3

rule run_chisq:
    input:
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_chisq.R",
        script_support="workflow/scripts/methods/support_functions.R"
    output:
        rds="results/02_activity_scores/cptac/chisq/{normalisation}/{dataset}-{PKN}.csv"
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
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv"
    output:
        rds="results/02_activity_scores/cptac/fgsea/{normalisation}/{dataset}-{PKN}.csv"
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
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv"
    output:
        rds="results/02_activity_scores/cptac/fisher/{normalisation}/{dataset}-{PKN}.csv"
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
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_INKA.R",
        script_support="workflow/scripts/methods/support_functions.R"
    output:
        rds="results/02_activity_scores/cptac/INKA/{normalisation}/{dataset}-{PKN}.csv"
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
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_KARP.R",
        script_support="workflow/scripts/methods/support_functions.R"
    output:
        rds="results/02_activity_scores/cptac/KARP/{normalisation}/{dataset}-{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_KARP.R"

rule run_KS:
    input:
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_erics_methods.R"
    output:
        rds="results/02_activity_scores/cptac/KS/{normalisation}/{dataset}-{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_KS.R"

rule run_KSEA:
    input:
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_zscore.R",
        script_support="workflow/scripts/methods/support_functions.R"
    output:
        rds="results/02_activity_scores/cptac/KSEA/{normalisation}/{dataset}-{PKN}.csv"
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
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_lm_rokai.R",
        script_support="workflow/scripts/methods/support_functions.R"
    output:
        rds="results/02_activity_scores/cptac/lmRoKAI/{normalisation}/{dataset}-{PKN}.csv"
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
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_erics_methods.R"
    output:
        rds="results/02_activity_scores/cptac/mean/{normalisation}/{dataset}-{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_mean.R"

rule run_median:
    input:
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_erics_methods.R"
    output:
        rds="results/02_activity_scores/cptac/median/{normalisation}/{dataset}-{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_median.R"

rule run_mlm:
    input:
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv"
    output:
        rds="results/02_activity_scores/cptac/mlm/{normalisation}/{dataset}-{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_mlm.R"

rule run_norm:
    input:
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv"
    output:
        rds="results/02_activity_scores/cptac/norm_mean/{normalisation}/{dataset}-{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_norm.R"

rule run_number:
    input:
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_erics_methods.R"
    output:
        rds="results/02_activity_scores/cptac/number_of_targets/{normalisation}/{dataset}-{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_number.R"

rule run_PCA:
    input:
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_erics_methods.R"
    output:
        rds="results/02_activity_scores/cptac/PCA/{normalisation}/{dataset}-{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_PCA.R"

rule run_ptmsea:
    input:
        file_dataset="results/01_processed_data/cptac/datasets/{dataset}_{normalisation}.gct",
        file_PKN="results/01_processed_data/cptac/mapped_priors/ptm-sea/{PKN}.gmt"
    output:
        log=temp("results/02_activity_scores/cptac/ptmsea/res/{normalisation}/{dataset}/log/{dataset}-{PKN}.log"),
        rds="results/02_activity_scores/cptac/ptmsea/{normalisation}/{dataset}-{PKN}.csv"
    params:
        output_folder=temp("results/02_activity_scores/cptac/ptmsea/res/{normalisation}/{dataset}"),
        minsize=minsize
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_ptmsea.R"

rule run_ulm:
    input:
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv"
    output:
        rds="results/02_activity_scores/cptac/ulm/{normalisation}/{dataset}-{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_ulm.R"

rule run_UQ:
    input:
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_erics_methods.R"
    output:
        rds="results/02_activity_scores/cptac/UQ/{normalisation}/{dataset}-{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_UQ.R"

rule run_viper:
    input:
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv"
    output:
        rds="results/02_activity_scores/cptac/viper/{normalisation}/{dataset}-{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_viper.R"

rule run_wilcox:
    input:
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_erics_methods.R"
    output:
        rds="results/02_activity_scores/cptac/wilcox/{normalisation}/{dataset}-{PKN}.csv"
    params:
        rm_auto=rm_auto,
        minsize=minsize
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_wilcox.R"

rule run_zscore:
    input:
        file_dataset= "results/01_processed_data/cptac/data/{normalisation}/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_zscore.R",
        script_support="workflow/scripts/methods/support_functions.R"
    output:
        rds="results/02_activity_scores/cptac/zscore/{normalisation}/{dataset}-{PKN}.csv"
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
        file_scores=expand("results/02_activity_scores/cptac/{method}/{{normalisation}}/{{dataset}}-{{PKN}}.csv", method = config["perturbation"]["methods"]),
    output:
        rds="results/02_activity_scores/cptac/scores/{normalisation}/{dataset}-{PKN}.rds"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/02_combine_scores.R"


## Run with autophosphorylation
rule run_zscore_auto:
    input:
        file_dataset= "results/01_processed_data/cptac/data/original/{dataset}.csv",
        file_PKN="results/01_processed_data/cptac/mapped_priors/{PKN}.tsv",
        scripts="workflow/scripts/methods/run_zscore.R",
        script_support="workflow/scripts/methods/support_functions.R"
    output:
        rds="results/02_activity_scores_auto/cptac/zscore/original/{dataset}-{PKN}.csv"
    params:
        rm_auto=False,
        minsize=minsize
    threads:
        threads
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/01_zscore.R"

rule combine_scores_auto:
    input:
        file_scores=expand("results/02_activity_scores_auto/cptac/{method}/original/{{dataset}}-{{PKN}}.csv", method = ["zscore"]),
    output:
        rds="results/02_activity_scores_auto/cptac/scores/original/{dataset}-{PKN}.rds"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_kinase_activity_estimation/02_combine_scores.R"