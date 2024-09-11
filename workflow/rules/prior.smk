# This rule is to prepare all kinase-substrate libraries and bring them in the
# same format for mapping to the phosphorylation sites in any dataset

# ------------------------------ PRIOR PREPARATION ------------------------------
rule prepare_omnipath:
    output:
        tsv = temp("results/00_prior/raw/omnipath.tsv")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/01.1_prepare_omnipath.R"

rule prepare_phosphositeplus:
    input:
        ppsp = "data/kinase_libraries/prior/phosphositeplus"
    output:
        tsv = temp("results/00_prior/raw/phosphositeplus.tsv")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/01.2_prepare_phosphositeplus.R"

rule prepare_ptmsigdb:
    input:
        ptmsig = "data/kinase_libraries/prior/ptm.sig.db.all.uniprot.human.v2.0.0.gmt"
    output:
        tsv = temp("results/00_prior/raw/ptmsigdb.tsv")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/01.3_prepare_ptmsigdb.R"

rule prepare_ikipdb:
    input:
        ikip = "data/kinase_libraries/prior/iKiP-DB-Table.tsv"
    output:
        tsv = temp("results/00_prior/raw/iKiPdb.tsv")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/01.4_prepare_ikipdb.R"

rule prepare_GPS:
    input:
        GPS_file = "data/kinase_libraries/prior/mmc4.xlsx"
    output:
        tsv = temp("results/00_prior/raw/GPS.tsv")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/01.5_prepare_GPS.R"

rule prepare_NetworKIN:
    input:
        networkin_file = "data/kinase_libraries/prior/networkin_human_predictions_3.1.tsv"
    output:
        tsv = temp("results/00_prior/raw/networkin.tsv")
    params:
    	score = 5
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/01.6_prepare_NetworKIN.R"

rule prepare_johnson:
    input:
        ppsp = "data/kinase_libraries/johnson_library/pps_fifteenmer_ser_thr_percent.tsv",
        tyr = "data/kinase_libraries/johnson_library/pps_fifteenmer_tyr_percent.tsv"
    params:
        perc = lambda w: w.perc
    output:
        tsv = temp("results/00_prior/raw/johnson{perc}.tsv")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/01.7_prepare_johnson.R"

rule prepare_phosformer:
    input:
        ppsp = "data/kinase_libraries/phosformer/phosformer_results_{nkin}.csv"
    output:
        tsv = temp("results/00_prior/raw/phosformer{nkin}.tsv")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/01.8_prepare_phosformer.R"

rule shuffle_net:
    input:
        tsv = "results/01_processed_data/{dataset}/mapped_priors/phosphositeplus.tsv"
    output:
        shuffled = "results/01_processed_data/{dataset}/mapped_priors/shuffled.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/01.9_shuffle_network.R"

rule shuffle_NetworKIN:
    input:
        tsv = "results/00_prior/merged/phosphositeplus_networkin.tsv"
    output:
        shuffled = "results/00_prior/shuffledNET.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/01.9_shuffle_network.R"

rule shuffle_iKiP:
    input:
        tsv = "results/00_prior/merged/phosphositeplus_iKiPdb.tsv"
    output:
        shuffled = "results/00_prior/shuffledKIP.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/01.9_shuffle_network.R"

# ------------------------------ MERGE PRIOR ------------------------------
rule merge_GPS_PPSP:
    input:
        gps = "results/00_prior/raw/GPS.tsv",
        ppsp = "results/00_prior/raw/phosphositeplus.tsv",
        ptmsig = "results/00_prior/raw/ptmsigdb.tsv"
    output:
        merged = temp("results/00_prior/raw/GSknown.tsv")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/02.1_merge_GPS_PPSP.R"

rule merge_GPS_PPSP_omni:
    input:
        gps = "results/00_prior/raw/GPS.tsv",
        ppsp = "results/00_prior/raw/phosphositeplus.tsv",
        ptmsig = "results/00_prior/raw/ptmsigdb.tsv",
        omnipath = "results/00_prior/raw/omnipath.tsv"
    output:
        merged = temp("results/00_prior/raw/combined.tsv")
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/02.1_merge_GPS_PPSP_omni.R"

rule merge_known_predicted:
    input:
        known_file = "results/00_prior/{known_targets}.tsv",
        predicted_file = "results/00_prior/{predicted_targets}.tsv"
    output:
        tsv = "results/00_prior/merged/{known_targets}_{predicted_targets}.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/02.2_merge_known_predicted.R"

# ------------------------------ FILTER PRIOR ------------------------------
rule kinase_list:
    input:
        kinhub = "data/misc/kinase_list_kinhub.txt",
        go = "data/misc/QuickGO_kinase_20231214.tsv"
    output:
        filter = "data/misc/kinase_list.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/03_generate_kinaseList.R"


rule filter_prior:
    input:
        raw = "results/00_prior/raw/{prior}.tsv",
        kin = "data/misc/kinase_list.tsv"
    output:
        filter = "results/00_prior/{prior}.tsv"
    wildcard_constraints:
        prior = '(?!shuffled)[a-zA-Z0-9]+'
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/00_prior_processing/03_filter_priors.R"


