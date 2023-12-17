# ------------------------------ PRIOR PREPARATION ------------------------------
rule prepare_omnipath:
    output:
        tsv = "results/prior/raw/omnipath.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/processing/01.1_prepare_omnipath.R"

rule prepare_phosphositeplus:
    input:
        ppsp = "data/prior/phosphositeplus"
    output:
        tsv = "results/prior/raw/phosphositeplus.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/processing/01.2_prepare_phosphositeplus.R"

rule prepare_ptmsigdb:
    input:
        ptmsig = "data/prior/ptm.sig.db.all.uniprot.human.v1.9.0.gmt"
    output:
        tsv = "results/prior/raw/ptmsigdb.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/processing/01.3_prepare_ptmsigdb.R"

rule prepare_ikipdb:
    input:
        ikip = "data/prior/iKiP-DB-Table.tsv"
    output:
        tsv = "results/prior/raw/iKiPdb.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/processing/01.4_prepare_ikipdb.R"

rule prepare_GPS:
    input:
        GPS_file = "data/prior/mmc4.xlsx"
    output:
        tsv = "results/prior/raw/GPS.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/processing/01.5_prepare_GPS.R"

rule prepare_NetworKIN:
    input:
        networkin_file = "data/prior/networkin_human_predictions_3.1.tsv"
    output:
        tsv = "results/prior/raw/networkin.tsv"
    params:
    	score = 5
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/processing/01.6_prepare_NetworKIN.R"

# ------------------------------ MERGE PRIOR ------------------------------
rule merge_GPS_PPSP:
    input:
        gps = "results/prior/raw/GPS.tsv",
        ppsp = "results/prior/raw/phosphositeplus.tsv",
        ptmsig = "results/prior/raw/ptmsigdb.tsv"
    output:
        merged = "results/prior/raw/GSknown.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/processing/02.1_merge_GPS_PPSP.R"

rule merge_known_predicted:
    input:
        known_file = "results/prior/{known_targets}.tsv",
        predicted_file = "results/prior/{predicted_targets}.tsv"
    output:
        tsv = "results/prior/merged/{known_targets}_{predicted_targets}.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/processing/02.2_merge_known_predicted.R"

# ------------------------------ FILTER PRIOR ------------------------------
rule filter_prior:
    input:
        raw = "results/prior/raw/{prior}.tsv",
        kinhub = "data/kinase_list_kinhub.txt",
        go = "data/QuickGO_kinase_20231214.tsv"
    output:
        filter = "results/prior/{prior}.tsv"
    conda:
        "../envs/phospho.yml"
    script:
        "../scripts/02_prior_mapping/processing/03_filter_priors.R"


