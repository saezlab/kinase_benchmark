# ------------------------------------ DATA PROCESSING ------------------------------------
rule protein_mapping:
    input:
        input_folder = "data/decryptm/{decryptm_set}",
        ref_proteome_file = "data/decryptm/uniprot_proteome_up000005640_03112020.fasta"
    output:
        output_file = "results/decryptm/protein_mapping/mapped_protein_{decryptm_set}.csv"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/01_data_processing/decryptm/01_protein_mapping.py"

rule phospho_preprocessing:
    input:
        input_folder = "data/decryptm/{decryptm_set}",
        mapped_file = "results/decryptm/protein_mapping/mapped_protein_{decryptm_set}.csv"
    output:
        phospho = "results/decryptm/phosphoproteome/EC50_{decryptm_set}.csv",
        meta_output = "results/decryptm/phosphoproteome/metadata_EC50_{decryptm_set}.csv"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/01_data_processing/decryptm/02_preprocessing_phospho.R"

rule merge_datasets:
    input:
        EC50_files = expand("results/decryptm/phosphoproteome/EC50_{decryptm_set}.csv", decryptm_set = config["decryptm"]["decryptm_set"]),
        meta_files = expand("results/decryptm/phosphoproteome/metadata_EC50_{decryptm_set}.csv", decryptm_set = config["decryptm"]["decryptm_set"]),
        targets = "data/decryptm/drug_targets.csv"
    output:
        output_meta = "results/decryptm/processed_data/meta_data.csv",
        output_EC50 = "results/decryptm/processed_data/pEC50.csv",
        output_R2_EC50 = "results/decryptm/processed_data/R2_pEC50.csv",
        output_drugs = "results/decryptm/processed_data/overview_drugs.csv"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/01_data_processing/decryptm/03_merge_datasets.R"

rule QC_decryptm:
    input:
        EC50_file = "results/decryptm/processed_data/{EC50}.csv"
    output:
        output_plots = "results/decryptm/processed_data/figures/QC_{EC50}.pdf"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/01_data_processing/decryptm/04_quality_control.R"


# ------------------------------ INPUT PREPARATION ------------------------------
# ------------------------------ Prior knowledge preparation ------------------------------
rule prepare_omnipath:
    input:
        file_dataset = expand("data/CPTAC_phospho/{dataset}_{normalisation}_medcent_30plus.tsv", dataset = config["activity_estimation"]["datasets"], normalisation = config["activity_estimation"]["normalisation"]),
        decryptm = expand("results/decryptm/processed_data/{decryptm_measurement}.csv", decryptm_measurement = config["decryptm"]["decryptm_measurement"])
    output:
        tsv = "results/prior/omnipath.tsv",
        out_decryptm = "results/decryptm/prior/omnipath.tsv"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/02_prior_mapping/01.1_prepare_omnipath.R"

rule prepare_phosphositeplus:
    input:
        ppsp = "data/prior/phosphositeplus",
        file_dataset = expand("data/CPTAC_phospho/{dataset}_{normalisation}_medcent_30plus.tsv", dataset = config["activity_estimation"]["datasets"], normalisation = config["activity_estimation"]["normalisation"]),
        decryptm = expand("results/decryptm/processed_data/{decryptm_measurement}.csv", decryptm_measurement = config["decryptm"]["decryptm_measurement"])
    output:
        tsv = "results/prior/phosphositeplus.tsv",
        out_decryptm = "results/decryptm/prior/phosphositeplus.tsv"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/02_prior_mapping/01.2_prepare_phosphositeplus.R"

rule prepare_ptmsigdb:
    input:
        ptmsig_file = "data/prior/ptm.sig.db.all.uniprot.human.v1.9.0.gmt",
        file_dataset = expand("data/CPTAC_phospho/{dataset}_{normalisation}_medcent_30plus.tsv", dataset = config["activity_estimation"]["datasets"], normalisation = config["activity_estimation"]["normalisation"]),
        decryptm = expand("results/decryptm/processed_data/{decryptm_measurement}.csv", decryptm_measurement = config["decryptm"]["decryptm_measurement"])
    output:
        tsv = "results/prior/ptmsigdb.tsv",
        out_decryptm = "results/decryptm/prior/ptmsigdb.tsv"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/02_prior_mapping/01.3_prepare_ptmsigdb.R"

rule prepare_ikipdb:
    input:
        ptmsig_file = "data/prior/iKiP-DB-Table.tsv",
        file_dataset = expand("data/CPTAC_phospho/{dataset}_{normalisation}_medcent_30plus.tsv", dataset = config["activity_estimation"]["datasets"], normalisation = config["activity_estimation"]["normalisation"]),
        decryptm = expand("results/decryptm/processed_data/{decryptm_measurement}.csv", decryptm_measurement = config["decryptm"]["decryptm_measurement"])
    output:
        tsv = "results/prior/iKiPdb.tsv",
        out_decryptm = "results/decryptm/prior/iKiPdb.tsv"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/02_prior_mapping/01.4_prepare_ikipdb.R"

rule prepare_GPS:
    input:
        GPS_file = "data/prior/mmc4.xlsx",
        file_dataset = expand("data/CPTAC_phospho/{dataset}_{normalisation}_medcent_30plus.tsv", dataset = config["activity_estimation"]["datasets"], normalisation = config["activity_estimation"]["normalisation"]),
        decryptm = expand("results/decryptm/processed_data/{decryptm_measurement}.csv", decryptm_measurement = config["decryptm"]["decryptm_measurement"])
    output:
        tsv = "results/prior/GPS.tsv",
        out_decryptm = "results/decryptm/prior/GPS.tsv"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/02_prior_mapping/01.5_prepare_GPS.R"

rule prepare_NetworKIN:
    input:
        networkin_file = "data/prior/networkin_human_predictions_3.1.tsv",
        file_dataset = expand("data/CPTAC_phospho/{dataset}_{normalisation}_medcent_30plus.tsv", dataset = config["activity_estimation"]["datasets"], normalisation = config["activity_estimation"]["normalisation"]),
        decryptm = expand("results/decryptm/processed_data/{decryptm_measurement}.csv", decryptm_measurement = config["decryptm"]["decryptm_measurement"])
    output:
        tsv = "results/prior/networkin.tsv",
        tsv_decryptm = "results/decryptm/prior/networkin.tsv"
    params:
    	score = 5
    conda:
        "envs/phospho.yml"
    script:
        "scripts/02_prior_mapping/01.6_prepare_NetworKIN.R"

rule prepare_jhonson:
    input:
        ppsp = "data/prior/simplified_jhonson_with_psite.csv",
        file_dataset = expand("data/CPTAC_phospho/{dataset}_{normalisation}_medcent_30plus.tsv", dataset = config["activity_estimation"]["datasets"], normalisation = config["activity_estimation"]["normalisation"]),
        decryptm = expand("results/decryptm/processed_data/{decryptm_measurement}.csv", decryptm_measurement = config["decryptm"]["decryptm_measurement"])
    output:
        tsv = "results/prior/jhonson.tsv",
        out_decryptm = "results/decryptm/prior/jhonson.tsv"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/02_prior_mapping/01.9_prepare_jhonson.R"

rule merge_GPS_PPSP:
    input:
        gps = "results/prior/GPS.tsv",
        ppsp = "results/prior/phosphositeplus.tsv",
        gps_decryptm = "results/decryptm/prior/GPS.tsv",
        ppsp_dectyptm = "results/decryptm/prior/phosphositeplus.tsv"
    output:
        tsv = "results/prior/GPSppsp.tsv",
        tsv_decryptm = "results/decryptm/prior/GPSppsp.tsv"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/02_prior_mapping/01.7_merge_GPS_PPSP.R"

rule merge_known_predicted:
    input:
        known_file = "results/prior/{known_targets}.tsv",
        predicted_file = "results/prior/{predicted_targets}.tsv",
        known_decryptm = "results/decryptm/prior/{known_targets}.tsv",
        predicted_dectyptm = "results/decryptm/prior/{predicted_targets}.tsv"
    output:
        tsv = "results/prior/{known_targets}_{predicted_targets}.tsv",
        tsv_decryptm = "results/decryptm/prior/{known_targets}_{predicted_targets}.tsv"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/02_prior_mapping/01.8_merge_known_predicted.R"

rule map_kinase_ids:
    input:
        prior_files = expand("results/prior/{PKN}.tsv", PKN = config["activity_estimation"]["PKNs"])
    output:
        tsv = "resources/kinase_mapping.tsv"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/02_prior_mapping/02_convert_kinase_ids.R"

# ------------------------------- PTM-SEA input preparation -------------------------------
rule ptmsea_datasets:
    input:
        file_dataset = "data/CPTAC_phospho/{dataset}_{normalisation}_medcent_30plus.tsv"
    output:
        gct = "results/datasets/{dataset}_{normalisation}.gct"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/02_prior_mapping/02.1_prepare_datasets_ptm-sea.R"

rule ptmsea_prior:
    input:
        file_PKN = "results/prior/{PKN}.tsv"
    output:
        gmt = "results/prior/ptm-sea/{PKN}.gmt"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/02_prior_mapping/02.2_prepare_prior_ptm-sea.R"


# ---------------------------------- ACTIVITY ESTIMATION ----------------------------------
rule activity_estimation_decryptm:
    input:
        file_dataset = "results/decryptm/processed_data/{decryptm_measurement}.csv",
        file_PKN = "results/decryptm/prior/{PKN}.tsv",
        scripts = expand("workflow/scripts/methods/run_{method}.R", method = ["INKA", "KARP", "lm_rokai", "zscore", "erics_methods"]),
        script_support = "workflow/scripts/methods/support_functions.R"
    output:
        rds = "results/decryptm/activity_scores/{decryptm_measurement}-{PKN}.rds"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/03_kinase_activity_estimation/decryptm/03_activity_estimation_decryptm.R"


# -------------------------------------- BENCHMARK ---------------------------------------
rule prepare_decryptm_benchmark:
    input:
        rds = "results/decryptm/activity_scores/{decryptm_measurement}-{PKN}.rds",
        meta = "results/decryptm/processed_data/meta_data.csv",
        kinome = "data/decryptm/decryptm_processed_targets.csv"
    output:
        output = "results/decryptm/benchmark_scores/{decryptm_methods}-{decryptm_measurement}-{PKN}.csv",
        meta_out = "results/decryptm/benchmark_scores/obs_{decryptm_methods}-{decryptm_measurement}-{PKN}.csv"
    params:
    	  perturb = "kinomebeads"
    conda:
        "envs/phospho.yml"
    script:
        "scripts/04_decryptm_benchmark/01_prepare_bench_input.R"

rule run_decryptm_benchmark:
    input:
        scores = "results/decryptm/benchmark_scores/{decryptm_methods}-{decryptm_measurement}-{PKN}.csv",
        meta = "results/decryptm/benchmark_scores/obs_{decryptm_methods}-{decryptm_measurement}-{PKN}.csv"
    output:
        output = "results/decryptm/benchmark/{PKN}/bench_{decryptm_methods}-{decryptm_measurement}-{PKN}.csv"
    conda:
        "envs/benchmark.yml"
    script:
        "scripts/04_decryptm_benchmark/02_decouple_bench.py"