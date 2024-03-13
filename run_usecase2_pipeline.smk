"""
This file is a SnakeMake Script to automate generating KG-based metagenomic sample embeddings as use case2.

Usage:
    snakemake --cores 16 -s run_usecase2_pipeline.smk targets
"""
## Import Config Files
ROOT_PATH = os.getcwd()
configfile: f"{ROOT_PATH}/config.yml"

## Import Python standard libraries
import os, sys
import subprocess

DATA_PATH = os.path.join(ROOT_PATH, 'data')
USECASE2_DATA_PATH = os.path.join(DATA_PATH, 'usecase2')
USECASE2_SCRIPT_PATH = os.path.join(ROOT_PATH, 'usecase2_embeddings')
USECASE2_INPUT_PATH = os.path.join(USECASE2_DATA_PATH, "input_data")
USECASE2_RESULTS_PATH = os.path.join(USECASE2_DATA_PATH, "results")

## Create Required Folders
if not os.path.exists(USECASE2_INPUT_PATH):
    os.makedirs(USECASE2_INPUT_PATH)
if not os.path.exists(USECASE2_RESULTS_PATH):
    os.makedirs(USECASE2_RESULTS_PATH)

## Build Rules
rule targets:
    input:
        os.path.join(USECASE2_RESULTS_PATH, "sample_embeddings.csv"),
        os.path.join(USECASE2_RESULTS_PATH, "PCoA_results.png")

# Generate KG-based metagenomic sample embeddings
rule generate_kg_sample_embeddings:
    input:
        script = ancient(os.path.join(USECASE2_SCRIPT_PATH, "generate_kg_sample_embeddings.py")),
        existing_KG_nodes = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v6.tsv')),
        existing_KG_edges = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v6.tsv')),
        abund_metagenomic_samples = ancient(os.path.join(USECASE2_INPUT_PATH, "subsample_100_funiqweighted_fillna0.csv")),
        output_dir = ancient(USECASE2_RESULTS_PATH)
    params:
        alpha = 0.95,
        n_jobs = 32
    output:
        os.path.join(USECASE2_RESULTS_PATH, "sample_embeddings.csv")
    run:
        shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --abund_metagenomic_samples {input.abund_metagenomic_samples} --alpha {params.alpha} --n_jobs {params.n_jobs} --output_dir {input.output_dir}")

# Use tSNE to visualize the embeddings
rule visualize_kg_sample_embeddings:
    input:
        script = ancient(os.path.join(USECASE2_SCRIPT_PATH, "visualize_kg_sample_embeddings.py")),
        metadata = ancient(os.path.join(USECASE2_INPUT_PATH, "metadata_bodysite_subsample.tsv")),
        raw_abundance = ancient(os.path.join(USECASE2_INPUT_PATH, "subsample_100_funiqweighted_fillna0.csv")),
        embeddings = os.path.join(USECASE2_RESULTS_PATH, "sample_embeddings.csv"),
        output_dir = ancient(USECASE2_RESULTS_PATH)
    params:
        top_values = 1000
    output:
        os.path.join(USECASE2_RESULTS_PATH, "PCoA_results.png")
    run:
        shell("python {input.script} --metadata {input.metadata} --raw_abundance {input.raw_abundance} --embeddings {input.embeddings} --top_values {params.top_values} --output_dir {input.output_dir}")
