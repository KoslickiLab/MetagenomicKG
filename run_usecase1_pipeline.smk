"""
This file is a SnakeMake Script to automate model training for pathogen identification as use case1.

Usage:
    snakemake --cores 16 -s run_usecase1_pipeline.smk targets
"""
## Import Config Files
configfile: "./config.yaml"

## Import Python standard libraries
import os, sys
import subprocess

## Define Some Global Variables
CURRENT_PATH = '/scratch/shared_data_new/metagenomics_graph/model_training'
main_dir = '/scratch/shared_data_new/metagenomics_graph'
EXPERIMENT_NAME = "pathogen_detection"

ROOT_PATH = os.getcwd()
DATA_PATH = os.path.join(ROOT_PATH, 'data')
GNN_PROCESSED_DATA_PATH = os.path.join(DATA_PATH, "usecase1_processed_data", 'GNN'+'_'+EXPERIMENT_NAME)
GNN_RESULTS_PATH = os.path.join(DATA_PATH, "usecase1_results", 'GNN'+'_'+EXPERIMENT_NAME)
GNN_SCRIPT_PATH = os.path.join(ROOT_PATH, 'usecase1_pathogen_identification', 'GNN')
EXPERIMENT_NAME = config['USECASE1_VARIABLES']['EXPERIMENT_NAME']

## Create Required Folders
if not os.path.exists(GNN_PROCESSED_DATA_PATH):
    os.makedirs(GNN_PROCESSED_DATA_PATH)
if not os.path.exists(GNN_RESULTS_PATH):
    os.makedirs(GNN_RESULTS_PATH)

## Build Rules
rule targets:
    input:
        os.path.join(GNN_PROCESSED_DATA_PATH, "node_info.tsv"),
        os.path.join(GNN_PROCESSED_DATA_PATH, "node_to_index.tsv"),
        os.path.join(GNN_PROCESSED_DATA_PATH, "pathogen_info.tsv"),
        os.path.join(GNN_PROCESSED_DATA_PATH, "text_embedding", "id2index.json"),
        os.path.join(GNN_PROCESSED_DATA_PATH, "text_embedding", "index2id.json"),
        os.path.join(GNN_PROCESSED_DATA_PATH, "text_embedding", "embedding_biobert_namecat.pkl"),
        os.path.join(GNN_RESULTS_PATH, "performance_summary_test.tsv")

##############################  GNN Model Training Pipeline  ##############################

# Prepare data for model training
rule step0_prepare_data:
    input:
        script = ancient(os.path.join(GNN_SCRIPT_PATH, "model_training_preparation.py")),
        existing_KG_nodes = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v6.tsv')),
        existing_KG_edges = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v6.tsv')),
        micobial_hierarchy_dir = ancient(os.path.join(DATA_PATH, "Micobial_hierarchy")),
        output_dir = ancient(GNN_PROCESSED_DATA_PATH)
    output:
        os.path.join(GNN_PROCESSED_DATA_PATH, "node_info.tsv"),
        os.path.join(GNN_PROCESSED_DATA_PATH, "node_to_index.tsv"),
        os.path.join(GNN_PROCESSED_DATA_PATH, "pathogen_info.tsv"),
        os.path.join(GNN_PROCESSED_DATA_PATH, "edge_info.tsv")
    run:
        shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --micobial_hierarchy_dir {input.micobial_hierarchy_dir} --output_dir {input.output_dir}")

# Convert node description and node type to embeddings
rule step1_convert_node_description_to_embeddings:
    input:
        script = ancient(os.path.join(GNN_SCRIPT_PATH, "calculate_node_embedding.py")),
        node_info = ancient(os.path.join(GNN_PROCESSED_DATA_PATH, "node_info.tsv")),
        node_to_index = ancient(os.path.join(GNN_PROCESSED_DATA_PATH, "node_to_index.tsv")),
        output_dir = ancient(GNN_PROCESSED_DATA_PATH)
    params:
        gpu = 0,
        final_embedding_dim = 200,
        batch_size = 100
    output:
        os.path.join(GNN_PROCESSED_DATA_PATH, "text_embedding", "id2index.json"),
        os.path.join(GNN_PROCESSED_DATA_PATH, "text_embedding", "index2id.json"),
        os.path.join(GNN_PROCESSED_DATA_PATH, "text_embedding", "embedding_biobert_namecat.pkl")
    run:
        shell("python {input.script} --node_info {input.node_info} --node_to_index {input.node_to_index} --gpu {params.gpu} --use_gpu --final_embedding_dim {params.final_embedding_dim} --batch_size {params.batch_size} --output_dir {input.output_dir}")

# Train the model
rule step2_train_model:
    input:
        script = ancient(os.path.join(GNN_SCRIPT_PATH, "run_model.py")),
        data_path = ancient(GNN_PROCESSED_DATA_PATH),
        output_dir = ancient(GNN_RESULTS_PATH)
    params:
        experiment_name = EXPERIMENT_NAME,
        gpu = 0,
        use_gpu = True,
        ratio = 1.0,
        cv_num = 10,
        lr = 0.001,
        num_neighbors = [100, 100, 100],
        num_workers = 100,
        num_epochs = 500,
        train_batch_size = 256,
        eval_batch_size = 128,
        emb_size = 512,
        early_stop_n = 500,
        num_layers = 3
    output:
        os.path.join(GNN_RESULTS_PATH, "performance_summary_test.tsv")
    run:
        if params.use_gpu:
            shell("python {input.script} --data_path {input.data_path} --output_dir {input.output_dir} --experiment_name {params.experiment_name} --gpu {params.gpu} --use_gpu --ratio {params.ratio} --cv_num {params.cv_num} --lr {params.lr} --num_neighbors {params.num_neighbors} --num_workers {params.num_workers} --num_epochs {params.num_epochs} --train_batch_size {params.train_batch_size} --eval_batch_size {params.eval_batch_size} --emb_size {params.emb_size} --num_layers {params.num_layers} --early_stop_n {params.early_stop_n}")
        else:
            shell("python {input.script} --data_path {input.data_path} --output_dir {input.output_dir} --experiment_name {params.experiment_name} --gpu {params.gpu} --ratio {params.ratio} --cv_num {params.cv_num} --lr {params.lr} --num_neighbors {params.num_neighbors} --num_workers {params.num_workers} --num_epochs {params.num_epochs} --train_batch_size {params.train_batch_size} --eval_batch_size {params.eval_batch_size} --emb_size {params.emb_size} --num_layers {params.num_layers} --early_stop_n {params.early_stop_n}")

##############################  Other baseline Model Training Pipeline  ##############################