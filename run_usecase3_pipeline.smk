"""
This file is a SnakeMake Script to automate model training for pathogen identification as use case3.

Usage:
    snakemake --cores 16 -s run_usecase3_pipeline.smk targets
"""
## Import Config Files
ROOT_PATH = os.getcwd()
configfile: f"{ROOT_PATH}/config.yml"

## Import Python standard libraries
import os, sys
import subprocess

## Define Some Global Variables
EXPERIMENT_NAME = config['USECASE3_VARIABLES']['EXPERIMENT_NAME']
DATA_PATH = os.path.join(ROOT_PATH, 'data')
USECASE3_DATA_PATH = os.path.join(DATA_PATH, 'usecase3')
GNN_PROCESSED_DATA_PATH = os.path.join(USECASE3_DATA_PATH, "processed_data", 'GNN'+'_'+EXPERIMENT_NAME)
GNN_RESULTS_PATH = os.path.join(USECASE3_DATA_PATH, "results", 'GNN'+'_'+EXPERIMENT_NAME)
GNN_SCRIPT_PATH = os.path.join(ROOT_PATH, 'usecase3_pathogen_identification', 'GNN')

## Create Required Folders
if not os.path.exists(GNN_PROCESSED_DATA_PATH):
    os.makedirs(GNN_PROCESSED_DATA_PATH)
if not os.path.exists(GNN_RESULTS_PATH):
    os.makedirs(GNN_RESULTS_PATH)
for index in range(10):
    if not os.path.exists(f"{GNN_RESULTS_PATH}_fold{index+1}"):
        os.makedirs(f"{GNN_RESULTS_PATH}_fold{index+1}")

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

# # Train the model with 10-fold cross-validation
# rule step3_train_model_10CV:
#     input:
#         script = ancient(os.path.join(GNN_SCRIPT_PATH, "run_model.py")),
#         data_path = ancient(GNN_PROCESSED_DATA_PATH),
#         output_dir = ancient(GNN_RESULTS_PATH)
#     params:
#         experiment_name = EXPERIMENT_NAME,
#         gpu = 0,
#         use_gpu = True,
#         ratio = 1.0,
#         cv_num = 10,
#         lr = 0.001,
#         num_workers = 100,
#         num_epochs = 500,
#         train_batch_size = 256,
#         eval_batch_size = 128,
#         emb_size = 512,
#         early_stop_n = 500,
#         num_layers = 3
#     output:
#         os.path.join(f"{GNN_RESULTS_PATH}_fold1", "performance_summary_test.tsv"),
#         os.path.join(f"{GNN_RESULTS_PATH}_fold2", "performance_summary_test.tsv"),
#         os.path.join(f"{GNN_RESULTS_PATH}_fold3", "performance_summary_test.tsv"),
#         os.path.join(f"{GNN_RESULTS_PATH}_fold4", "performance_summary_test.tsv"),
#         os.path.join(f"{GNN_RESULTS_PATH}_fold5", "performance_summary_test.tsv"),
#         os.path.join(f"{GNN_RESULTS_PATH}_fold6", "performance_summary_test.tsv"),
#         os.path.join(f"{GNN_RESULTS_PATH}_fold7", "performance_summary_test.tsv"),
#         os.path.join(f"{GNN_RESULTS_PATH}_fold8", "performance_summary_test.tsv"),
#         os.path.join(f"{GNN_RESULTS_PATH}_fold9", "performance_summary_test.tsv"),
#         os.path.join(f"{GNN_RESULTS_PATH}_fold10", "performance_summary_test.tsv")
#     run:
#         for index in range(10):
#             if params.use_gpu:
#                 shell(f"python {input.script} --data_path {input.data_path}/random_10cv/{index+1} --output_dir {input.output_dir}_fold{index+1} --experiment_name {params.experiment_name} --gpu {params.gpu} --use_gpu --ratio {params.ratio} --cv_num {params.cv_num} --lr {params.lr} --num_workers {params.num_workers} --num_epochs {params.num_epochs} --train_batch_size {params.train_batch_size} --eval_batch_size {params.eval_batch_size} --emb_size {params.emb_size} --num_layers {params.num_layers} --early_stop_n {params.early_stop_n}")
#             else:
#                 shell(f"python {input.script} --data_path {input.data_path}/random_10cv/{index+1} --output_dir {input.output_dir}_fold{index+1} --experiment_name {params.experiment_name} --gpu {params.gpu} --ratio {params.ratio} --cv_num {params.cv_num} --lr {params.lr} --num_workers {params.num_workers} --num_epochs {params.num_epochs} --train_batch_size {params.train_batch_size} --eval_batch_size {params.eval_batch_size} --emb_size {params.emb_size} --num_layers {params.num_layers} --early_stop_n {params.early_stop_n}")

##############################  Other baseline Model Training Pipeline  ##############################