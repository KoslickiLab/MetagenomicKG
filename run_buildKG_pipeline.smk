"""
This file is a SnakeMake Script to automate the MetagenomicKG construction pipeline.

Usage:
    snakemake --cores 16 -s run_buildKG_pipeline.smk targets
"""
## Import Config Files
configfile: "./config.yaml"

## Import Python standard libraries
import os, sys
import subprocess

## Define Some Global Variables
ROOT_PATH = os.getcwd()
DATA_PATH = os.path.join(ROOT_PATH, 'data')
SCRIPT_PATH = os.path.join(ROOT_PATH, 'build_KG')
GTDB_version = config['BUILD_KG_VARIABLES']['GTDB_VERSION']
GTDB_DOWNLOAD_URL = f'https://data.gtdb.ecogenomic.org/releases/release{GTDB_version}/{GTDB_version}.0/'
if 'UMLS_API_KEY' in os.environ and os.environ['UMLS_API_KEY']:
    umls_apikey = os.environ['UMLS_API_KEY']
elif 'UMLS_API_KEY' in config['BUILD_KG_VARIABLES'] and config['BUILD_KG_VARIABLES']['UMLS_API_KEY']:
    umls_apikey = config['BUILD_KG_VARIABLES']['UMLS_API_KEY']
else:
    raise ValueError("UMLS_API_KEY is not set in the environment or the config file. Please set it in 'config.yaml' file before running the pipeline.")
node_synonymizer_dbname = 'node_synonymizer_v1.0_KG2.8.4.sqlite'
neo4j_dbname = 'MetagenomicsKG'

## Create Required Folders
if not os.path.exists(os.path.join(DATA_PATH, "KEGG_data")):
    os.makedirs(os.path.join(DATA_PATH, "KEGG_data"))
if not os.path.exists(os.path.join(DATA_PATH, "GTDB_data")):
    os.makedirs(os.path.join(DATA_PATH, "GTDB_data"))
if not os.path.exists(os.path.join(DATA_PATH, "Micobial_hierarchy")):
    os.makedirs(os.path.join(DATA_PATH, "Micobial_hierarchy"))
if not os.path.exists(os.path.join(DATA_PATH, "merged_KG")):
    os.makedirs(os.path.join(DATA_PATH, "merged_KG"))
if not os.path.exists(os.path.join(DATA_PATH, "neo4j")):
    os.makedirs(os.path.join(DATA_PATH, "neo4j"))

## Download GTDB Data
# Call the wget command using subprocess
if not os.path.exists(os.path.join(DATA_PATH, "GTDB_data", f"bac120_taxonomy_r{GTDB_version}.tsv")) or \
    not os.path.exists(os.path.join(DATA_PATH, "GTDB_data", f"bac120_metadata_r{GTDB_version}.tar.gz")) or \
    not os.path.exists(os.path.join(DATA_PATH, "GTDB_data", f"ar53_taxonomy_r{GTDB_version}.tsv")) or \
    not os.path.exists(os.path.join(DATA_PATH, "GTDB_data", f"ar53_metadata_r{GTDB_version}.tar.gz")):
    for cur_file in [f'bac120_taxonomy_r{GTDB_version}.tsv', f'bac120_metadata_r{GTDB_version}.tar.gz', f'ar53_taxonomy_r{GTDB_version}.tsv', f'ar53_metadata_r{GTDB_version}.tar.gz',]:
        destination_path = os.path.join(DATA_PATH, "GTDB_data", cur_file)
        result = subprocess.run(["wget", "-O", destination_path, GTDB_DOWNLOAD_URL + cur_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # Check if the download was successful
        if result.returncode == 0:
            if cur_file.endswith('.tar.gz'):
                result = subprocess.run(["tar", "-xvzf", destination_path, '-C', os.path.join(DATA_PATH, "GTDB_data")], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if result.returncode == 0:
                    print("File downloaded and extracted successfully!")
                else:
                    print("Error extracting file.")
            else:
                print("File downloaded successfully!")
            print("File downloaded successfully!")
        else:
            print("Error downloading file.")

## Download RTX-KG2 data
if not os.path.exists(os.path.join(DATA_PATH, "RTX_KG2", "kg2c-tsv.tar.gz")):
    raise ValueError("the kg2c-tsv.tar.gz file is not present in the RTX_KG2 folder. Please contact authors to get access and download it to ./data/RTX_KG2 folder.")
else:
    if not os.path.exists(os.path.join(DATA_PATH, "RTX_KG2", "nodes_c_header.tsv")) or \
        not os.path.exists(os.path.join(DATA_PATH, "RTX_KG2", "edges_c_header.tsv")) or \
        not os.path.exists(os.path.join(DATA_PATH, "RTX_KG2", "nodes_c.tsv")) or \
        not os.path.exists(os.path.join(DATA_PATH, "RTX_KG2", "edges_c.tsv")):
        result = subprocess.run(["tar", "-xvzf", os.path.join(DATA_PATH, "RTX_KG2", "kg2c-tsv.tar.gz"), '-C', os.path.join(DATA_PATH, "RTX_KG2")], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode == 0:
            print("File extracted successfully!")
        else:
            print("Error extracting file.")

## Check required files from Zenodo and download if not present
# Download NodeSynonymizer database
if not os.path.exists(os.path.join(DATA_PATH, "Zenodo_data", node_synonymizer_dbname)):
    # add code to download the file from Zenodo
    pass

# Download pathogen database tar file
if not os.path.exists(os.path.join(DATA_PATH, "Zenodo_data", "pathogen_database.tar.gz")):
    # add code to download the file from Zenodo
    pass
else:
    if not os.path.exists(os.path.join(DATA_PATH, "pathogen_database")):
        result = subprocess.run(["tar", "-xvzf", os.path.join(DATA_PATH, "Zenodo_data", "pathogen_database.tar.gz"), '-C', os.path.join(DATA_PATH, "Zenodo_data")], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode == 0:
            print("File extracted successfully!")
        else:
            print("Error extracting file for pathogen_database.")

# Download GTDB_tk taxonomy assignment
if not os.path.exists(os.path.join(DATA_PATH, "Zenodo_data", "taxonomy_assignment_by_GTDB_tk.tar.gz")):
    # add code to download the file from Zenodo
    pass
else:
    if not os.path.exists(os.path.join(DATA_PATH, "Zenodo_data", "taxonomy_assignment_by_GTDB_tk")):
        result = subprocess.run(["tar", "-xvzf", os.path.join(DATA_PATH, "Zenodo_data", "taxonomy_assignment_by_GTDB_tk.tar.gz"), '-C', os.path.join(DATA_PATH, "Zenodo_data")], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode == 0:
            print("File extracted successfully!")
        else:
            print("Error extracting file for taxonomy_assignment_by_GTDB_tk.")

# Download AMRFinderResults
if not os.path.exists(os.path.join(DATA_PATH, "Zenodo_data", "AMRFinderResults.tar.gz")):
    # add code to download the file from Zenodo
    pass
else:
    if not os.path.exists(os.path.join(DATA_PATH, "Zenodo_data", "AMRFinderResults")):
        result = subprocess.run(["tar", "-xvzf", os.path.join(DATA_PATH, "Zenodo_data", "AMRFinderResults.tar.gz"), '-C', os.path.join(DATA_PATH, "Zenodo_data")], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode == 0:
            print("File extracted successfully!")
        else:
            print("Error extracting file for AMRFinderResults.")


## Build Rules
rule targets:
    input:
        os.path.join(DATA_PATH, "Micobial_hierarchy", "archaea_hierarchy.tsv"),
        os.path.join(DATA_PATH, "Micobial_hierarchy", "bacteria_hierarchy.tsv"),  
        os.path.join(DATA_PATH, "Micobial_hierarchy", "fungi_hierarchy.tsv"),
        os.path.join(DATA_PATH, "Micobial_hierarchy", "viruses_hierarchy.tsv"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_compounds.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_pathways.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_modules.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_glycans.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_reactions.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_enzymes.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_networks.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_koids.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_diseases.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_drugs.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_rclasses.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_dgroups.txt"),
        os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v1.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v1.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v2.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v2.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v3.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v3.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v4.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v4.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v5.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v5.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v6.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v6.tsv')
        os.path.join(ROOT_PATH, "neo4j", "input_files", 'nodes.tsv'),
        os.path.join(ROOT_PATH, "neo4j", "input_files", 'edges.tsv')


# Get Taxonomy hierarchy for Archeaa and Bacteria from GTDB as well as Fungi and Viruses from NCBI Taxonomy
rule step0_get_taxonomy:
    input:
        script = ancient(os.path.join(SCRIPT_PATH, "get_hierarchy.py")),
        bacteria_taxonomy = ancient(os.path.join(DATA_PATH, "GTDB_data", f"bac120_taxonomy_r{GTDB_version}.tsv")),
        archaea_taxonomy = ancient(os.path.join(DATA_PATH,  "GTDB_data", f"ar53_taxonomy_r{GTDB_version}.tsv")),
        bacteria_metadata = ancient(os.path.join(DATA_PATH, "GTDB_data", f"bac120_metadata_r{GTDB_version}.tsv")),
        archaea_metadata = ancient(os.path.join(DATA_PATH, "GTDB_data", f"ar53_metadata_r{GTDB_version}.tsv")),
        output_dir = ancient(os.path.join(DATA_PATH, "Micobial_hierarchy"))
    output:
        os.path.join(DATA_PATH, "Micobial_hierarchy", "archaea_hierarchy.tsv"),
        os.path.join(DATA_PATH, "Micobial_hierarchy", "bacteria_hierarchy.tsv"),        
        os.path.join(DATA_PATH, "Micobial_hierarchy", "fungi_hierarchy.tsv"),
        os.path.join(DATA_PATH, "Micobial_hierarchy", "viruses_hierarchy.tsv")
    run:
        shell("python {input.script} --output_dir {input.output_dir} --bacteria_taxonomy {input.bacteria_taxonomy} --archaea_taxonomy {input.archaea_taxonomy} --bacteria_metadata {input.bacteria_metadata} --archaea_metadata {input.archaea_metadata}")

# Process the KEGG FTP data and download KEGG link info from APIs
rule step1_process_kegg_data:
    input:
        script = ancient(os.path.join(SCRIPT_PATH, "kegg_utils", "extract_KEGG_data.py")),
        kegg_data_dir = ancient(config['BUILD_KG_VARIABLES']['KEGG_FTP_DATA_DIR']),
        output_dir = ancient(os.path.join(DATA_PATH, "KEGG"))
    params:
        microb_only = True
    output:
        os.path.join(DATA_PATH, "KEGG_data", "kegg_compounds.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_pathways.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_modules.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_glycans.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_reactions.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_enzymes.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_networks.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_koids.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_diseases.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_drugs.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_rclasses.txt"),
        os.path.join(DATA_PATH, "KEGG_data", "kegg_dgroups.txt")
    run:
        if params.microb_only:
            shell("python {input.script} --kegg_data_dir {input.kegg_data_dir} --microb_only --output_dir {input.output_dir}")
        else:
            shell("python {input.script} --kegg_data_dir {input.kegg_data_dir} --output_dir {input.output_dir}")

# Integrate microbial hierarchy into a KG
rule step2_integrate_microbial_hierarchy:
    input:
        script = ancient(os.path.join(SCRIPT_PATH, "integrate_microbial_hierarchy.py")),
        data_dir = ancient(os.path.join(DATA_PATH, "Micobial_hierarchy")),
        bacteria_metadata = ancient(os.path.join(DATA_PATH, "GTDB_data", "bac120_metadata_r214.tsv")),
        archaea_metadata = ancient(os.path.join(DATA_PATH, "GTDB_data", "ar53_metadata_r214.tsv")),
        archaea_hierarchy_unused = ancient(os.path.join(DATA_PATH, "Micobial_hierarchy", "archaea_hierarchy.tsv")),
        output_dir = ancient(os.path.join(DATA_PATH, "merged_KG"))
    output:
        os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v1.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v1.tsv')
    run:
        shell("python {input.script} --data_dir {input.data_dir} --bacteria_metadata {input.bacteria_metadata} --archaea_metadata {input.archaea_metadata} --output_dir {input.output_dir}")


# Integrate knowledge from KEGG into a KG
rule step3_integrate_kegg_data:
    input:
        script = ancient(os.path.join(SCRIPT_PATH, "integrate_KEGG.py")),
        kegg_data_dir = ancient(config['BUILD_KG_VARIABLES']['KEGG_FTP_DATA_DIR']),
        kegg_processed_data_dir = ancient(os.path.join(DATA_PATH, "KEGG_data")),
        gtdb_assignment = ancient(os.path.join(DATA_PATH, "Zenodo_data", "taxonomy_assignment_by_GTDB_tk", "merged_KEGG_assignment.tsv")),
        existing_KG_nodes = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v1.tsv')),
        existing_KG_edges = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v1.tsv')),
        output_dir = ancient(os.path.join(DATA_PATH, "merged_KG")),
        unused_file = ancient(os.path.join(DATA_PATH, "KEGG_data", "kegg_compounds.txt"))
    params:
        microb_only = True,
        ani_threshold = 99.5,
        af_threshold = 0.0
    output:
        os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v2.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v2.tsv')
    run:
        if params.microb_only:
            shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --kegg_data_dir {input.kegg_data_dir} --kegg_processed_data_dir {input.kegg_processed_data_dir} --gtdb_assignment {input.gtdb_assignment} --microb_only --ANI_threshold {params.ani_threshold} --AF_threshold {params.af_threshold} --output_dir {input.output_dir}")
        else:
            shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --kegg_data_dir {input.kegg_data_dir} --kegg_processed_data_dir {input.kegg_processed_data_dir} --gtdb_assignment {input.gtdb_assignment} --ANI_threshold {params.ani_threshold} --AF_threshold {params.af_threshold} --output_dir {input.output_dir}")

# Integrate KG2 data into into a KG
rule step4_integrate_kg2_data:
    input:
        script = ancient(os.path.join(SCRIPT_PATH, "integrate_KG2.py")),
        existing_KG_nodes = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v2.tsv')),
        existing_KG_edges = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v2.tsv')),
        data_dir = ancient(os.path.join(DATA_PATH, "RTX_KG2")),
        output_dir = ancient(os.path.join(DATA_PATH, "merged_KG"))
    output:
        os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v3.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v3.tsv')
    run:
        shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --data_dir {input.data_dir} --output_dir {input.output_dir}")

# Integrate BVBRC data into the a KG
rule step5_integrate_bvbrc_data:
    input:
        script = ancient(os.path.join(SCRIPT_PATH, "integrate_BVBRC.py")),
        existing_KG_nodes = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v3.tsv')),
        existing_KG_edges = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v3.tsv')),
        data_dir = ancient(os.path.join(DATA_PATH, "Zenodo_data", "pathogen_database", "BV-BRC")),
        gtdb_assignment = ancient(os.path.join(DATA_PATH, "Zenodo_data", "taxonomy_assignment_by_GTDB_tk", "merged_BVBRC_assignment.tsv")),
        synonymizer_dir = ancient(os.path.join(DATA_PATH, "Zenodo_data")),
        output_dir = ancient(os.path.join(DATA_PATH, "merged_KG"))
    params:
        umls_api_key = umls_apikey,
        synonymizer_dbname = node_synonymizer_dbname,
        ani_threshold = 99.5,
        af_threshold = 0.0
    output:
        os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v4.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v4.tsv')
    run:
        shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --data_dir {input.data_dir} --gtdb_assignment {input.gtdb_assignment} --synonymizer_dir {input.synonymizer_dir} --synonymizer_dbname {params.synonymizer_dbname} --umls_api_key {params.umls_api_key} --ANI_threshold {params.ani_threshold} --AF_threshold {params.af_threshold} --output_dir {input.output_dir}")


# Integarte MicroPhenoDB data into a KG
rule step5_integrate_micropheno_data:
    input:
        script = ancient(os.path.join(SCRIPT_PATH, "integrate_MicroPhenoDB.py")),
        existing_KG_nodes = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v4.tsv')),
        existing_KG_edges = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v4.tsv')),
        data_dir = ancient(os.path.join(DATA_PATH, "Zenodo_data", 'pathogen_database', 'MicroPhenoDB')),
        synonymizer_dir = ancient(os.path.join(DATA_PATH, "Zenodo_data")),
        output_dir = ancient(os.path.join(DATA_PATH, "merged_KG"))
    params:
        umls_api_key = umls_apikey,
        synonymizer_dbname = node_synonymizer_dbname
    output:
        os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v5.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v5.tsv')
    run:
        shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --data_dir {input.data_dir} --synonymizer_dir {input.synonymizer_dir} --synonymizer_dbname {params.synonymizer_dbname} --umls_api_key {params.umls_api_key} --output_dir {input.output_dir}")

# Integrate AMR data into a KG
rule step6_integrate_amr_data:
    input:
        script = ancient(os.path.join(SCRIPT_PATH, "integrate_AMR.py")),
        existing_KG_nodes = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v5.tsv')),
        existing_KG_edges = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v5.tsv')),
        amr_result = ancient(os.path.join(DATA_PATH, "Zenodo_data", 'AMRFinderResults', 'all_AMR_result.tsv')),
        amr_metadata = ancient(os.path.join(DATA_PATH, "Zenodo_data", 'AMRFinderResults', 'ReferenceGeneCatalog.txt')),
        output_dir = ancient(os.path.join(DATA_PATH, "merged_KG"))
    params:
        coverage_threshold = 80,
        identity_threshold = 90
    output:
        os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v6.tsv'),
        os.path.join(DATA_PATH, "merged_KG", 'KG_edges_v6.tsv')
    run:
        shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --amr_result {input.amr_result} --amr_metadata {input.amr_metadata} --coverage_threshold {params.coverage_threshold} --identity_threshold {params.identity_threshold} --output_dir {input.output_dir}")


# Prepare neo4j input files
rule step7_prepare_neo4j_inputs:
    input:
        script = ancient(os.path.join(SCRIPT_PATH, "neo4j_utils", "prepare_neo4j_inputs.py")),
        existing_KG_nodes = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v6.tsv')),
        existing_KG_edges = ancient(os.path.join(DATA_PATH, "merged_KG", 'KG_nodes_v6.tsv')),
        kg_dir = ancient(os.path.join(DATA_PATH, "merged_KG"))
    params:
        output_dir = os.path.join(ROOT_PATH, "neo4j", "input_files")
    output:
        os.path.join(ROOT_PATH, "neo4j", "input_files", 'nodes.tsv'),
        os.path.join(ROOT_PATH, "neo4j", "input_files", 'edges.tsv')
    run:
        shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --kg_dir {input.kg_dir} --output_dir {params.output_dir}")
