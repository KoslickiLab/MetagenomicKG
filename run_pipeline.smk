"""
This file is a SnakeMake Script to automate the KEGG-based Pathogen KG construction pipeline.

Usage:
    snakemake --cores 16 -s run_pipeline.smk targets
"""

## Import Python standard libraries
import os, sys
import subprocess

## Define Some Global Variables
CURRENT_PATH = os.getcwd()
GTDB_DOWNLOAD_URL = 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/'
umls_apikey = os.environ['UMLS_API_KEY']
node_synonymizer_dbname = 'node_synonymizer_v1.1_KG2.8.0.1.sqlite'
neo4j_version = '3.5.26'
neo4j_dbname = 'MetagenomicsKG'

## Create Required Folders
if not os.path.exists(os.path.join(CURRENT_PATH, "KEGG_data")):
    os.makedirs(os.path.join(CURRENT_PATH, "KEGG_data"))
if not os.path.exists(os.path.join(CURRENT_PATH, "GTDB_data")):
    os.makedirs(os.path.join(CURRENT_PATH, "GTDB_data"))
if not os.path.exists(os.path.join(CURRENT_PATH, "Micobial_hierarchy")):
    os.makedirs(os.path.join(CURRENT_PATH, "Micobial_hierarchy"))
if not os.path.exists(os.path.join(CURRENT_PATH, "KG")):
    os.makedirs(os.path.join(CURRENT_PATH, "KG"))
if not os.path.exists(os.path.join(CURRENT_PATH, "NodeSynonymizer")):
    os.makedirs(os.path.join(CURRENT_PATH, "NodeSynonymizer"))
if not os.path.exists(os.path.join(CURRENT_PATH, "neo4j")):
    os.makedirs(os.path.join(CURRENT_PATH, "neo4j"))

## Download GTDB Data
# Call the wget command using subprocess
if not os.path.exists(os.path.join(CURRENT_PATH, "GTDB_data", "bac120_taxonomy_r207.tsv")) or \
    not os.path.exists(os.path.join(CURRENT_PATH, "GTDB_data", "bac120_metadata_r207.tar.gz")) or \
    not os.path.exists(os.path.join(CURRENT_PATH, "GTDB_data", "ar53_taxonomy_r207.tsv")) or \
    not os.path.exists(os.path.join(CURRENT_PATH, "GTDB_data", "ar53_metadata_r207.tar.gz")):
    for cur_file in ['bac120_taxonomy_r207.tsv', 'bac120_metadata_r207.tar.gz', 'ar53_taxonomy_r207.tsv', 'ar53_metadata_r207.tar.gz',]:
        destination_path = os.path.join(CURRENT_PATH, "GTDB_data", cur_file)
        result = subprocess.run(["wget", "-O", destination_path, GTDB_DOWNLOAD_URL + cur_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # Check if the download was successful
        if result.returncode == 0:
            if cur_file.endswith('.tar.gz'):
                result = subprocess.run(["tar", "-xvzf", destination_path, '-C', os.path.join(CURRENT_PATH, "GTDB_data")], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if result.returncode == 0:
                    print("File downloaded and extracted successfully!")
                else:
                    print("Error extracting file.")
            else:
                print("File downloaded successfully!")
            print("File downloaded successfully!")
        else:
            print("Error downloading file.")

# Download Neo4j
if not os.path.exists(os.path.join(CURRENT_PATH, "neo4j", f"neo4j-community-{neo4j_version}-unix.tar.gz")):
    result = subprocess.run(["wget", "-O", os.path.join(CURRENT_PATH, "neo4j", f"neo4j-community-{neo4j_version}-unix.tar.gz"), "https://dist.neo4j.org/neo4j-community-{neo4j_version}-unix.tar.gz"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Check if the download was successful
    if result.returncode == 0:
        result = subprocess.run(["tar", "-xvzf", os.path.join(CURRENT_PATH, "neo4j", f"neo4j-community-{neo4j_version}-unix.tar.gz"), '-C', CURRENT_PATH], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode == 0:
            print("File downloaded and extracted successfully!")
        else:
            print("Error extracting file.")
    else:
        print("Error downloading file.")


## Build Rules
rule targets:
    input:
        os.path.join(CURRENT_PATH, "Micobial_hierarchy", "archaea_hierarchy.tsv"),
        os.path.join(CURRENT_PATH, "Micobial_hierarchy", "bacteria_hierarchy.tsv"),  
        os.path.join(CURRENT_PATH, "Micobial_hierarchy", "fungi_hierarchy.tsv"),
        os.path.join(CURRENT_PATH, "Micobial_hierarchy", "viruses_hierarchy.tsv"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_compounds.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_pathways.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_modules.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_glycans.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_reactions.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_enzymes.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_networks.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_koids.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_diseases.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_drugs.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_rclasses.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_dgroups.txt"),
        os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v1.tsv'),
        os.path.join(CURRENT_PATH, "KG", 'KG_edges_v1.tsv'),
        os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v2.tsv'),
        os.path.join(CURRENT_PATH, "KG", 'KG_edges_v2.tsv'),
        os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v3.tsv'),
        os.path.join(CURRENT_PATH, "KG", 'KG_edges_v3.tsv'),
        os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v4.tsv'),
        os.path.join(CURRENT_PATH, "KG", 'KG_edges_v4.tsv'),
        os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v5.tsv'),
        os.path.join(CURRENT_PATH, "KG", 'KG_edges_v5.tsv'),
        os.path.join(CURRENT_PATH, "neo4j", "input_files", 'nodes_header.tsv'),
        os.path.join(CURRENT_PATH, "neo4j", "input_files", 'edges_header.tsv'),
        os.path.join(CURRENT_PATH, "neo4j", "input_files", 'nodes.tsv'),
        os.path.join(CURRENT_PATH, "neo4j", "input_files", 'edges.tsv'),
        os.path.join(CURRENT_PATH, "KG", "import_to_neo4j_done.txt")


# Get Taxonomy hierarchy for Archeaa and Bacteria from GTDB as well as Fungi and Viruses from NCBI Taxonomy
rule step0_get_taxonomy:
    input:
        script = ancient(os.path.join(CURRENT_PATH, "scripts", "get_hierarchy.py")),
        bacteria_taxonomy = ancient(os.path.join(CURRENT_PATH, "GTDB_data", "bac120_taxonomy_r207.tsv")),
        archaea_taxonomy = ancient(os.path.join(CURRENT_PATH,  "GTDB_data", "ar53_taxonomy_r207.tsv")),
        bacteria_metadata = ancient(os.path.join(CURRENT_PATH, "GTDB_data", "bac120_metadata_r207.tsv")),
        archaea_metadata = ancient(os.path.join(CURRENT_PATH, "GTDB_data", "ar53_metadata_r207.tsv")),
        output_dir = ancient(os.path.join(CURRENT_PATH, "Micobial_hierarchy"))
    output:
        os.path.join(CURRENT_PATH, "Micobial_hierarchy", "archaea_hierarchy.tsv"),
        os.path.join(CURRENT_PATH, "Micobial_hierarchy", "bacteria_hierarchy.tsv"),        
        os.path.join(CURRENT_PATH, "Micobial_hierarchy", "fungi_hierarchy.tsv"),
        os.path.join(CURRENT_PATH, "Micobial_hierarchy", "viruses_hierarchy.tsv")
    run:
        shell("python {input.script} --output_dir {input.output_dir} --bacteria_taxonomy {input.bacteria_taxonomy} --archaea_taxonomy {input.archaea_taxonomy} --bacteria_metadata {input.bacteria_metadata} --archaea_metadata {input.archaea_metadata}")

# Process the KEGG FTP data and download KEGG link info from APIs
rule step1_process_kegg_data:
    input:
        script = ancient(os.path.join(CURRENT_PATH, "scripts", "extract_KEGG_data.py")),
        kegg_data_dir = ancient("/data/shared_data/KEGG_FTP/kegg"),
        output_dir = ancient(os.path.join(CURRENT_PATH, "KEGG"))
    params:
        microb_only = True
    output:
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_compounds.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_pathways.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_modules.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_glycans.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_reactions.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_enzymes.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_networks.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_koids.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_diseases.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_drugs.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_rclasses.txt"),
        os.path.join(CURRENT_PATH, "KEGG_data", "kegg_dgroups.txt")
    run:
        if params.microb_only:
            shell("python {input.script} --kegg_data_dir {input.kegg_data_dir} --microb_only --output_dir {input.output_dir}")
        else:
            shell("python {input.script} --kegg_data_dir {input.kegg_data_dir} --output_dir {input.output_dir}")

# Integrate microbial hierarchy into a KG
rule step2_integrate_microbial_hierarchy:
    input:
        script = ancient(os.path.join(CURRENT_PATH, "scripts", "integrate_microbial_hierarchy.py")),
        data_dir = ancient(os.path.join(CURRENT_PATH, "Micobial_hierarchy")),
        bacteria_metadata = ancient(os.path.join(CURRENT_PATH, "GTDB_data", "bac120_metadata_r207.tsv")),
        archaea_metadata = ancient(os.path.join(CURRENT_PATH, "GTDB_data", "ar53_metadata_r207.tsv")),
        archaea_hierarchy_unused = ancient(os.path.join(CURRENT_PATH, "Micobial_hierarchy", "archaea_hierarchy.tsv")),
        output_dir = ancient(os.path.join(CURRENT_PATH, "KG"))
    output:
        os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v1.tsv'),
        os.path.join(CURRENT_PATH, "KG", 'KG_edges_v1.tsv')
    run:
        shell("python {input.script} --data_dir {input.data_dir} --bacteria_metadata {input.bacteria_metadata} --archaea_metadata {input.archaea_metadata} --output_dir {input.output_dir}")


# Integrate knowledge from KEGG into a KG
rule step3_integrate_kegg_data:
    input:
        script = ancient(os.path.join(CURRENT_PATH, "scripts", "integrate_KEGG.py")),
        kegg_data_dir = ancient("/data/shared_data/KEGG_FTP/kegg"),
        kegg_processed_data_dir = ancient(os.path.join(CURRENT_PATH, "KEGG_data")),
        gtdb_assignment = ancient(os.path.join(CURRENT_PATH, "taxonomy_assignment_by_GTDB_tk", "merged_KEGG_assignment.tsv")),
        existing_KG_nodes = ancient(os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v1.tsv')),
        existing_KG_edges = ancient(os.path.join(CURRENT_PATH, "KG", 'KG_edges_v1.tsv')),
        output_dir = ancient(os.path.join(CURRENT_PATH, "KG")),
        unused_file = ancient(os.path.join(CURRENT_PATH, "KEGG_data", "kegg_compounds.txt"))
    params:
        microb_only = True,
        ani_threshold = 0.0,
        af_threshold = 0.0
    output:
        os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v2.tsv'),
        os.path.join(CURRENT_PATH, "KG", 'KG_edges_v2.tsv')
    run:
        if params.microb_only:
            shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --kegg_data_dir {input.kegg_data_dir} --kegg_processed_data_dir {input.kegg_processed_data_dir} --gtdb_assignment {input.gtdb_assignment} --microb_only --ANI_threshold {params.ani_threshold} --AF_threshold {params.af_threshold} --output_dir {input.output_dir}")
        else:
            shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --kegg_data_dir {input.kegg_data_dir} --kegg_processed_data_dir {input.kegg_processed_data_dir} --gtdb_assignment {input.gtdb_assignment} --ANI_threshold {params.ani_threshold} --AF_threshold {params.af_threshold} --output_dir {input.output_dir}")

# Integrate KG2 data into into a KG
rule step4_integrate_kg2_data:
    input:
        script = ancient(os.path.join(CURRENT_PATH, "scripts", "integrate_KG2.py")),
        existing_KG_nodes = ancient(os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v2.tsv')),
        existing_KG_edges = ancient(os.path.join(CURRENT_PATH, "KG", 'KG_edges_v2.tsv')),
        data_dir = ancient(os.path.join(CURRENT_PATH, "KG2")),
        output_dir = ancient(os.path.join(CURRENT_PATH, "KG"))
    output:
        os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v3.tsv'),
        os.path.join(CURRENT_PATH, "KG", 'KG_edges_v3.tsv')
    run:
        shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --data_dir {input.data_dir} --output_dir {input.output_dir}")

# Integrate BVBRC data into the a KG
rule step5_integrate_bvbrc_data:
    input:
        script = ancient(os.path.join(CURRENT_PATH, "scripts", "integrate_BVBRC.py")),
        existing_KG_nodes = ancient(os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v3.tsv')),
        existing_KG_edges = ancient(os.path.join(CURRENT_PATH, "KG", 'KG_edges_v3.tsv')),
        data_dir = ancient(os.path.join(CURRENT_PATH, "pathogen_database", "BV-BRC", "raw_data")),
        gtdb_assignment = ancient(os.path.join(CURRENT_PATH, "taxonomy_assignment_by_GTDB_tk", "merged_BVBRC_assignment.tsv")),
        synonymizer_dir = ancient(os.path.join(CURRENT_PATH, "NodeSynonymizer")),
        output_dir = ancient(os.path.join(CURRENT_PATH, "KG"))
    params:
        umls_api_key = umls_apikey,
        synonymizer_dbname = node_synonymizer_dbname,
        ani_threshold = 0.0,
        af_threshold = 0.0
    output:
        os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v4.tsv'),
        os.path.join(CURRENT_PATH, "KG", 'KG_edges_v4.tsv')
    run:
        shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --data_dir {input.data_dir} --gtdb_assignment {input.gtdb_assignment} --synonymizer_dir {input.synonymizer_dir} --synonymizer_dbname {params.synonymizer_dbname} --umls_api_key {params.umls_api_key} --ANI_threshold {params.ani_threshold} --AF_threshold {params.af_threshold} --output_dir {input.output_dir}")


# Integarte MicroPhenoDB data into a KG
rule step5_integrate_micropheno_data:
    input:
        script = ancient(os.path.join(CURRENT_PATH, "scripts", "integrate_MicroPhenoDB.py")),
        existing_KG_nodes = ancient(os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v4.tsv')),
        existing_KG_edges = ancient(os.path.join(CURRENT_PATH, "KG", 'KG_edges_v4.tsv')),
        data_dir = ancient(os.path.join(CURRENT_PATH, 'pathogen_database', 'MicroPhenoDB')),
        synonymizer_dir = ancient(os.path.join(CURRENT_PATH, "NodeSynonymizer")),
        output_dir = ancient(os.path.join(CURRENT_PATH, "KG"))
    params:
        umls_api_key = umls_apikey,
        synonymizer_dbname = node_synonymizer_dbname
    output:
        os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v5.tsv'),
        os.path.join(CURRENT_PATH, "KG", 'KG_edges_v5.tsv')
    run:
        shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --data_dir {input.data_dir} --synonymizer_dir {input.synonymizer_dir} --synonymizer_dbname {params.synonymizer_dbname} --umls_api_key {params.umls_api_key} --output_dir {input.output_dir}")

# Prepare neo4j input files
rule step6_prepare_neo4j_inputs:
    input:
        script = ancient(os.path.join(CURRENT_PATH, "scripts", "prepare_neo4j_inputs.py")),
        existing_KG_nodes = ancient(os.path.join(CURRENT_PATH, "KG", 'KG_nodes_v5.tsv')),
        existing_KG_edges = ancient(os.path.join(CURRENT_PATH, "KG", 'KG_edges_v5.tsv')),
        kg_dir = ancient(os.path.join(CURRENT_PATH, "KG"))
    params:
        output_dir = os.path.join(CURRENT_PATH, "neo4j", "input_files")
    output:
        os.path.join(CURRENT_PATH, "neo4j", "input_files", 'nodes_header.tsv'),
        os.path.join(CURRENT_PATH, "neo4j", "input_files", 'edges_header.tsv'),
        os.path.join(CURRENT_PATH, "neo4j", "input_files", 'nodes.tsv'),
        os.path.join(CURRENT_PATH, "neo4j", "input_files", 'edges.tsv')
    run:
        shell("python {input.script} --existing_KG_nodes {input.existing_KG_nodes} --existing_KG_edges {input.existing_KG_edges} --kg_dir {input.kg_dir} --output_dir {params.output_dir}")

# Import KG into Neo4j
rule step7_import_kg_into_neo4j:
    input:
        script = ancient(os.path.join(CURRENT_PATH, "scripts", "tsv_to_neo4j.sh")),
        nodes_header = ancient(os.path.join(CURRENT_PATH, "neo4j", "input_files", 'nodes_header.tsv')),
        edges_header = ancient(os.path.join(CURRENT_PATH, "neo4j", "input_files", 'edges_header.tsv')),
        nodes = ancient(os.path.join(CURRENT_PATH, "neo4j", "input_files", 'nodes.tsv')),
        edges = ancient(os.path.join(CURRENT_PATH, "neo4j", "input_files", 'edges.tsv'))
    params:
        tsv_dir = os.path.join(CURRENT_PATH, "neo4j", "input_files"),
        dbname = neo4j_dbname.lower(),
        neo4j_config = os.path.join(CURRENT_PATH, "neo4j", f"neo4j-community-{neo4j_version}", "conf", "neo4j.conf")
    output:
        touch(os.path.join(CURRENT_PATH, "KG", "import_to_neo4j_done.txt"))
    run:
        shell("bash {input.script} {params.dbname} {params.neo4j_config} {params.tsv_dir}")
