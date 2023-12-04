"""
This script is used to integrate all microbial hierarchical data into a knolwedge graph.
"""

## Import standard libraries
import os
import sys
from tqdm import tqdm, trange
from glob import glob
import pandas as pd
import argparse
import logging

## Import custom libraries
from utils import get_logger, check_files, read_tsv_file, Node, Edge, KnowledgeGraph

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Integrate all microbial hierarchical data into a knolwedge graph')
    parser.add_argument('--data_dir', type=str, help='path of the Microbial hierarchy data directory')
    parser.add_argument('--bacteria_metadata', type=str, help='path of the GTDB metadata file for archaea (e.g., bac120_metadata_r207.tsv)')
    parser.add_argument('--archaea_metadata', type=str, help='path of the GTDB metadata file for archaea (e.g., ar53_metadata_r207.tsv)')
    parser.add_argument('--output_dir', type=str, help='path of the output directory')
    args = parser.parse_args()

    # Create a logger object
    logger = get_logger()
    logger.setLevel(logging.INFO)
    
    name_mapping = {}
    # Read bacteria metadata
    logger.info("Reading bacteria metadata...")
    if check_files(args.bacteria_metadata, logger):
        temp_data = read_tsv_file(args.bacteria_metadata)
        bacteria_meta_df = pd.DataFrame(temp_data[1:], columns=temp_data[0])
        name_mapping.update({row[0]:row[1] for row in bacteria_meta_df[['accession', 'ncbi_organism_name']].to_numpy()})
    else:
        sys.exit(1)
        
    # Read archaea metadata
    logger.info("Reading archaea metadata...")
    if check_files(args.archaea_metadata, logger):
        temp_data = read_tsv_file(args.archaea_metadata)
        archaea_meta_df = pd.DataFrame(temp_data[1:], columns=temp_data[0])
        name_mapping.update({row[0]:row[1] for row in archaea_meta_df[['accession', 'ncbi_organism_name']].to_numpy()})
    else:
        sys.exit(1)
    
    # Create a knowledge graph object
    logger.info("Creating a knowledge graph object...")
    kg = KnowledgeGraph(logger)
    
    ## Read all microbial hierarchy data
    logger.info("Reading all microbial hierarchy data...")
    microbial_hierarchy_files = os.listdir(args.data_dir)
    for microbial_hierarchy_file in microbial_hierarchy_files:
        microbial_hierarchy_file_path = os.path.join(args.data_dir, microbial_hierarchy_file)
        logger.info("Reading {}...".format(microbial_hierarchy_file_path))
        temp_data = read_tsv_file(microbial_hierarchy_file_path)
        temp_df = pd.DataFrame(temp_data[1:], columns=temp_data[0])
        for row in tqdm(temp_df.to_numpy(), desc="Integrating {}...".format(microbial_hierarchy_file), total=len(temp_df)):
            if microbial_hierarchy_file.split('_')[0] in ['viruses', 'fungi']:
                kg_source = 'NCBI'
                source_synonyms = [f"NCBI:{row[0]}"]
                source_name = [row[0]]
                if row[2] == '':
                    source_link = [f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name={row[0]}"]
                else:
                    source_link = [f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={row[2]}"]
                target_synonyms = [f"NCBI:{row[3]}"]
                target_name = [row[3]]
                if row[5] == '':
                    target_link = [f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name={row[3]}"]
                else:
                    target_link = [f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={row[5]}"]
            elif microbial_hierarchy_file.split('_')[0] in ['bacteria', 'archaea']:
                kg_source = 'GTDB'
                if row[0].split('_')[0] in ['GB', 'RS']:
                    source_name = [name_mapping[row[0]]]
                    row[0] = row[0].replace('GB_', '').replace('RS_', '')
                    source_link = [f"https://gtdb.ecogenomic.org/genome?gid={row[0]}"]
                    source_link += [f"https://www.ncbi.nlm.nih.gov/assembly/{row[0]}"]
                    if row[2] != '':
                        source_link += [f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={row[2]}"]
                else:
                    source_name = [row[0]]
                    if row[2] != '':
                        source_link = [f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={row[2]}"]
                    else:
                        source_link = []
                source_synonyms = [f"GTDB:{row[0]}"]
                
                if row[3].split('_')[0] in ['GB', 'RS']:
                    target_name = [name_mapping[row[3]]]
                    row[3] = row[3].replace('GB_', '').replace('RS_', '')
                    target_link = [f"https://gtdb.ecogenomic.org/genome?gid={row[3]}"]
                    target_link += [f"https://www.ncbi.nlm.nih.gov/assembly/{row[3]}"]
                    if row[5] != '':
                        target_link += [f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={row[5]}"]
                else:
                    target_name = [row[3]]
                    if row[5] != '':
                        target_link = [f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={row[5]}"]
                    else:
                        target_link = []
                target_synonyms = [f"GTDB:{row[3]}"]

            source_node = Node(node_type="Microbe", all_names=source_name, description=[('rank',row[1]), ('taxid', row[2])], knowledge_source=[kg_source], synonyms=source_synonyms, link=source_link)
            target_node = Node(node_type="Microbe", all_names=target_name, description=[('rank',row[4]), ('taxid', row[5])], knowledge_source=[kg_source], synonyms=target_synonyms, link=target_link)

            # Add nodes to the knowledge graph
            kg.add_node(source_node)
            kg.add_node(target_node)
            # Add edge to the knowledge graph
            kg.add_edge(Edge(source_node=source_synonyms[0], target_node=target_synonyms[0], predicate='biolink:has_part', knowledge_source=[kg_source]))
            kg.add_edge(Edge(source_node=target_synonyms[0], target_node=source_synonyms[0], predicate='biolink:part_of', knowledge_source=[kg_source]))


    # Save the knowledge graph
    logger.info("Saving the knowledge graph...")
    kg.save_graph(save_dir = args.output_dir, node_filename = 'KG_nodes_v1.tsv', edge_filename = 'KG_edges_v1.tsv')
    logger.info("KG node is saved to {}".format(os.path.join(args.output_dir, 'KG_nodes_v1.tsv')))
    logger.info("KG edge is saved to {}".format(os.path.join(args.output_dir, 'KG_edges_v1.tsv')))
    
