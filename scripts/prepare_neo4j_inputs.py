"""
This script is used to prepare the input files for Neo4j.
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
from utils import get_logger, KnowledgeGraph


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Integrate all MicroPhenoDB into a knolwedge graph')
    parser.add_argument('--existing_KG_nodes', type=str, help='path of the existing knowledge graph nodes')
    parser.add_argument('--existing_KG_edges', type=str, help='path of the existing knowledge graph edges')
    parser.add_argument('--kg_dir', type=str, help='path of the knowledge graph directory')
    parser.add_argument('--output_dir', type=str, help='path of the output directory')
    args = parser.parse_args()
    
    # Create a logger object
    logger = get_logger()
    logger.setLevel(logging.INFO)

    # Create a knowledge graph object
    logger.info("Creating a knowledge graph object...")
    kg = KnowledgeGraph(logger)
    # Load existing knowledge graph nodes and edges
    logger.info("Loading existing knowledge graph nodes and edges...")
    node_filename = args.existing_KG_nodes.split('/')[-1]
    edge_filename = args.existing_KG_edges.split('/')[-1]
    kg.load_graph(load_dir=args.kg_dir, node_filename=node_filename, edge_filename=edge_filename)

    # Convert knowledge graph nodes and edges to dataframes
    logger.info("Converting knowledge graph nodes and edges to dataframes...")
    kg_node_df = kg.nodes_to_dataframe()
    kg_edge_df = kg.edges_to_dataframe()
    
    ## Construct the headers
    logger.info("Constructing the headers...")
    node_headers = list(kg_node_df.columns)
    node_headers[node_headers.index('node_id')] = node_headers[node_headers.index('node_id')] + ':ID'
    node_headers[node_headers.index('all_names')] = node_headers[node_headers.index('all_names')] + ':string[]'
    node_headers[node_headers.index('knowledge_source')] = node_headers[node_headers.index('knowledge_source')] + ':string[]'
    node_headers[node_headers.index('link')] = node_headers[node_headers.index('link')] + ':string[]'
    node_headers[node_headers.index('synonyms')] = node_headers[node_headers.index('synonyms')] + ':string[]'
    node_headers[node_headers.index('is_pathogen')] = node_headers[node_headers.index('is_pathogen')]
    node_headers += [':LABEL']
    
    edge_headers = list(kg_edge_df.columns)
    edge_headers[edge_headers.index('source_node')] = edge_headers[edge_headers.index('source_node')]
    edge_headers[edge_headers.index('target_node')] = edge_headers[edge_headers.index('target_node')]
    edge_headers[edge_headers.index('knowledge_source')] = edge_headers[edge_headers.index('knowledge_source')] + ':string[]'
    edge_headers += [':TYPE', ':START_ID', ':END_ID']

    ## Write headers to a TSV file
    logger.info("Writing headers to a TSV file...")
    with open(os.path.join(args.output_dir, 'nodes_header.tsv'), 'w') as f:
        f.write('\t'.join(node_headers)+'\n')

    with open(os.path.join(args.output_dir, 'edges_header.tsv'), 'w') as f:
        f.write('\t'.join(edge_headers)+'\n')

    # Format nodes for neo4j
    logger.info("Formatting nodes for neo4j...")
    kg_node_temp = []
    for row in tqdm(kg_node_df.to_numpy()):
        node_id, node_type, all_names, description, knowledge_source, link, synonyms, is_pathogen = row
        ## remove KO related genes for better neo4j visualization
        description = [x for x in description if x[0] != 'KO_related_genes']
        all_names = 'ǂ'.join(all_names)
        description = '; '.join([pair[0]+":"+pair[1] for pair in description if pair[1]])
        knowledge_source = 'ǂ'.join(knowledge_source)
        link = 'ǂ'.join(link)
        synonyms = 'ǂ'.join(synonyms)
        kg_node_temp += [[node_id, node_type, all_names, description, knowledge_source, link, synonyms, is_pathogen, node_type]]
    kg_node_df = pd.DataFrame(kg_node_temp)
    
    ## Save KG nodes with neo4j format
    logger.info("Saving KG nodes with neo4j format...")
    kg_node_df.to_csv(os.path.join(args.output_dir, 'nodes.tsv'), sep='\t', index=None, header=False)
    
    # Format edges for neo4j
    logger.info("Formatting edges for neo4j...")
    kg_edge_temp = []
    for row in tqdm(kg_edge_df.to_numpy()):
        source_node, target_node, predicate, knowledge_source = row
        knowledge_source = 'ǂ'.join(knowledge_source)
        kg_edge_temp += [[source_node, target_node, predicate, knowledge_source, predicate, source_node, target_node]]
    kg_edge_df = pd.DataFrame(kg_edge_temp)
    
    ## Save KG edges with neo4j format
    logger.info("Saving KG edges with neo4j format...")
    kg_edge_df.to_csv(os.path.join(args.output_dir, 'edges.tsv'), sep='\t', index=None, header=False)
