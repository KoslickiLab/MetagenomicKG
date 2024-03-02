"""
This script is used to generate KG-based metagenomic sample embeddings using personalized PageRank.
"""

## Import standard libraries
import os
import sys
from tqdm import tqdm, trange
from glob import glob
import pandas as pd
import networkx as nx
import argparse
import logging
from multiprocessing import Pool

## Import custom libraries
from utils import get_logger, read_tsv_file

def run_personalized_pagerank(sample):
    # Set the personalized vector
    temp_initial_dict = initial_dict.copy()
    temp_initial_dict.update(dict(zip(genome_idx, abund_metagenomic_samples[sample].to_list())))
    
    # Run personalized PageRank
    pr = nx.pagerank(G, alpha=args.alpha, personalization=temp_initial_dict)
    
    # Convert embeddings to dataframe
    embeddings = pd.DataFrame(pr.items(), columns=['node_id', 'embedding'])
    embeddings['node_id'] = embeddings['node_id'].map(idx_to_node)
    embeddings = embeddings.set_index('node_id')
    embeddings.columns = [sample]
    
    return embeddings 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run personalized PageRank on the knowledge graph to generate embeddings for metagenomic samples')
    parser.add_argument('--existing_KG_nodes', type=str, help='path of the existing knowledge graph nodes')
    parser.add_argument('--existing_KG_edges', type=str, help='path of the existing knowledge graph edges')
    parser.add_argument('--abund_metagenomic_samples', type=str, help='path of a CSV file with relative abundance where row is genome and column is sample')
    parser.add_argument('--alpha', type=float, help='Damping parameter for personalized PageRank, default=0.9', default=0.9)
    parser.add_argument('--n_jobs', type=int, help='Number of jobs to run in parallel', default=16)
    parser.add_argument('--output_dir', type=str, help='path of the output directory')
    args = parser.parse_args()

    # Create a logger object
    logger = get_logger()
    logger.setLevel(logging.DEBUG)
    
    # Read KG nodes and edges
    logger.info('Reading KG nodes and edges')
    KG_nodes = read_tsv_file(args.existing_KG_nodes)
    KG_nodes = pd.DataFrame(KG_nodes[1:], columns=KG_nodes[0])
    KG_edges = read_tsv_file(args.existing_KG_edges)
    KG_edges = pd.DataFrame(KG_edges[1:], columns=KG_edges[0])
    
    # Read metagenomic sample file
    logger.info('Reading metagenomic sample file')
    abund_metagenomic_samples = pd.read_csv(args.abund_metagenomic_samples, sep=',', header=0)
    
    # Normalize the relative abundance in each column
    logger.info('Normalizing the relative abundance in each column')
    abund_metagenomic_samples = abund_metagenomic_samples.set_index('name')
    abund_metagenomic_samples = abund_metagenomic_samples.div(abund_metagenomic_samples.sum(axis=0), axis=1)
    
    # generate node mappping (map node to index and index to node, as well as synonyms to node)
    logger.info('Generating node mappings')
    node_to_idx = {node: idx for idx, node in enumerate(KG_nodes['node_id'])}
    idx_to_node = {idx: node for idx, node in enumerate(KG_nodes['node_id'])}
    synonym_to_node = {synonym:node for node, synonym_list in KG_nodes[['node_id', 'synonyms']].to_numpy() for synonym in eval(synonym_list)}
    
    # generate edge dataframe and remove duplicates
    logger.info('Generating edge dataframe and removing duplicates')
    edges_df = pd.DataFrame([(node_to_idx[row[0]], node_to_idx[row[1]]) for row in KG_edges[['source_node','target_node']].to_numpy()])
    edges_df.columns = ['source', 'target']
    edges_df_no_duplicates = edges_df.drop_duplicates().reset_index(drop=True)
    
    # Initialize a directed graph
    G = nx.DiGraph()

    # Add nodes to the graph
    logger.info('Adding nodes to the graph')
    for node in tqdm(node_to_idx.values()):
        G.add_node(node)

    # Add edges to the graph
    logger.info('Adding edges to the graph')
    G.add_edges_from([(row[0], row[1]) for row in edges_df_no_duplicates.to_numpy()])
    
    # Initiating the personalized vector
    logger.info('Initiating the personalized vector')
    initial_dict = {idx: 0 for idx in idx_to_node.keys()}
    
    # Run personalized PageRank for each metagenomic sample
    logger.info('Running personalized PageRank for each metagenomic sample')
    # map genome to idex
    genome_idx = [node_to_idx[synonym_to_node[f'GTDB:{x}']] for x in abund_metagenomic_samples.index] 
    
    # List of sample names to iterate over
    samples = list(abund_metagenomic_samples.columns)
    
    # Execute the function in parallel and collect the results
    logger.info('Executing the function in parallel and collecting the results')
    with Pool(processes=args.n_jobs) as pool:
        results = list(tqdm(pool.imap(run_personalized_pagerank, samples), total=len(samples)))
    
    # Combine the results into a single dataframe
    combined_embeddings = pd.concat(results, axis=1)
    
    # Save the embeddings
    logger.info('Saving the embeddings')
    combined_embeddings.to_csv(os.path.join(args.output_dir, 'sample_embeddings.csv'), sep=',', index=True, header=True)