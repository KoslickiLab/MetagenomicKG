## Import standard libraries
import os
import sys
from tqdm import tqdm, trange
from glob import glob
import pandas as pd
import random
import logging
import argparse

## Import custom libraries
from utils import get_logger, read_tsv_file, set_seed

def find_specific_parent(microbe_id, ranks, microbe_info, find_parent):
    """
    Find the specific parent of a microbe
    """
    
    def _recursive_find_parent(microbe_id, ranks, microbe_info, find_parent, res, count=0):
        """
        iterate all the way up to the top of the hierarchy and find all parents belong to the specific ranks
        """
        if count == 8:
            return res
        
        count += 1
        if microbe_id not in find_parent:
            return res
        elif len(microbe_info[microbe_id][0].intersection(set(ranks))) > 0:
            res += [microbe_id]
            return _recursive_find_parent(find_parent[microbe_id], ranks, microbe_info, find_parent, res, count)
        else:
            return _recursive_find_parent(find_parent[microbe_id], ranks, microbe_info, find_parent, res, count)
    
    res = _recursive_find_parent(microbe_id, ranks, microbe_info, find_parent, [], 0)
    if len(res) == 0:
        return [microbe_id]
    else:
        return res + [microbe_id]


def extract_info(sysnonyms, microbe_info_dict, find_parent):
    """
    Get the rank and parent of a microbe based on its id
    """
    if isinstance(sysnonyms, str):
        sysnonyms = eval(sysnonyms)

    rank_list = []
    taxon_id_list = []
    parent_list = []
    for synonym in sysnonyms:
        prefix = synonym.split(':')[0]
        if prefix not in ['NCBI', 'GTDB']:
            continue
        if synonym not in microbe_info_dict:
            continue
        rank_list += list(microbe_info_dict[synonym][0])
        taxon_id_list += list(microbe_info_dict[synonym][1])
        parent_list.append(find_parent[synonym])

    parent_rank_list = []
    parent_taxon_id_list = []
    for parent in parent_list:
        parent_rank_list += list(microbe_info_dict[parent][0])
        parent_taxon_id_list += list(microbe_info_dict[parent][1])
    
    sysnonyms = str(sysnonyms)
    rank_list = str(rank_list)
    taxon_id_list = str(taxon_id_list)
    parent_list = str(parent_list)
    parent_rank_list = str(parent_rank_list)
    parent_taxon_id_list = str(parent_taxon_id_list)
    
    return sysnonyms, rank_list, taxon_id_list, parent_list, parent_rank_list, parent_taxon_id_list


def get_args():
    """
    Parse the arguments from command line
    """
    parser = argparse.ArgumentParser(description='Prepare data for the model training')
    parser.add_argument('--existing_KG_nodes', type=str, help='path of the existing knowledge graph nodes')
    parser.add_argument('--existing_KG_edges', type=str, help='path of the existing knowledge graph edges')
    parser.add_argument('--micobial_hierarchy_dir', type=str, help='path of the microbial hierarchy directory')
    parser.add_argument('--random_seed', type=int, default=100, help='random seed (default: 100)')
    parser.add_argument('--output_dir', type=str, help='path of the output directory')
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    
    # Create a logger object
    logger = get_logger()
    logger.setLevel(logging.DEBUG)
    
    # set the seed for reproducibility
    set_seed(args.random_seed)
    
    # Load existing knowledge graph nodes and edges
    logger.info("Loading existing knowledge graph nodes and edges...")
    KG_nodes_data = read_tsv_file(args.existing_KG_nodes)
    KG_nodes_header = KG_nodes_data[0]
    KG_edges_data = read_tsv_file(args.existing_KG_edges)
    KG_edges_header = KG_edges_data[0]

    # Load microbial hierarchy data
    find_parent = dict()
    find_child = dict()
    microbe_info_dict = dict()
    microbial_hierarchy_files = os.listdir(args.micobial_hierarchy_dir)
    for microbial_hierarchy_file in tqdm(microbial_hierarchy_files, desc='Reading microbial hierarchy files'):
        if microbial_hierarchy_file.split('_')[0] in ['viruses', 'fungi']:
            prefix = 'NCBI:'
        else:
            prefix = 'GTDB:'
        microbial_hierarchy_file_path = os.path.join(args.micobial_hierarchy_dir, microbial_hierarchy_file)
        logger.info("Reading {}...".format(microbial_hierarchy_file_path))
        temp_data = read_tsv_file(microbial_hierarchy_file_path)
        for row in temp_data[1:]:
            row[0] = row[0].replace('GB_', '').replace('RS_', '')
            row[3] = row[3].replace('GB_', '').replace('RS_', '')
            if f"{prefix}{row[0]}" not in microbe_info_dict:
                microbe_info_dict[f"{prefix}{row[0]}"] = [set(), set()]
            microbe_info_dict[f"{prefix}{row[0]}"][0].add(row[1])
            microbe_info_dict[f"{prefix}{row[0]}"][1].add(row[2])              
            if f"{prefix}{row[3]}" not in microbe_info_dict:
                microbe_info_dict[f"{prefix}{row[3]}"] = [set(), set()]
            microbe_info_dict[f"{prefix}{row[3]}"][0].add(row[4])
            microbe_info_dict[f"{prefix}{row[3]}"][1].add(row[5])
            if f"{prefix}{row[3]}" != f"{prefix}{row[0]}":
                find_parent[f"{prefix}{row[3]}"] = f"{prefix}{row[0]}"
                if f"{prefix}{row[0]}" not in find_child:
                    find_child[f"{prefix}{row[0]}"] = set()
                find_child[f"{prefix}{row[0]}"].add(f"{prefix}{row[3]}")

    ## Extract data from MKG for GNN model training 
    logger.info("Extracting MKG for GNN model training...")
    node_info = pd.DataFrame([(row[KG_nodes_header.index('node_id')], row[KG_nodes_header.index('node_type')], \
                                row[KG_nodes_header.index('all_names')], row[KG_nodes_header.index('synonyms')], row[KG_nodes_header.index('is_pathogen')]) \
                                for row in KG_nodes_data[1:]])
    node_info.columns = ['node_id', 'node_type', 'all_names', 'synonyms', 'is_pathogen']
    edge_info = pd.DataFrame([(row[KG_edges_header.index('source_node')], row[KG_edges_header.index('target_node')], row[KG_edges_header.index('predicate')]) for row in KG_edges_data[1:]])
    edge_info.columns = ['source_node', 'target_node', 'predicate']

    # Filter out the microbe-disease edges
    edge_info_temp = []
    filtered_type = ['Microbe', 'Disease']
    for row in tqdm(edge_info.to_numpy(), desc="Filtering out the microbe-disease edges"):
        source_node, target_node, predicate = row
        source_type = source_node.split(':')[0]
        target_type = target_node.split(':')[0]
        if source_type == target_type:
            edge_info_temp.append(row)
        else:
            if source_type in filtered_type and target_type in filtered_type:
                continue
            else:
                edge_info_temp.append(row)
    edge_info = pd.DataFrame(edge_info_temp)
    edge_info.columns = ['source_node', 'target_node', 'predicate']
    
    ## Covert node and edge to index
    # save node to index mapping
    all_nodes = node_info['node_id'].to_list()
    random.shuffle(all_nodes)
    node_to_index = pd.DataFrame(all_nodes).reset_index()[[0, 'index']]
    node_to_index.columns = ['node_id', 'node_index']
    node_to_index.to_csv(os.path.join(args.output_dir, 'node_to_index.tsv'), sep='\t', index=False)
    node_info.to_csv(os.path.join(args.output_dir, 'node_info.tsv'), sep='\t', index=False)
    
    # save pathogen info
    pathogen_info = node_info.loc[node_info['is_pathogen']=='True', :].reset_index(drop=True)
    # find pathogen rank and its parent info
    part1 = pathogen_info[['node_id', 'node_type', 'all_names', 'is_pathogen']]
    part2 = pd.DataFrame(pathogen_info['synonyms'].apply(lambda x: extract_info(x, microbe_info_dict, find_parent)).to_list(), columns=['synonyms', 'rank_list', 'taxon_id_list', 'parent_list', 'parent_rank_list', 'parent_taxon_id_list'])
    pathogen_info = pd.concat([part1, part2], axis=1)
    pathogen_info = pathogen_info.loc[pathogen_info['synonyms'].str.contains('GTDB'),:].reset_index(drop=True)
    pathogen_info = pathogen_info.loc[pathogen_info['rank_list'].str.contains('species|strain'),:].reset_index(drop=True)
    pathogen_info.to_csv(os.path.join(args.output_dir, 'pathogen_info.tsv'), sep='\t', index=False)

    # find all potential pathogens (all strains under the known species-level pathogens that have at least one pathogen or themself are known as pathogens )
    species1 = set(pathogen_info.loc[pathogen_info['rank_list'].str.contains('strain'),'parent_list'].to_list())
    species2 = set(pathogen_info.loc[pathogen_info['rank_list'].str.contains('species'),'synonyms'].to_list())
    merged_species = species1 | species2
    all_potential_pathogenic_strain = set()
    for species in tqdm(merged_species):
        species_synonyms_list = eval(species)
        for species_synonyms in species_synonyms_list:
            all_potential_pathogenic_strain |= find_child[species_synonyms]
    microbe_nodes = node_info.loc[node_info['node_id'].str.contains('Microbe'), :].reset_index(drop=True)
    all_gtdb_id_to_microbe_id = dict()
    for row in tqdm(microbe_nodes.to_numpy()):
        node_id, node_type, all_names, synonyms, is_pathogen = row
        if 'GTDB:' not in synonyms:
            continue
        all_gtdb_id_to_microbe_id
        synonyms = eval(synonyms)
        for synonym in synonyms:
            if 'GTDB:' in synonym:
                all_gtdb_id_to_microbe_id[synonym] = (synonym, node_id, node_type, all_names, synonyms, is_pathogen)
    gtdbid_df = pd.DataFrame(all_gtdb_id_to_microbe_id.values(), columns=['gtdb_id', 'node_id', 'node_type', 'all_names', 'synonyms', 'is_pathogen'])
    # gtdbid_df.loc[gtdbid_df['gtdb_id'].isin(all_potential_pathogenic_strain),:].reset_index(drop=True).to_csv(os.path.join(args.output_dir, 'all_potential_pathogenic_strain.tsv'), sep='\t', index=False)
    all_potential_pathogenic_strain_df = gtdbid_df.loc[gtdbid_df['gtdb_id'].isin(all_potential_pathogenic_strain),:].reset_index(drop=True)
    gtdbid_to_nodeid = {row[0]:row[1] for row in gtdbid_df.to_numpy()}
    temp_df = pd.DataFrame(all_potential_pathogenic_strain_df['gtdb_id'].apply(lambda x: (find_parent[x], gtdbid_to_nodeid[x])).to_list())
    temp_df.columns = ['species_parent', 'node_id']
    all_potential_pathogenic_strain_df = pd.concat([all_potential_pathogenic_strain_df, temp_df], axis=1)
    all_potential_pathogenic_strain_df.to_csv(os.path.join(args.output_dir, 'all_potential_pathogenic_strain.tsv'), sep='\t', index=False)

    # randomly hold out 10% of species that have at least one strain-level pathogen
    species1_list = list(species1)
    random.shuffle(species1_list)
    hold_out_species = species1_list[:int(len(species1_list)*0.1)]
    pathogen_info.loc[pathogen_info['parent_list'].isin(hold_out_species),:].to_csv(os.path.join(args.output_dir, 'hold_out_pathogen_info.tsv'), sep='\t', index=False)

    # save edge info
    edge_info.to_csv(os.path.join(args.output_dir, 'edge_info.tsv'), sep='\t', index=False)    

    logger.info("Finished")

if __name__ == "__main__":
    main()
