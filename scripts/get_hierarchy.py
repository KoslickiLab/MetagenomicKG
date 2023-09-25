## Import standard libraries
import os
import sys
from tqdm import tqdm, trange
import pandas as pd
import argparse
import pytaxonkit
from typing import List, Dict, Tuple, Union, Any, Optional
import logging

## Import custom libraries
from utils import get_logger, check_files

def parse_taxonomy1(taxonomy, parent=None):
    """
    Parse the NCBI taxonomy file and return a list of rows.
    :param taxonomy: a dictionary of the taxonomy generated frm pytaxonkit.list()
    :param parent: a tuple of parent taxon id, rank, and name (default: None)
    :return: a list of rows
    """
    
    rows = []

    for key, children in taxonomy.items():
        taxon_id, rank_and_name = key.split(' [', 1)
        rank, taxon_name = rank_and_name.split('] ', 1)
        taxon_name = taxon_name.strip()

        if taxon_name == '':
            if parent:
                parent_id, parent_rank, parent_name = parent
                rows.extend(parse_taxonomy1(children, (parent_id, parent_rank, parent_name)))
        else:
            pass

        if parent:
            parent_id, parent_rank, parent_name = parent
            row = f"{parent_name}\t{parent_rank}\t{parent_id}\t{taxon_name}\t{rank}\t{taxon_id}"
            rows.append(row)

        rows.extend(parse_taxonomy1(children, (taxon_id, rank, taxon_name)))

    return rows

def parse_taxonomy2(taxonomy_file, taxonomy_metadata_file, logger, type='bacteria'):
    """
    Parse the GTDB taxonomy file and return a dataframe with the parent-child relationship.
    :param taxonomy_file: path of the taxonomy file
    :param taxonomy_metadata_file: path of the taxonomy metadata file
    :param logger: logger object
    :param type: type of the taxonomy, either 'bacteria' or 'archaea'
    :return: a dataframe with the parent-child relationship
    """
    logger.info(f"Parsing taxonomy for {type}...")
    
    table_rows = []
    mapping = {}
    abb_mapping = {'d': 'superkingdom', 'p': 'phylum', 'c': 'class', 'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species'}
    ## FIXME: ADD below line to fix a bug from GTDB
    def __fix_gtdb_bug1(x):
        return x.replace('Moorea','Moorena').replace('Nocardioides flavus','Nocardioides flavus Wang et al. 2016').replace('Sphingobacterium composti','Sphingobacterium composti Yoo et al. 2007 non Ten et al. 2007').replace('Fuerstia', 'Fuerstiella')
    def __fix_gtdb_bug2(x):
        return x.replace('Moorea','Moorena').replace('Fuerstia', 'Fuerstiella')

    ## read taxonomy file
    taxonomy = pd.read_csv(taxonomy_file, sep='\t', header=None)
    ## get all parent taxids
    all_parent_names = list(set([__fix_gtdb_bug1(y.split('__')[1]) for x in list(taxonomy[1]) for y in x.split(';')]))
    taxids = pytaxonkit.name2taxid(all_parent_names)
    ## distinguish which is is correct if there are multiple taxids
    duplicates = taxids.loc[taxids['Name'].isin(taxids['Name'][taxids['Name'].duplicated()]),:].reset_index(drop=True)
    if len(duplicates) > 0:
        taxid_fiter_out_duplicates = taxids[~taxids['Name'].isin(duplicates['Name'])].reset_index(drop=True)
        temp = pytaxonkit.lineage(list(duplicates['TaxID']))
        duplicates = duplicates.merge(temp, on=['Name','TaxID'], how='left').reset_index(drop=True)
        if type == 'bacteria':
            duplicates = duplicates.loc[(duplicates['Lineage'].str.contains('Bacteria;')) & (~duplicates['Lineage'].isna()),:].reset_index(drop=True)
        elif type == 'archaea':
            duplicates = duplicates.loc[(duplicates['Lineage'].str.contains('Archaea;')) & (~duplicates['Lineage'].isna()),:].reset_index(drop=True)
        else:
            logger.error(f"Type {type} is not supported.")
            sys.exit(1)
        taxids = pd.concat([taxid_fiter_out_duplicates[['Name','TaxID']], duplicates[['Name', 'TaxID']]]).reset_index(drop=True)
    else:
        taxids = taxids[['Name','TaxID']]

    ## convert taxids dataframe to dictionary
    temp_parent_mapping = {name:ncbi_taxid for name, ncbi_taxid in taxids.to_numpy()}
    mapping.update(temp_parent_mapping)
    
    ## read metadata file
    metadata = pd.read_csv(taxonomy_metadata_file, sep='\t', header=0)
    ## convert metadata dataframe to dictionary
    temp_mapping = {gtdb_id:ncbi_taxid for gtdb_id, ncbi_taxid in metadata[['accession','ncbi_taxid']].to_numpy()}    
    mapping.update(temp_mapping)
    
    for gtdb_strain, lineage in tqdm(taxonomy.to_numpy(), desc='Parsing taxonomy file'):
        lineage_list = [(x.split('__')[1], abb_mapping[x.split('__')[0]]) for x in lineage.split(';')] + [(gtdb_strain, 'strain')]
        # existing_taxids = {}
        remove_duplicates = {}
        for i in range(len(lineage_list)-1):
            parent = lineage_list[i]
            child = lineage_list[i+1]
            parent_taxid, parent_rank = mapping[__fix_gtdb_bug1(parent[0])], parent[1]
            # if isinstance(parent_taxid, int):
            #     if parent_taxid in existing_taxids and existing_taxids[parent_taxid] == parent[0]:
            #         pass
            #     elif parent_taxid in existing_taxids:
            #         parent_taxid = ''
            #     else:
            #         existing_taxids[parent_taxid] = parent[0]
            # else:
            #     parent_taxid = ''
            child_taxid, child_rank = mapping[__fix_gtdb_bug1(child[0])], child[1]
            # if isinstance(child_taxid, int):
            #     if child_taxid in existing_taxids and existing_taxids[child_taxid] == child[0]:
            #         pass
            #     elif child_taxid in existing_taxids:
            #         child_taxid = ''
            #     else:
            #         existing_taxids[child_taxid] = child[0]
            # else:
            #     child_taxid = ''
            if __fix_gtdb_bug2(parent[0]) == __fix_gtdb_bug2(child[0]):
                continue
            
            if __fix_gtdb_bug2(parent[0]) in remove_duplicates:
                parent_name, parent_rank, parent_taxid = remove_duplicates[__fix_gtdb_bug2(parent[0])]
            else:
                remove_duplicates[__fix_gtdb_bug2(parent[0])] = [__fix_gtdb_bug2(parent[0]), parent_rank, parent_taxid]
                parent_name, parent_rank, parent_taxid = __fix_gtdb_bug2(parent[0]), parent_rank, parent_taxid
            
            if __fix_gtdb_bug2(child[0]) in remove_duplicates:
                child_name, child_rank, child_taxid = remove_duplicates[__fix_gtdb_bug2(child[0])]
            else:
                remove_duplicates[__fix_gtdb_bug2(child[0])] = [__fix_gtdb_bug2(child[0]), child_rank, child_taxid]
                child_name, child_rank, child_taxid = __fix_gtdb_bug2(child[0]), child_rank, child_taxid
            
            table_rows.append((parent_name, parent_rank, parent_taxid, child_name, child_rank, child_taxid))
        
    
    out_df = pd.DataFrame(table_rows, columns=['Parent', 'Parent_rank', 'Parent_taxon_id', 'Child', 'Child_rank', 'Child_taxon_id']).drop_duplicates().reset_index(drop=True)

    return out_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get the taxonomy hierarchy of archaea, bacteria, fungi and viruses')
    parser.add_argument('--bacteria_taxonomy', type=str, help='path of the GTDB taxonomy file for bacteria (e.g., bac120_taxonomy_r207.tsv)', required=True)
    parser.add_argument('--bacteria_metadata', type=str, help='path of the GTDB metadata file for archaea (e.g., bac120_metadata_r207.tsv)', required=True)
    parser.add_argument('--archaea_taxonomy', type=str, help='path of the GTDB taxonomy file for archaea (e.g., ar53_taxonomy_r207.tsv)', required=True)
    parser.add_argument('--archaea_metadata', type=str, help='path of the GTDB metadata file for archaea (e.g., ar53_metadata_r207.tsv)', required=True)
    parser.add_argument('--output_dir', type=str, help='path of the output directory')
    args = parser.parse_args()

    # Create a logger object
    logger = get_logger()
    logger.setLevel(logging.INFO)

    # Generate a edge list for bacteria hierarchy
    logger.info('Parsing the taxonomy file for bacteria')
    if check_files(args.bacteria_taxonomy, logger) and check_files(args.bacteria_metadata, logger):
        bacteria_df = parse_taxonomy2(args.bacteria_taxonomy, args.bacteria_metadata, logger, 'bacteria')
    else:
        sys.exit(1)
    # Save the edge list to a file
    logger.info('Save the edge list to a file')
    bacteria_df.to_csv(os.path.join(args.output_dir, 'bacteria_hierarchy.tsv'), sep='\t', index=False)
    
    # Generate a edge list for archaea hierarchy
    logger.info('Generate a edge list for archaea hierarchy')
    if check_files(args.archaea_taxonomy, logger) and check_files(args.archaea_metadata, logger):
        archaea_df = parse_taxonomy2(args.archaea_taxonomy, args.archaea_metadata, logger, 'archaea')
    else:
        sys.exit(1)
    # Save the edge list to a file
    logger.info('Save the edge list to a file')
    archaea_df.to_csv(os.path.join(args.output_dir, 'archaea_hierarchy.tsv'), sep='\t', index=False)

    # Generate a edge list for fungi hierarchy
    logger.info('Generate a edge list for fungi hierarchy')
    result = pytaxonkit.list([4751], raw=True)
    table_rows = table_data = [row.split('\t') for row in parse_taxonomy1(result)]
    # Create a DataFrame from the table_data list
    logger.info('Create a DataFrame from the table_data list')
    fungi_df = pd.DataFrame(table_rows, columns=['Parent', 'Parent_rank', 'Parent_taxon_id', 'Child', 'Child_rank', 'Child_taxon_id']).drop_duplicates().reset_index(drop=True)
    # Save the edge list to a file
    logger.info('Save the edge list to a file')
    fungi_df.to_csv(os.path.join(args.output_dir, 'fungi_hierarchy.tsv'), sep='\t', index=False)

    # Generate a edge list for viruses hierarchy
    logger.info('Generate a edge list for viruses hierarchy')
    result = pytaxonkit.list([10239], raw=True)
    table_rows = table_data = [row.split('\t') for row in parse_taxonomy1(result)]
    # Create a DataFrame from the table_data list
    logger.info('Create a DataFrame from the table_data list')
    viruses_df = pd.DataFrame(table_rows, columns=['Parent', 'Parent_rank', 'Parent_taxon_id', 'Child', 'Child_rank', 'Child_taxon_id']).drop_duplicates().reset_index(drop=True)
    # Save the edge list to a file
    logger.info('Save the edge list to a file')
    viruses_df.to_csv(os.path.join(args.output_dir, 'viruses_hierarchy.tsv'), sep='\t', index=False)

