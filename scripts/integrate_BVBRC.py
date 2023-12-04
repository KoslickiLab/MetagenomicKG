"""
This script is used to integrate all BV-BRC data into a knolwedge graph.
"""

## Import standard libraries
import os
import sys
from tqdm import tqdm, trange
from glob import glob
import pandas as pd
import argparse
import pytaxonkit
import requests
import logging

## Import custom libraries
from utils import get_logger, read_tsv_file, Node, Edge, KnowledgeGraph, extract_disease_synonyms, UMLSMapping


def get_new_name(disease_name):
    
    if disease_name == 'Seticemic plague':
        disease_name = 'Septicemic plague'
    elif disease_name == 'empiem':
        disease_name = 'empyema'
    elif disease_name == "Legionaires' disease":
        disease_name = "Legionnaires' Disease"
    elif disease_name == 'Bone and joint infection':
        disease_name = 'Bone and joint infections'
    elif disease_name == 'meniigitis':
        disease_name = 'meningitis'
    elif disease_name == 'Peridontitis':
        disease_name = 'Periodontitis'
    elif disease_name == 'Gastrointeritis':
        disease_name = 'Gastroenteritis'
    elif disease_name == 'Wide range of infections':
        disease_name = 'Infection'
    elif disease_name == 'severe pneumonia':    
        disease_name = 'Pneumonia'
    elif disease_name == 'bacteremic pneumonia':    
        disease_name = 'Pneumonia'
    elif disease_name == 'streptococcal toxic shock syndrome (STSS)':
        disease_name = 'streptococcal toxic shock syndrome'
    elif disease_name == 'Various infections':
        disease_name = 'Infection'
    elif disease_name == 'bacteriemia':
        disease_name = 'Bacteremia'
    elif disease_name == 'Chorioamnioitis':
        disease_name = 'Chorioamnionitis'
    elif disease_name == 'Gaslrointestinal perforation':
        disease_name = 'Gastrointestinal perforation'
    elif disease_name == 'respritory tract infection':
        disease_name = 'Respiratory Tract Infection'
    elif disease_name == 'inflammatory Diarrheal disease':
        disease_name = 'Diarrheal disorder'
    elif disease_name == 'Sepsis of The Newbornl':
        disease_name = 'Sepsis of The Newborn'
    elif disease_name == 'haematological malignancies':
        disease_name = 'hematological malignancies'
    elif disease_name == 'Hemoptoic pneumonia':
        disease_name = 'Pneumonia'
    elif disease_name == 'Necrotizing faciitis':
        disease_name = 'Necrotizing fasciitis'
    elif disease_name == 'bacterimia':
        disease_name = 'Bacteremia'
    else:
        disease_name = disease_name
                    
    return disease_name

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Integrate all BV-BRC data into a knolwedge graph')
    parser.add_argument('--existing_KG_nodes', type=str, help='path of the existing knowledge graph nodes')
    parser.add_argument('--existing_KG_edges', type=str, help='path of the existing knowledge graph edges')
    parser.add_argument('--data_dir', type=str, help='path of the BV-BRC data directory')
    parser.add_argument('--umls_api_key', type=str, help='UMLS API key')
    parser.add_argument('--gtdb_assignment', type=str, help='path of the gtdb assignment using GTDB-Tk')
    parser.add_argument('--synonymizer_dir', type=str, help='path of the synonymizer directory')
    parser.add_argument('--synonymizer_dbname', type=str, help='name of the synonymizer database')
    parser.add_argument('--ANI_threshold', type=float, help='ANI threshold to identify the same strain (default 0 for no filtering)', default=0.0)
    parser.add_argument('--AF_threshold', type=float, help='AF threshold to dentify the same strain (default 0 for no filtering)', default=0.0)
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
    kg.load_graph(load_dir=args.output_dir, node_filename=node_filename, edge_filename=edge_filename)

    # Import Node Synonymizer
    logger.info("Importing Node Synonymizer...")
    sys.path.append(args.synonymizer_dir)
    from node_synonymizer import NodeSynonymizer
    nodesynonymizer = NodeSynonymizer(args.synonymizer_dir, args.synonymizer_dbname)
    
    # set up umls mapping
    logger.info('Setting up UMLS mapping')
    umls_class = UMLSMapping(args.umls_api_key)
    
    # Import GTDB assignment data
    logger.info("Importing GTDB assignment data...")
    gtdb_assignment = read_tsv_file(args.gtdb_assignment)
    temp_header = gtdb_assignment[0]
    mapping_gn_to_microbe_id = {}
    for row in tqdm(gtdb_assignment[1:]):
        if (row[temp_header.index('fastani_ani')] != 'N/A' and float(row[temp_header.index('fastani_ani')]) >= args.ANI_threshold) and (row[temp_header.index('fastani_af')] != 'N/A' and float(row[temp_header.index('fastani_af')]) >= args.AF_threshold):
            mapping_gn_to_microbe_id[f"BVBRC:gn_{row[temp_header.index('user_genome')]}"] = [f"GTDB:{row[temp_header.index('fastani_reference')]}"]
        else:
            mapping_gn_to_microbe_id[f"BVBRC:gn_{row[temp_header.index('user_genome')]}"] = [f"BVBRC:gn_{row[temp_header.index('user_genome')]}"]
        mapping_gn_to_microbe_id[f"BVBRC:gn_{row[temp_header.index('user_genome')]}"] += [row[temp_header.index('classification')], f"ANI_reference_radius:{row[temp_header.index('fastani_reference_radius')]}; ANI:{row[temp_header.index('fastani_ani')]}; AF:{row[temp_header.index('fastani_af')]}; classification_method:{row[temp_header.index('classification_method')]}; MSA_percent:{row[temp_header.index('msa_percent')]}"]

    # Import BV-BRC data
    logger.info("Importing BV-BRC data...")
    bvbrc_relation_data = read_tsv_file(os.path.join(args.data_dir, 'output_subset_human_related_cleaned_metadata.tsv'))
    temp_header = bvbrc_relation_data[0]
    bvbrc_relation_data = pd.DataFrame(bvbrc_relation_data[1:], columns=temp_header)
    bvbrc_relation_data = bvbrc_relation_data.loc[bvbrc_relation_data['disease'] != '',:].reset_index(drop=True)
    bvbrc_gn_meta_data = read_tsv_file(os.path.join(args.data_dir, 'genome_metadata.txt'))
    temp_header = bvbrc_gn_meta_data[0]
    bvbrc_gn_meta_data = pd.DataFrame(bvbrc_gn_meta_data[1:], columns=temp_header)
    bvbrc_gn_meta_data = bvbrc_gn_meta_data.loc[bvbrc_gn_meta_data['genome_id'].isin(list(bvbrc_relation_data['genome_id'])),:].reset_index(drop=True)
    ## FIXME: add more metainfo
    bvbrc_relation_data = bvbrc_relation_data.merge(bvbrc_gn_meta_data[['genome_id', 'genome_name', 'organism_name', 'taxon_id']], on='genome_id', how='left').reset_index(drop=True)
    ## Add NCBI lineage to distinguish between vrus/fungi and bacteria/archaea
    result = pytaxonkit.lineage(list(set(bvbrc_relation_data['taxon_id'])))
    result['TaxID'] = result['TaxID'].astype(str)
    bvbrc_relation_data = bvbrc_relation_data.merge(result[['TaxID','Rank','Lineage']], left_on='taxon_id', right_on='TaxID', how='left').reset_index(drop=True)
    bvbrc_relation_data.drop(columns=['TaxID'], inplace=True)
    bvbrc_relation_data = bvbrc_relation_data.loc[bvbrc_relation_data['Lineage'].str.contains('Archaea|Bacteria|Viruses'),:].reset_index(drop=True)
    # if len(bvbrc_relation_data.loc[~bvbrc_relation_data['Lineage'].str.contains('Archaea|Bacteria|Viruses'),:]) > 0:
    #     raise ValueError(f"Found non-bacterial/archaeal genomes")

    ## Extract Disease Synonyms
    mapping_disease_to_synonyms = {}
    logger.info("Extracting Disease Synonyms...")
    unknown = []
    for disease_name_list1 in tqdm(set(bvbrc_relation_data['disease'].tolist())):
        for disease_name_list2 in disease_name_list1.split(';'):
            for disease_name_list3 in disease_name_list2.split(','):
                for disease_name in disease_name_list3.split('/'):
                    if disease_name == '':
                        continue
                    disease_name = disease_name.strip()
                    if disease_name in ['other', 'Not connected', 'agn', 'and sepsis', 'and soft tissue infections.', 'healthy', 'healthy controls', 'healty', 'lrti', 'pcs']:
                        continue
                    # Fix dirty data:
                    disease_name = get_new_name(disease_name)
                    if disease_name == 'Urinary tract and respiratory infections':
                        for x in ['Respiratory Tract Infections', 'Urinary tract infection']:
                            synonym_res = extract_disease_synonyms(x, nodesynonymizer, umls_class, args.umls_api_key)
                            if len(synonym_res) != 0:
                                mapping_disease_to_synonyms[disease_name] = synonym_res
                            else:
                                logger.warning(f"Cannot find synonyms for {disease_name}")
                                unknown += [disease_name]
                        continue
                    synonym_res = extract_disease_synonyms(disease_name, nodesynonymizer, umls_class, args.umls_api_key)
                    if len(synonym_res) != 0:
                        mapping_disease_to_synonyms[disease_name] = synonym_res
                    else:
                        logger.warning(f"Cannot find synonyms for {disease_name}")
                        unknown += [disease_name]
                
    ## Map the BVBRC genomes to KEGG nodes
    parent_dict = {}
    for row in tqdm(bvbrc_relation_data.to_numpy(), desc="Map the BVBRC genomes to KEGG nodes"):
        genome_id, _, _, _, _, _, _, assembly_accession, _, disease, genome_name, _, taxon_id, ncbi_rank, lineage = row
        if assembly_accession != '' and kg.find_node_by_synonym(f"GTDB:{assembly_accession}"):
            node_id = kg.find_node_by_synonym(f"GTDB:{assembly_accession}")
            existing_node = kg.get_node_by_id(node_id)
            kg.map_synonym_to_node_id[f"BVBRC:gn_{genome_id}"] = node_id
            existing_node.synonyms = list(set(existing_node.synonyms + [f"BVBRC:gn_{genome_id}"]))
            description_dict = dict(existing_node.description)
            if 'taxid' in description_dict and description_dict['taxid'] == '':
                description_dict['taxid'] = taxon_id
                description_dict['rank'] = ncbi_rank
            elif 'taxid' not in description_dict:
                description_dict['taxid'] = taxon_id
                description_dict['rank'] = ncbi_rank
            if assignment_info[2].split(';')[1].split(':')[1] != 'N/A':
                description_dict.update(dict([tuple(item.strip().split(':')) for item in assignment_info[2].split(';')]))
            existing_node.description = list(description_dict.items())
            existing_node.link = list(set(existing_node.link + [f"https://gtdb.ecogenomic.org/genome?gid={assembly_accession}", f"https://www.ncbi.nlm.nih.gov/assembly/{assembly_accession}", f"https://www.bv-brc.org/view/Genome/{genome_id}"]))
            existing_node.all_names = list(set(existing_node.all_names + [genome_name]))
            existing_node.knowledge_source = list(set(existing_node.knowledge_source + ['BVBRC']))
            existing_node.is_pathogen = True
        elif f"BVBRC:gn_{genome_id}" in mapping_gn_to_microbe_id:
            assignment_info = mapping_gn_to_microbe_id[f"BVBRC:gn_{genome_id}"]
            node_id = kg.find_node_by_synonym(assignment_info[0])
            if node_id:
                existing_node = kg.get_node_by_id(node_id)
                kg.map_synonym_to_node_id[f"BVBRC:gn_{genome_id}"] = node_id
                existing_node.synonyms = list(set(existing_node.synonyms + [f"BVBRC:gn_{genome_id}"]))
                description_dict = dict(existing_node.description)
                if 'taxid' in description_dict and description_dict['taxid'] == '':
                    description_dict['taxid'] = taxon_id
                    description_dict['rank'] = ncbi_rank
                elif 'taxid' not in description_dict:
                    description_dict['taxid'] = taxon_id
                    description_dict['rank'] = ncbi_rank
                if assignment_info[2].split(';')[1].split(':')[1] != 'N/A':
                    description_dict.update(dict([tuple(item.strip().split(':')) for item in assignment_info[2].split(';')]))
                existing_node.description = list(description_dict.items())
                existing_node.link = list(set(existing_node.link + [f"https://www.bv-brc.org/view/Genome/{genome_id}"]))
                if assembly_accession != '':
                    existing_node.link += [f"https://gtdb.ecogenomic.org/genome?gid={assembly_accession}", f"https://www.ncbi.nlm.nih.gov/assembly/{assembly_accession}"]
                existing_node.all_names = list(set(existing_node.all_names + [genome_name]))
                existing_node.knowledge_source = list(set(existing_node.knowledge_source + ['BVBRC']))
                existing_node.is_pathogen = True
            elif 'Unclassified' in assignment_info[1]:
                # This genome can't be unclassified based on the GTDB hierarchy
                continue
            else:
                temp_all_names = [genome_name]
                temp_description_dict = {'taxid': taxon_id, 'rank': ncbi_rank}
                if assignment_info[2].split(';')[1].split(':')[1] != 'N/A':
                    temp_description_dict.update(dict([tuple(item.strip().split(':')) for item in assignment_info[2].split(';')]))
                temp_knowledge_source = ['BVBRC']
                temp_synonyms = [f"BVBRC:gn_{genome_id}"]
                temp_link = [f"https://www.bv-brc.org/view/Genome/{genome_id}"]
                if assembly_accession != '':
                    temp_synonyms += [f"GTDB:{assembly_accession}"]
                    temp_link += [f"https://gtdb.ecogenomic.org/genome?gid={assembly_accession}", f"https://www.ncbi.nlm.nih.gov/assembly/{assembly_accession}"]
                temp_node = Node(node_type="Microbe", all_names=temp_all_names, description=list(temp_description_dict.items()), knowledge_source=temp_knowledge_source, synonyms=temp_synonyms, link=temp_link, is_pathogen=True)
                kg.add_node(temp_node)
                gtdb_classification = assignment_info[1]              
                parent_dict[f"BVBRC:gn_{genome_id}"] = f"GTDB:{[x for x in gtdb_classification.split(';')[::-1] if x.split('__')[1] != ''][0].split('__')[1]}"
        elif 'Viruses;' in lineage:            
            node_id = kg.find_node_by_synonym(f"NCBI:{genome_name}")
            existing_node = kg.get_node_by_id(node_id)
            kg.map_synonym_to_node_id[f"BVBRC:gn_{genome_id}"] = node_id
            existing_node.synonyms = list(set(existing_node.synonyms + [f"BVBRC:gn_{genome_id}"]))
            description_dict = dict(existing_node.description)
            if 'taxid' in description_dict and description_dict['taxid'] == '':
                description_dict['taxid'] = taxon_id
                description_dict['rank'] = ncbi_rank
            elif 'taxid' not in description_dict:
                description_dict['taxid'] = taxon_id
                description_dict['rank'] = ncbi_rank
            existing_node.description = list(description_dict.items())
            existing_node.link = list(set(existing_node.link + [f"https://www.ncbi.nlm.nih.gov/assembly/{assembly_accession}", f"https://www.bv-brc.org/view/Genome/{genome_id}"]))
            existing_node.all_names = list(set(existing_node.all_names + [genome_name]))
            existing_node.knowledge_source = list(set(existing_node.knowledge_source + ['BVBRC']))
            existing_node.is_pathogen = True
        else:
            logger.warning(f"Cannot find {genome_id} in the GTDB assignment data")
            
    ## Connect genome hierarchy
    logger.info('Connecting genome hierarchy')
    for genome_child in parent_dict:
        genome_parent = parent_dict[genome_child]
        kg_source = 'BVBRC'
        
        # Add edge to the knowledge graph
        kg.add_edge(Edge(source_node=genome_parent, target_node=genome_child, predicate='biolink:has_part', knowledge_source=[kg_source]))
        kg.add_edge(Edge(source_node=genome_child, target_node=genome_parent, predicate='biolink:part_of', knowledge_source=[kg_source]))
    
    
    ## Connect the BVBRC genomes to BVBRC disease
    for row in tqdm(bvbrc_relation_data.to_numpy(), desc="Connecting BVBRC genomes to BVBRC disease"):
        genome_id, _, _, _, _, _, _, _, _, disease_list, _, _, _, _, _ = row
        kg_source = 'BVBRC'
        for disease_name_list1 in disease_list.split(';'):
            for disease_name_list2 in disease_name_list1.split(','):
                for disease_name in disease_name_list2.split('/'):
                    disease_name = disease_name.strip()
                    genome_node_id = kg.find_node_by_synonym(f"BVBRC:gn_{genome_id}")
                    if not genome_node_id:
                        continue
                    if disease_name not in mapping_disease_to_synonyms:
                        continue
                    disease_ids = list(set([kg.find_node_by_synonym(x) for x in mapping_disease_to_synonyms[disease_name] if kg.find_node_by_synonym(x)]))
                    if len(disease_ids) == 0:
                        continue
                    if len(disease_ids) == 1:
                        disease_node_id = disease_ids[0]
                        # Add edge to the knowledge graph
                        kg.add_edge(Edge(source_node=genome_node_id, target_node=disease_node_id, predicate='biolink:associated_with', knowledge_source=[kg_source]))
                        kg.add_edge(Edge(source_node=disease_node_id, target_node=genome_node_id, predicate='biolink:associated_with', knowledge_source=[kg_source]))
                    else:
                        logger.warning(f"Multiple disease nodes found for {disease_name}")
                        for disease_node_id in disease_ids:
                            # Add edge to the knowledge graph
                            kg.add_edge(Edge(source_node=genome_node_id, target_node=disease_node_id, predicate='biolink:associated_with', knowledge_source=[kg_source]))
                            kg.add_edge(Edge(source_node=disease_node_id, target_node=genome_node_id, predicate='biolink:associated_with', knowledge_source=[kg_source]))

    # Save the knowledge graph
    logger.info("Saving the knowledge graph...")
    kg.save_graph(save_dir = args.output_dir, node_filename = 'KG_nodes_v4.tsv', edge_filename = 'KG_edges_v4.tsv')
    logger.info("KG node is saved to {}".format(os.path.join(args.output_dir, 'KG_nodes_v4.tsv')))
    logger.info("KG edge is saved to {}".format(os.path.join(args.output_dir, 'KG_edges_v4.tsv')))
    
    logger.info(f'Done!')