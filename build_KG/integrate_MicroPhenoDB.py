"""
This script is used to integrate all MicroPhenoDB data into a knolwedge graph.
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
import json
import itertools

## Import custom libraries
from utils import get_logger, read_tsv_file, Node, Edge, KnowledgeGraph, change_prefix, extract_disease_synonyms, UMLSMapping, OxOMapping, MappingLink
from kg2_utils.node_synonymizer import NodeSynonymizer

def map_umls_to_taxon_id(ncit, umls_id, apikey):
    
    def _parse_response(response_json, preferred_name):
        temp = list(set([int(res['code'].split('/')[-1]) for res in response_json['result'] if res['rootSource'] == 'NCBI' and res['name'].lower() == preferred_name.lower()]))
        if len(temp) == 0:
            temp = list(set([int(res['code'].split('/')[-1]) for res in response_json['result'] if res['rootSource'] == 'NCBI']))
            return temp
        else:
            return temp
    
    if isinstance(ncit, str):        
        umls_api = f"https://uts-ws.nlm.nih.gov/rest/content/current/source/NCI/{ncit.split('_')[1]}/atoms"
        response = requests.get(umls_api, params={'apiKey':apikey})
        if response.status_code == 200:
            umls_id = response.json()['result'][0]['concept'].split('/')[-1]
            umls_api = f"https://uts-ws.nlm.nih.gov/rest/content/current/CUI/{umls_id}"
            response = requests.get(umls_api, params={'apiKey':apikey})
            if response.status_code == 200:
                preferred_name = response.json()['result']['name']
            else:
                return None

            umls_api = f"https://uts-ws.nlm.nih.gov/rest/content/current/CUI/{umls_id}/atoms"
            response = requests.get(umls_api, params={'apiKey':apikey})
            if response.status_code == 200:
                res = _parse_response(response.json(), preferred_name)
                if len(res) > 0:
                    return res
                else:
                    return None
            else:
                return None
        else:
            return None

    elif umls_id and umls_id.split(':')[0] == 'UMLS':
        umls_api = f"https://uts-ws.nlm.nih.gov/rest/content/current/CUI/{umls_id.split(':')[1]}"
        response = requests.get(umls_api, params={'apiKey':apikey})
        if response.status_code == 200:
            preferred_name = response.json()['result']['name']
        else:
            return None

        umls_api = f"https://uts-ws.nlm.nih.gov/rest/content/current/CUI/{umls_id.split(':')[1]}/atoms"
        response = requests.get(umls_api, params={'apiKey':apikey})
        if response.status_code == 200:
            res = _parse_response(response.json(), preferred_name)
            if len(res) > 0:
                return res
            else:
                return None
        else:
            return None
    else:
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Integrate all MicroPhenoDB into a knolwedge graph')
    parser.add_argument('--existing_KG_nodes', type=str, help='path of the existing knowledge graph nodes')
    parser.add_argument('--existing_KG_edges', type=str, help='path of the existing knowledge graph edges')
    # parser.add_argument('--bacteria_metadata', type=str, help='path of the GTDB metadata file for archaea (e.g., bac120_metadata_r207.tsv)')
    # parser.add_argument('--archaea_metadata', type=str, help='path of the GTDB metadata file for archaea (e.g., ar53_metadata_r207.tsv)')
    parser.add_argument('--data_dir', type=str, help='path of the MicroPhenoDB data directory')
    parser.add_argument('--umls_api_key', type=str, help='UMLS API key')
    parser.add_argument('--synonymizer_dir', type=str, help='path of the synonymizer directory')
    parser.add_argument('--synonymizer_dbname', type=str, help='name of the synonymizer database')
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
    nodesynonymizer = NodeSynonymizer(args.synonymizer_dir, args.synonymizer_dbname)

    # set up oxo mappingf
    logger.info('Setting up OxO mapping')
    oxo_class = OxOMapping()
    
    # set up umls mapping
    logger.info('Setting up UMLS mapping')
    umls_class = UMLSMapping(args.umls_api_key)

    # # Read bacteria metadata
    # logger.info("Reading bacteria metadata...")
    # temp_data = read_tsv_file(args.bacteria_metadata)
    # bacteria_meta_df = pd.DataFrame(temp_data[1:], columns=temp_data[0])
        
    # # Read archaea metadata
    # logger.info("Reading archaea metadata...")
    # temp_data = read_tsv_file(args.archaea_metadata)
    # archaea_meta_df = pd.DataFrame(temp_data[1:], columns=temp_data[0])

    # gtdb_meta_df = pd.concat([bacteria_meta_df, archaea_meta_df])[['accession','ncbi_taxid']].reset_index(drop=True)

    # Load MicroPhenoDB data
    logger.info('Loading MicroPhenoDB data')
    if not os.path.exists(os.path.join(args.data_dir,'core_table.txt')):
        logger.error(f'Could not find MicroPhenoDB core table at {os.path.join(args.data_dir,"core_table.txt")}')
        sys.exit(1)
    else:
        MDBcorrelation = read_tsv_file(os.path.join(args.data_dir,'core_table.txt'))
        temp_header = MDBcorrelation[0]
        MDBcorrelation = pd.DataFrame(MDBcorrelation[1:], columns=temp_header)

    if not os.path.exists(os.path.join(args.data_dir,'EFO.txt')):
        logger.error(f'Could not find MicroPhenoDB EFO table at {os.path.join(args.data_dir,"EFO.txt")}')
        sys.exit(1)
    else:
        # disease_table = pd.read_csv(os.path.join(args.data_dir,'EFO.txt'), sep='\t', header=0, encoding = "ISO-8859-1")
        disease_table = read_tsv_file(os.path.join(args.data_dir,'EFO.txt'), encoding = "ISO-8859-1")
        temp_header = disease_table[0]
        disease_table = pd.DataFrame(disease_table[1:], columns=temp_header)
        
    if not os.path.exists(os.path.join(args.data_dir,'NCIT.txt')):
        logger.error(f'Could not find MicroPhenoDB NCIT table at {os.path.join(args.data_dir,"NCIT.txt")}')
        sys.exit(1)
    else:
        # species_table = pd.read_csv(os.path.join(args.data_dir,'NCIT.txt'), sep='\t', header=0)
        species_table = read_tsv_file(os.path.join(args.data_dir,'NCIT.txt'))
        temp_header = species_table[0]
        species_table = pd.DataFrame(species_table[1:], columns=temp_header)

    ## Generate disease and species mapping dict
    logger.info('Generating disease and species mapping dict')
    disease_dict = {}
    species_dict = {}
    MDBcorrelation = MDBcorrelation.merge(disease_table[['Scientific_disease_name', 'EFO_id', 'Disease_annotation']], left_on='Disease', right_on='Scientific_disease_name', how='left')
    MDBcorrelation = MDBcorrelation.merge(species_table[['Scientific_species_name', 'NCIT_id', 'Species_annotation']], left_on='Microbe', right_on='Scientific_species_name', how='left')
    
    
    # find taxon id for species
    taxon_names = list(set(MDBcorrelation[['Microbe']]['Microbe']))
    name_to_ncit = {name:ncit_id for name, ncit_id in MDBcorrelation[['Microbe', 'NCIT_id']].to_numpy()}
    umls_ids = {name:umls_class.get_mapping(name) for name in taxon_names}
    name_umls = pd.DataFrame([(key, value) for key, value in umls_ids.items() if value], columns=['name', 'umls_id'])
    name_umls_taxonid = pd.DataFrame([(name, umls, map_umls_to_taxon_id(name_to_ncit[name], umls, args.umls_api_key)) for name, umls in name_umls.to_numpy()])
    # use node synonymizer to find the taxon ids
    temp = []
    for name, umls, taxon_id in name_umls_taxonid.to_numpy():
        if taxon_id:
            if len(taxon_id) == 1:
                temp.append((name, umls, str(taxon_id[0])))
            else:
                logger.warning(f'Found multiple taxon ids for {name} with UMLS id {umls}: {taxon_id}')
                temp.append((name, umls, str(taxon_id[0])))
        else:
            normalizer = nodesynonymizer.get_equivalent_nodes(names=name)[name]
            if normalizer:
                taxon_id = [int(x.split(':')[1]) for x in normalizer if x.split(':')[0] == 'NCBITaxon']
                if len(taxon_id) == 0:
                    temp.append((name, umls, None))
                elif len(taxon_id) == 1:
                    temp.append((name, umls, str(taxon_id[0])))
                else:
                    logger.warning(f'Found multiple taxon ids for {name} with UMLS id {umls}: {taxon_id}')
                    temp.append((name, umls, str(taxon_id[0])))
    name_umls_taxonid = pd.DataFrame(temp, columns=['name', 'umls_id', 'taxon_id'])
    name_umls_taxonid = name_umls_taxonid.loc[~name_umls_taxonid['taxon_id'].isna(),:].reset_index(drop=True)
    temp_result = pytaxonkit.lineage([int(x) for x in name_umls_taxonid['taxon_id']])
    temp_result = temp_result[['TaxID','Name','Lineage','Rank','FullLineage']]
    temp_result['TaxID'] = temp_result['TaxID'].astype(str)
    name_umls_taxonid = name_umls_taxonid.merge(temp_result, left_on='taxon_id', right_on='TaxID', how='left')
    name_umls_taxonid_dict = {row[0]: row[1:] for row in name_umls_taxonid.to_numpy()}
    
    
    for row in tqdm(MDBcorrelation.to_numpy(), desc='Generating disease and species mapping dict'):
        _, microbe, disease, _, pmid, _, _, _, disease_name, efo_id, disease_annotation, species_name, ncit_id, species_annotation = row
        
        
        ## map disease to the node in kg
        if disease in ['Disease', 'Disease-free', 'Null', 'Not foundthogenic']:
            continue
        elif disease in ['Pelvic inflamm atory disease', 'Pelvic in铿俛mmatory disease']:
            disease = 'Pelvic inflammatory disease'

        if not isinstance(efo_id, str):
            continue

        efo_id = str(efo_id)
        efo_id = efo_id.replace('_',':')
        synonyms = oxo_class.get_mappings(efo_id, distance=1)    
        if synonyms and len(synonyms[efo_id]) > 0:
            synonyms = synonyms[efo_id]
        else:
            continue
        
        synonyms = [efo_id] + [change_prefix(synonym) for synonym in synonyms if change_prefix(synonym)]
        disease_node_ids = list(set([kg.find_node_by_synonym(x) for x in synonyms if kg.find_node_by_synonym(x)]))
        if len(disease_node_ids) == 0:
            synonyms = extract_disease_synonyms(disease, nodesynonymizer, umls_class, args.umls_api_key)
            if len(synonyms) != 0:
                disease_node_ids = list(set([kg.find_node_by_synonym(x) for x in synonyms if kg.find_node_by_synonym(x)]))
                if len(disease_node_ids) == 0:
                    logger.warning(f"Cannot find synonyms for {disease}")
                    continue
                elif len(disease_node_ids) == 1:
                    node_id = disease_node_ids[0]
                    existing_node = kg.get_node_by_id(node_id)
                    for synonym in synonyms:
                        kg.map_synonym_to_node_id[synonym] = node_id
                    existing_node.synonyms = list(set(existing_node.synonyms + synonyms))
                    description_dict = dict(existing_node.description)
                    if isinstance(disease_annotation, str) and disease_annotation !='':
                        description_dict['MicroPhenoDB Description'] = disease_annotation
                    existing_node.description = list(description_dict.items())
                    temp_link = [MappingLink(x) for x in synonyms if MappingLink(x)]
                    existing_node.link = list(set(existing_node.link + temp_link))
                    existing_node.all_names = list(set(existing_node.all_names + [disease_name]))
                    existing_node.knowledge_source = list(set(existing_node.knowledge_source + ['MicroPhenoDB']))
                else:
                    logger.warning(f"Multiple nodes found for {disease}")
            else:
                logger.warning(f"Cannot find synonyms for {disease}")
                continue
        elif len(disease_node_ids) == 1:
            node_id = disease_node_ids[0]
            existing_node = kg.get_node_by_id(node_id)
            for synonym in synonyms:
                kg.map_synonym_to_node_id[synonym] = node_id
            existing_node.synonyms = list(set(existing_node.synonyms + synonyms))
            description_dict = dict(existing_node.description)
            if isinstance(disease_annotation, str) and disease_annotation !='':
                description_dict['MicroPhenoDB Description'] = disease_annotation
            existing_node.description = list(description_dict.items())
            temp_link = [MappingLink(x) for x in synonyms if MappingLink(x)]
            existing_node.link = list(set(existing_node.link + temp_link))
            existing_node.all_names = list(set(existing_node.all_names + [disease_name]))
            existing_node.knowledge_source = list(set(existing_node.knowledge_source + ['MicroPhenoDB']))
        else:
            logger.warning(f"Multiple disease nodes found for {disease}")
    

        ## map species to the node in kg
        if microbe not in name_umls_taxonid_dict:
            continue
        else:
            umls_id, taxon_id, _, ncbi_name, _, rank, fullLineage = name_umls_taxonid_dict[microbe]
            if 'Bacteria' in fullLineage or 'Archaea' in fullLineage:
                species_node_ids = kg.find_node_by_synonym(f"GTDB:{ncbi_name}")
                if not species_node_ids:
                    logger.warning(f"Cannot find node for {microbe}")
                    continue
                existing_node = kg.get_node_by_id(species_node_ids)
                description_dict = dict(existing_node.description)
                if isinstance(species_annotation, str) and species_annotation !='':
                    description_dict['MicroPhenoDB Description'] = species_annotation
                existing_node.description = list(description_dict.items())
                existing_node.knowledge_source = list(set(existing_node.knowledge_source + ['MicroPhenoDB']))
                existing_node.is_pathogen = True
            elif 'Viruses' in fullLineage or 'Fungi' in fullLineage:
                species_node_ids = kg.find_node_by_synonym(f"NCBI:{ncbi_name}")
                if not species_node_ids:
                    logger.warning(f"Cannot find node for {microbe}")
                    continue
                existing_node = kg.get_node_by_id(species_node_ids)
                description_dict = dict(existing_node.description)
                if isinstance(species_annotation, str) and species_annotation !='':
                    description_dict['MicroPhenoDB Description'] = species_annotation
                existing_node.description = list(description_dict.items())
                existing_node.knowledge_source = list(set(existing_node.knowledge_source + ['MicroPhenoDB']))
                existing_node.is_pathogen = True
            else:
                continue

        ## add edge to KG
        kg_source = 'MicroPhenoDB'
        for pair in itertools.product([species_node_ids], disease_node_ids):
            kg.add_edge(Edge(source_node=pair[0], target_node=pair[1], predicate='biolink:associated_with', knowledge_source=[kg_source]))
            kg.add_edge(Edge(source_node=pair[1], target_node=pair[0], predicate='biolink:associated_with', knowledge_source=[kg_source]))


    # Save the knowledge graph
    logger.info("Saving the knowledge graph...")
    kg.save_graph(save_dir = args.output_dir, node_filename = 'KG_nodes_v5.tsv', edge_filename = 'KG_edges_v5.tsv')
    logger.info("KG node is saved to {}".format(os.path.join(args.output_dir, 'KG_nodes_v5.tsv')))
    logger.info("KG edge is saved to {}".format(os.path.join(args.output_dir, 'KG_edges_v5.tsv')))
    
    logger.info(f'Done!')

