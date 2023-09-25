"""
This script is used to integrate all KG2 data into a knolwedge graph.
"""

## Import standard libraries
import os
import sys
from tqdm import tqdm, trange
from glob import glob
import pandas as pd
import argparse
import itertools
import logging

## Import custom libraries
from utils import get_logger, read_tsv_file, Node, Edge, KnowledgeGraph, MappingLink

def change_prefix(synonym):
    prefix = synonym.split(':')[0]
    value = synonym.split(':')[1]
    
    if prefix == 'ICD9':
        return 'ICD-9:'+ value
    elif prefix == 'ICD10':
        return 'ICD-10:'+value
    elif prefix == 'ICD11':
        return 'ICD-11:'+value
    elif prefix == 'MESH':
        return 'MeSH:'+value
    elif prefix == 'PUBCHEM.COMPOUND':
        return 'PubChem:'+value
    elif prefix == 'CHEBI':
        return 'ChEBI:'+value
    else:
        return synonym


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Integrate all KG2 data into a knolwedge graph')
    parser.add_argument('--existing_KG_nodes', type=str, help='path of the existing knowledge graph nodes')
    parser.add_argument('--existing_KG_edges', type=str, help='path of the existing knowledge graph edges')
    parser.add_argument('--data_dir', type=str, help='path of the KG2 data directory')
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

    # Load KG2 data
    logger.info('loading KG2 node tsv file')
    # Read the header of the KG2 edge tsv file
    logger.info('reading the header of the KG2 edge tsv file')
    file_path = os.path.join(args.data_dir, "edges_c_header.tsv")
    if not os.path.exists(file_path):
        logger.error(f'Could not find the header of the KG2 edge tsv file at {file_path}')
        sys.exit(1)
    else:
        edges_c_header = read_tsv_file(file_path)[0]
    # Read the KG2 edge tsv file
    logger.info('reading the KG2 edge tsv file')
    file_path = os.path.join(args.data_dir, "edges_c.tsv")
    if not os.path.exists(file_path):
        logger.error(f'Could not find the KG2 edge tsv file at {file_path}')
        sys.exit(1)
    else:
        edges_c = read_tsv_file(file_path)
    # Read the header of the KG2 node tsv file
    logger.info('reading the header of the KG2 node tsv file')
    file_path = os.path.join(args.data_dir, "nodes_c_header.tsv")
    if not os.path.exists(file_path):
        logger.error(f'Could not find the header of the KG2 node tsv file at {file_path}')
        sys.exit(1)
    else:
        nodes_c_header = read_tsv_file(file_path)[0] 
    # Read the KG2 node tsv file
    logger.info('reading the KG2 node tsv file')
    file_path = os.path.join(args.data_dir, "nodes_c.tsv")
    if not os.path.exists(file_path):
        logger.error(f'Could not find the KG2 node tsv file at {file_path}')
        sys.exit(1)
    else:
        nodes_c = read_tsv_file(file_path)

    # Select the nodes with 'biolink:Disease' and 'biolink:PhenotypicFeature', as well as
    # with synonyms of KEGG concepts
    selected_node_categories = ['biolink:Disease', 'biolink:PhenotypicFeature']
    selected_node_synonyms = ['KEGG.COMPOUND', 'KEGG.DRUG', 'KEGG.ENZYME', 'KEGG.GLYCAN', 'KEGG.REACTION']
    node_category_index = nodes_c_header.index('category')
    node_equivalent_curies_index = nodes_c_header.index('equivalent_curies:string[]')
    selected_nodes = []
    selected_nodes = []
    for node in tqdm(nodes_c, desc="Selecting nodes based on node categories and synonyms"):
        if node[node_category_index] in selected_node_categories:
            selected_nodes.append(node)
        elif len([x for x in node[node_equivalent_curies_index].split('ǂ') if x.split(":")[0] in selected_node_synonyms]) > 0:
            selected_nodes.append(node)
    # de-duplicate selected nodes
    temp_check = {}
    for node in selected_nodes:
        if node[0] not in temp_check:
            temp_check[node[0]] = node
    selected_nodes = list(temp_check.values())

    # change synonyms prefix to match KEGG nodes
    temp_selected_nodes = []
    for node in tqdm(selected_nodes, desc="Changing synonyms prefix to match KEGG nodes"):
        node_equivalent_curies = node[node_equivalent_curies_index].split('ǂ')
        temp_synonyms = []
        for key in node_equivalent_curies:
            temp_synonyms.append(change_prefix(key))  
        node[node_equivalent_curies_index] = 'ǂ'.join([x for x in temp_synonyms if x])
        temp_selected_nodes.append(node)
    selected_nodes_c_df = pd.DataFrame(temp_selected_nodes, columns=nodes_c_header)

    # Filter edges based on the selected nodes and data sources
    edges_c_df = pd.DataFrame(edges_c, columns=edges_c_header)
    selected_indexes = edges_c_df['subject'].isin(selected_nodes_c_df['id:ID']) & edges_c_df['object'].isin(selected_nodes_c_df['id:ID'])
    selected_edge_c_df = edges_c_df.loc[selected_indexes,:].reset_index(drop=True)
    ## merge the edges with same subject, object, predicate
    temp_dict = dict()
    for row in tqdm(selected_edge_c_df.to_numpy(), desc="Merging edges with same subject, object, predicate"):
        temp_subject, temp_object, temp_predicate, primary_knowledge_source, temp_publications = row[:5]
        if (temp_subject, temp_predicate, temp_object) in temp_dict:
            existing_temp_ks = temp_dict[(temp_subject, temp_predicate, temp_object)]['knowledge_source']
            temp_dict[(temp_subject, temp_predicate, temp_object)]['knowledge_source'] = list(set(existing_temp_ks + primary_knowledge_source.split('; ')))
            temp_dict[(temp_subject, temp_predicate, temp_object)]['publications'] = list(set(temp_publications.split('ǂ') + temp_dict[(temp_subject, temp_predicate, temp_object)]['publications']))
        else:
            temp_dict[(temp_subject, temp_predicate, temp_object)] = dict()
            temp_dict[(temp_subject, temp_predicate, temp_object)]['knowledge_source'] = primary_knowledge_source.split('; ')
            temp_dict[(temp_subject, temp_predicate, temp_object)]['publications'] = list(set(temp_publications.split('ǂ')))
    selected_edge_c_df = pd.DataFrame([[x[0][0], x[0][2], x[0][1], x[1]['knowledge_source'], x[1]['publications']] for x in temp_dict.items()], columns=edges_c_header[:5])
    knowledge_source_index = edges_c_header.index('primary_knowledge_source')
    temp_selected_edges = []
    for row in tqdm(selected_edge_c_df.to_numpy(), desc="Filtering out edges that are only from SemMedDB"):
        knowledge_sources = row[knowledge_source_index]
        if not (len(knowledge_sources) == 1 and knowledge_sources[0] == 'infores:semmeddb'):
            temp_selected_edges.append((row[0], row[1], row[2], row[3], row[4]))
    selected_edge_c_df = pd.DataFrame(temp_selected_edges, columns=edges_c_header[:5])

    # Map the selected KG2 nodes to KG nodes
    kegg_nodes_matched_dict = {}
    kegg_prefix_mapping = {'KEGG.COMPOUND':'KEGG:cpd', 'KEGG.DRUG':'KEGG:dr', 'KEGG.ENZYME':'KEGG:ec', 'KEGG.GLYCAN':'KEGG:gl', 'KEGG.REACTION':'KEGG:rn'}
    node_type_mapping = {'biolink:Disease': 'Disease', 'biolink:PhenotypicFeature': 'Phenotypic_Feature', 'KEGG:dr': 'Drug', 'KEGG:cpd': 'Compound', 'KEGG:ec': 'Enzyme', 'KEGG:gl': 'Glycan', 'KEGG:rn': 'Reaction'}
    reliable_prefix = {
        'Disease': ['MONDO', 'OMIM', 'LOINC', 'RXNORM', 'DOID', 'ORPHANET', 'ICD-9', 'ICD-10','MeSH', 'UMLS'],
        'Phenotypic_Feature': ['HP', 'NBO', 'SYMP', 'PSY', 'UMLS'],
        'Drug': ['DRUGBANK', 'KEGG.DRUG', 'DrugCentral', 'VANDF', 'RXNORM'],
        'Compound': ['KEGG.COMPOUND', 'PathWhiz.Compound', 'HMDB', 'CHEMBL.COMPOUND', 'PubChem', 'ChEBI', 'RXNORM'],
        'Enzyme': ['KEGG.ENZYME', 'PathWhiz.ProteinComplex', 'UniProtKB'],
        'Glycan': ['KEGG.GLYCAN'],
        'Reaction': ['KEGG.REACTION', 'GO']
    }
    
    kg2_index = selected_nodes_c_df.columns.get_loc('id:ID')
    kg2_category_index = selected_nodes_c_df.columns.get_loc('category')
    kg2_name_index = selected_nodes_c_df.columns.get_loc('name')
    kg2_description_index = selected_nodes_c_df.columns.get_loc('description')
    kg2_iri_index = selected_nodes_c_df.columns.get_loc('iri')
    kg2node_to_synonyms = {}
    for kg2_node in tqdm(selected_nodes_c_df.to_numpy(), desc="Mapping the selected KG2 nodes to KEGG nodes"):
        equivalent_curies = kg2_node[node_equivalent_curies_index].split('ǂ')
        equivalent_curies = [kegg_prefix_mapping[x.split(":")[0]]+"_"+x.split(":")[1] if x.split(":")[0] in kegg_prefix_mapping else x for x in equivalent_curies]
        ## Find existing node in KG
        existing_nodes = list(set([kg.find_node_by_synonym(equivalent_curie) for equivalent_curie in equivalent_curies if kg.find_node_by_synonym(equivalent_curie)]))
        if len(existing_nodes) == 0:
            # No existing node found, create a new node
            # Determine node type
            temp_node_type = [node_type_mapping.get(kg2_node[kg2_category_index], None)] if node_type_mapping.get(kg2_node[kg2_category_index], None) else []
            temp_node_type += [node_type_mapping[x.split('_')[0]] for x in equivalent_curies if x.split('_')[0] in node_type_mapping]
            temp_node_type = list(set(temp_node_type))
            if len(temp_node_type) > 1:
                # Drop this KG2 node because of ambiguous node type
                continue
            else:
                temp_node_type = temp_node_type[0]
            allowable_prefixes = reliable_prefix[temp_node_type]
            temp_all_names = [kg2_node[kg2_name_index]]
            temp_synonyms = [kegg_prefix_mapping[y.split(":")[0]]+"_"+y.split(":")[1] if y.split(":")[0] in kegg_prefix_mapping else y for y in [x for x in kg2_node[node_equivalent_curies_index].split('ǂ') if x.split(":")[0] in allowable_prefixes]]
            temp_description_dict = {'RTX-KG2 Description': kg2_node[kg2_description_index]}
            temp_knowledge_source = [x.split(':')[0].upper() for x in temp_synonyms]
            temp_link = [MappingLink(x) for x in temp_synonyms if MappingLink(x)]
            if len(temp_synonyms) > 0:
                temp_node = Node(node_type=temp_node_type, all_names=temp_all_names, description=list(temp_description_dict.items()), knowledge_source=temp_knowledge_source, synonyms=temp_synonyms, link=temp_link, is_pathogen=False)
                kg.add_node(temp_node)
                kg2node_to_synonyms[kg2_node[kg2_index]] = temp_synonyms
        elif len(existing_nodes) == 1:
            ## kg has an unique id for this node
            node_id = existing_nodes[0]
            existing_node = kg.get_node_by_id(node_id)
            # Determine node type
            temp_node_type = [node_type_mapping.get(kg2_node[kg2_category_index], None)] if node_type_mapping.get(kg2_node[kg2_category_index], None) else []
            temp_node_type += [node_type_mapping[x.split('_')[0]] for x in equivalent_curies if x.split('_')[0] in node_type_mapping]
            temp_node_type = list(set(temp_node_type))
            if len(temp_node_type) > 1:
                # Drop this KG2 node because of ambiguous node type
                continue
            else:
                temp_node_type = temp_node_type[0]
            allowable_prefixes = reliable_prefix[temp_node_type]
            temp_all_names = [kg2_node[kg2_name_index]]
            temp_synonyms = [kegg_prefix_mapping[y.split(":")[0]]+"_"+y.split(":")[1] if y.split(":")[0] in kegg_prefix_mapping else y for y in [x for x in kg2_node[node_equivalent_curies_index].split('ǂ') if x.split(":")[0] in allowable_prefixes]]
            temp_description_dict = {'RTX-KG2 Description': kg2_node[kg2_description_index]}
            temp_knowledge_source = [x.split(':')[0].upper() for x in temp_synonyms]
            temp_link = [MappingLink(x) for x in temp_synonyms if MappingLink(x)]

            if len(temp_synonyms) > 0:
                existing_node.all_names = list(set(existing_node.all_names + temp_all_names))
                description_dict = dict(existing_node.description)
                description_dict.update(temp_description_dict)
                existing_node.description = list(description_dict.items())
                existing_node.knowledge_source = list(set(existing_node.knowledge_source + temp_knowledge_source))
                for synonym in temp_synonyms:
                    kg.map_synonym_to_node_id[synonym] = node_id
                existing_node.synonyms = list(set(existing_node.synonyms + temp_synonyms))
                existing_node.link = list(set(existing_node.link + temp_link))
                kg2node_to_synonyms[kg2_node[kg2_index]] = temp_synonyms
        elif len(existing_nodes) > 1:
            ## kg has multiple ids for this node, we don't add any KG2 node info to this node because of ambiguity
            pairwise = list(itertools.combinations(existing_nodes, 2))
            for pair in pairwise:
                # Add edge to the knowledge graph
                kg.add_edge(Edge(source_node=pair[0], target_node=pair[1], predicate='biolink:chemically_similar_to', knowledge_source=['KG2']))
                kg.add_edge(Edge(source_node=pair[1], target_node=pair[0], predicate='biolink:chemically_similar_to', knowledge_source=['KG2']))
                    
    # Add KG2 edges into the existing KG
    edge_count = 0
    for row in tqdm(selected_edge_c_df.to_numpy(), desc="Adding KG2 edges into the existing KG"):
        temp_subject, temp_object, temp_predicate, temp_knowledge_source, temp_publications = row
        # Skip the KG2 edge if the subject or object is not in the KG
        if temp_subject not in kg2node_to_synonyms or temp_object not in kg2node_to_synonyms:
            continue
        temp_subject_synonyms = list(set([kg.find_node_by_synonym(x) for x in kg2node_to_synonyms[temp_subject]]))
        if len(temp_subject_synonyms) == 1:
            temp_subject = temp_subject_synonyms[0]
        else:
            raise ValueError(f"KG2 subject {temp_subject} has {len(temp_object_synonyms)} ids in the KG")
        temp_object_synonyms = list(set([kg.find_node_by_synonym(x) for x in kg2node_to_synonyms[temp_object]]))
        if len(temp_object_synonyms) == 1:
            temp_object = temp_object_synonyms[0]
        else:
            raise ValueError(f"KG2 object {temp_object} has {len(temp_object_synonyms)} ids in the KG")
        temp_knowledge_source = [x.replace('infores:','').upper() for x in temp_knowledge_source]
        if temp_subject != temp_object:
            kg.add_edge(Edge(source_node=temp_subject, target_node=temp_object, predicate=temp_predicate, knowledge_source=temp_knowledge_source))
            edge_count += 1
    logger.info(f"Added {edge_count} KG2 edges into the existing KG")

    # Save the knowledge graph
    logger.info("Saving the knowledge graph...")
    kg.save_graph(save_dir = args.output_dir, node_filename = 'KG_nodes_v3.tsv', edge_filename = 'KG_edges_v3.tsv')
    logger.info("KG node is saved to {}".format(os.path.join(args.output_dir, 'KG_nodes_v3.tsv')))
    logger.info("KG edge is saved to {}".format(os.path.join(args.output_dir, 'KG_edges_v3.tsv')))
    
    logger.info(f'Done!')
