"""
Somme Common Functions

This script provides a series of functions for common use.

Knowledge Graph Structure:
Node Attributes:
node_id:ID[str]    node_type[str]    all_names[list]    description[list of tuple]    knowledge_source[list]    link[list]    synonyms[list]    is_pathogen[bool]

Edge Attributes:
source_node[str]    target_node[str]    predicate[str]    knowledge_source[list]

"""

# Import Python libraries
import os
import sys
import logging
import pandas as pd
import csv
from tqdm import tqdm, trange
from typing import List, Dict, Tuple, Union, Any, Optional, Set
csv.field_size_limit(100000000)  # set the maximum field size limit to 100MB or any value you need
import requests
import json
import re
import ast

def get_logger():
    """
    Setup a logger object
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s  [%(levelname)s]  %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger

def read_tsv_file(file_path, encoding='utf-8'):
    """
    Read a tsv file and return a list of lists
    """
    
    with open(file_path, newline='', encoding=encoding) as f:
        reader = csv.reader(f, delimiter='\t')
        data = [row for row in reader]
    return data

def check_files(file_path: str, logger):
    """
    Check if file path exists.
    :param file_path: path of the file
    :param logger: logger object
    :return: True if the file exists, False otherwise
    """
    
    if not os.path.exists(file_path):
        logger.error(f"File {file_path} does not exist.")
        return False
    else:
        return True


class Node:
    def __init__(self, node_type: str, all_names: List[str] = [], description: List[Tuple] = [], knowledge_source: List[str] = [], link: List[str] = [], synonyms: List[str] = [], is_pathogen: bool = False):
        self.id = None
        self.node_type = node_type
        self.all_names = list(set(all_names))
        self.description = list(set(description))
        self.knowledge_source = list(set(knowledge_source))
        self.link = list(set(link))
        self.synonyms = list(set(synonyms))
        self.is_pathogen = is_pathogen

class Edge:
    def __init__(self, source_node: str, target_node: str, predicate: str, description: List[Tuple] = [], knowledge_source: List[str] = []):
        self.edge_id = None
        self.source_node = source_node
        self.target_node = target_node
        self.predicate = predicate
        self.description = list(set(description))
        self.knowledge_source = list(set(knowledge_source))


class KnowledgeGraph:
    def __init__(self, logger):
        self.logger = logger
        self.nodes = {}
        self.edges = {}
        self.node_type_count = {'Microbe': 0, 'Pathway': 0, 'KO': 0, 'Network': 0, 'Disease': 0 , 'Drug': 0, 'Drug_Group': 0, 'Module': 0, 'Compound': 0, 'Enzyme': 0, 'Glycan': 0, 'Reaction': 0, 'Phenotypic_Feature': 0, 'AMR': 0}
        self.mapping_nodetype_biolink = {'Microbe': 'biolink:OrganismTaxon', 'Pathway': 'biolink:Pathway', 'KO': 'biolink:BiologicalEntity', 'Network': 'biolink:NamedThing', 'Disease': 'biolink:Disease', 'Drug': 'biolink:Drug', 'Drug_Group': 'biolink:MolecularMixture', 'Module': 'biolink:BiologicalProcess', 'Compound': 'biolink:MolecularEntity', 'Enzyme': 'biolink:Polypeptide', 'Glycan': 'biolink:MacromolecularComplex', 'Reaction': 'biolink:MolecularActivity', 'Phenotypic_Feature': 'biolink:PhenotypicFeature', 'AMR': 'biolink:Protein'}
        self.all_node_type = list(self.mapping_nodetype_biolink.keys()) + list(self.mapping_nodetype_biolink.values())
        self.mapping_biolink_nodetype = {v: k for k, v in self.mapping_nodetype_biolink.items()}
        self.map_synonym_to_node_id = {}
        self.in_edge = {}
        self.out_edge = {}
        self.node_by_type = {}
        
    def add_node(self, node: Node):
        if not isinstance(node, Node):
            self.logger.error("Node must be an instance of the Node class.")
            return

        ## Check if the node has existed
        temp_synonyms = [synonym for synonym in node.synonyms if synonym in self.map_synonym_to_node_id]
        
        
        if len(temp_synonyms) == 0:
            # This is a new node
            pass
        elif len(temp_synonyms) == len(node.synonyms):
            if len(set([self.map_synonym_to_node_id[synonym] for synonym in temp_synonyms])) == 1:
                # This node has existed
                temp_node_id = self.map_synonym_to_node_id[node.synonyms[0]]
                existing_node = self.nodes[temp_node_id]
                existing_node.all_names = list(set(existing_node.all_names + node.all_names))
                existing_node.description = list(set(existing_node.description + node.description))
                existing_node.knowledge_source = list(set(existing_node.knowledge_source + node.knowledge_source))
                existing_node.link = list(set(existing_node.link + node.link))
                existing_node.synonyms = list(set(existing_node.synonyms + node.synonyms))
                existing_node.is_pathogen = existing_node.is_pathogen or node.is_pathogen
                return
            else:
                # Some sysnonyms have conflicted node ids
                raise ValueError("Some sysnonyms have conflicted node ids!")
        else:
            if len(set([self.map_synonym_to_node_id[synonym] for synonym in temp_synonyms])) == 1:
                # This node has existed
                temp_node_id = self.map_synonym_to_node_id[temp_synonyms[0]]
                existing_node = self.nodes[temp_node_id]
                existing_node.all_names = list(set(existing_node.all_names + node.all_names))
                existing_node.description = list(set(existing_node.description + node.description))
                existing_node.knowledge_source = list(set(existing_node.knowledge_source + node.knowledge_source))
                existing_node.link = list(set(existing_node.link + node.link))
                for synonym in node.synonyms:
                    self.map_synonym_to_node_id[synonym] = temp_node_id
                existing_node.synonyms = list(set(existing_node.synonyms + node.synonyms))
                existing_node.is_pathogen = existing_node.is_pathogen or node.is_pathogen
                return
            else:
                # Some sysnonyms have conflicted node ids
                raise ValueError("Some sysnonyms have conflicted node ids!")

        if node.node_type in self.node_type_count:
            self.node_type_count[node.node_type] += 1
            temp_node_id = node.node_type + ":" + str(self.node_type_count[node.node_type])
            if node.id is None:
                node.id = temp_node_id
            elif node.id.split(':')[0] != node.node_type:
                self.logger.warning("Node type {} does not match node id {}! Set a new node id: {}".format(node.node_type, node.id, temp_node_id))
                node.id = temp_node_id
        else:
            self.logger.warning("Node type {} is not supported!".format(node.node_type))
            return

        ## Add node synonyms to the map
        for synonym in node.synonyms:
            self.map_synonym_to_node_id[node.id] = node.id
            self.map_synonym_to_node_id[synonym] = node.id
            
            # if synonym not in self.map_synonym_to_node_id:

            # elif self.map_synonym_to_node_id[synonym] != node.id:
            #     self.logger.warning("Synonym {} is mapped to multiple nodes: {} and {}!".format(synonym, self.map_synonym_to_node_id[synonym], node.id))
            #     self.logger.warning("Remove this synonym from both nodes!")
            #     existing_node = self.nodes[self.map_synonym_to_node_id[synonym]]
            #     if synonym in existing_node.synonyms:
            #         existing_node.synonyms.remove(synonym)
            #     if len(existing_node.synonyms) == 0:
            #         raise ValueError("Node {} has no synonyms!".format(existing_node.id))
            #     if synonym in node.synonyms:
            #         node.synonyms.remove(synonym)
            #     if len(node.synonyms) == 0:
            #         raise ValueError("Node {} has no synonyms!".format(node.id))
            #     del self.map_synonym_to_node_id[synonym]

        # Add node to the graph
        node.node_type = self.mapping_nodetype_biolink[node.node_type]
        if node.node_type not in self.node_by_type:
            self.node_by_type[node.node_type] = []
        self.node_by_type[node.node_type] += [node.id]
        self.nodes[node.id] = node

    def add_edge(self, edge: Edge):
        if not isinstance(edge, Edge):
            self.logger.warning("Edge must be an instance of the Edge class.")
            return
        
        # Normalize the edge id
        edge_source_node = self.find_node_by_synonym(edge.source_node)
        edge_target_node = self.find_node_by_synonym(edge.target_node)
        if edge_source_node is None:
            self.logger.warning("{} is not in the graph!".format(edge.source_node))
            return
        if edge_target_node is None:
            self.logger.warning("{} is not in the graph!".format(edge.target_node))
            return
        if edge.predicate is None:
            self.logger.warning("Predicate cannot be None!")
            return
        
        edge.source_node = edge_source_node
        edge.target_node = edge_target_node
        edge.edge_id = edge.source_node + "_" + edge.predicate + "_" + edge.target_node
        
        # Add edge to the graph
        if edge.edge_id not in self.edges:
            self.edges[edge.edge_id] = edge
        else:
            # self.edges[edge.edge_id].predicate = list(set(self.edges[edge.edge_id].predicate + edge.predicate))
            temp_description_dict = dict(edge.description)
            old_temp_description_dict = dict(self.edges[edge.edge_id].description)
            for key in temp_description_dict:
                if key in old_temp_description_dict:
                    old_temp_description_dict[key] = '#####'.join(list(set([temp_description_dict[key]] + old_temp_description_dict[key].split('#####'))))
                else:
                    old_temp_description_dict[key] = temp_description_dict[key]
            self.edges[edge.edge_id].description = list(old_temp_description_dict.items())
            self.edges[edge.edge_id].knowledge_source = list(set(self.edges[edge.edge_id].knowledge_source + edge.knowledge_source))
        
        # Add edge to the in_edge and out_edge
        if edge.target_node not in self.in_edge:
            self.in_edge[edge.target_node] = set()
        self.in_edge[edge.target_node].add(edge.source_node)
        if edge.source_node not in self.out_edge:
            self.out_edge[edge.source_node] = set()
        self.out_edge[edge.source_node].add(edge.target_node)

    def get_node_by_type(self, node_type):
        if node_type not in self.all_node_type:
            self.logger.warning("Node type {} is not supported!".format(node_type))
            return []
        else:
            return self.node_by_type[node_type]
        

    def get_node_by_id(self, node_id):
        return self.nodes.get(self.find_node_by_synonym(node_id), None)

    def get_edge_by_id(self, edge_id):
        return self.edges.get(edge_id, None)

    def find_node_by_synonym(self, synonym):
        
        ## Consider GCF_XXXXX and GCA_XXXXX are the same
        synonyms = list(set([synonym] + [synonym.replace(":GCF_", ":GCA_"), synonym.replace(":GCA_", ":GCF_")]))
        
        result = [self.map_synonym_to_node_id.get(synonym, None) for synonym in synonyms if self.map_synonym_to_node_id.get(synonym, None)]
        if len(result) == 0:
            return None
        else:
            return result[0]
    
    def find_all_in_edges(self, node_id):
        node_id = self.find_node_by_synonym(node_id)
        in_node_id_list = self.in_edge.get(node_id, None)
        if in_node_id_list is None:
            return []
        else:
            return [self.get_edge_by_id(in_node_id + "_" + node_id) for in_node_id in in_node_id_list]
    
    def find_all_out_edges(self, node_id):
        node_id = self.find_node_by_synonym(node_id)
        out_node_id_list = self.out_edge.get(node_id, None)
        if out_node_id_list is None:
            return []
        else:
            return [self.get_edge_by_id(node_id + "_" + out_node_id) for out_node_id in out_node_id_list]

    def count_nodes(self):
        return len(self.nodes)

    def count_edges(self):
        return len(self.edges)

    def nodes_to_dataframe(self):
        data = []
        for node in tqdm(self.nodes.values(), desc="Converting nodes to dataframe"):
            data.append({
                "node_id": node.id,
                "node_type": node.node_type,
                "all_names": node.all_names,
                "description": node.description,
                "knowledge_source": node.knowledge_source,
                "link": node.link,
                "synonyms": node.synonyms,
                "is_pathogen": node.is_pathogen
            })
        return pd.DataFrame(data)

    def edges_to_dataframe(self):
        data = []
        for edge in tqdm(self.edges.values(), desc="Converting edges to dataframe"):
            data.append({
                "source_node": edge.source_node,
                "target_node": edge.target_node,
                "predicate": edge.predicate,
                "description": edge.description,
                "knowledge_source": edge.knowledge_source
            })
        return pd.DataFrame(data)

    def save_graph(self, save_dir: str, node_filename: str = 'KG_node.tsv', edge_filename: str = 'KG_edge.tsv'):
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        ## Save the graph to tsv files
        savepath_node_filename = os.path.join(save_dir, node_filename)
        self.nodes_to_dataframe().to_csv(savepath_node_filename, sep='\t', index=False)
        self.logger.info("Save node information to {}".format(savepath_node_filename))
        savepath_edge_filename = os.path.join(save_dir, edge_filename)
        self.edges_to_dataframe().to_csv(savepath_edge_filename, sep='\t', index=False)
        self.logger.info("Save edge information to {}".format(savepath_edge_filename))
        
    def load_graph(self, load_dir: str, node_filename: str = 'KG_node.tsv', has_node_header = True, edge_filename: str = 'KG_edge.tsv', has_edge_header = True):
        ## Load the graph from tsv files
        loadpath_node_filename = os.path.join(load_dir, node_filename)
        self.logger.info("Load node information from {}".format(loadpath_node_filename))
        node_data = read_tsv_file(loadpath_node_filename)
        if has_node_header:
            node_df = pd.DataFrame(node_data[1:], columns=node_data[0])
        else:
            node_df = pd.DataFrame(node_data)
        # convert node attribute to correct type
        node_df['all_names'] = node_df['all_names'].apply(lambda x: eval(x) if x != '[nan]' else [])
        node_df['description'] = node_df['description'].apply(lambda x: eval(x))
        node_df['knowledge_source'] = node_df['knowledge_source'].apply(lambda x: eval(x))
        node_df['link'] = node_df['link'].apply(lambda x: eval(x))
        node_df['synonyms'] = node_df['synonyms'].apply(lambda x: eval(x))
        node_df['is_pathogen'] = node_df['is_pathogen'].apply(lambda x: eval(x))
        
        loadpath_edge_filename = os.path.join(load_dir, edge_filename)
        self.logger.info("Load edge information from {}".format(loadpath_edge_filename))
        edge_data = read_tsv_file(loadpath_edge_filename)
        if has_edge_header:
            edge_df = pd.DataFrame(edge_data[1:], columns=edge_data[0])
        else:
            edge_df = pd.DataFrame(edge_data)
        # convert edge attribute to correct type
        # edge_df['predicate'] = edge_df['predicate'].apply(lambda x: eval(x))
        if 'description' in edge_df.columns:
            edge_df['description'] = edge_df['description'].apply(lambda x: eval(x))
        else:
            edge_df['description'] = [[]] * edge_df.shape[0]
        edge_df['knowledge_source'] = edge_df['knowledge_source'].apply(lambda x: eval(x))

        ## Add nodes and edges to the graph
        for _, row in tqdm(node_df.iterrows(), desc="Add nodes to the graph", total=node_df.shape[0]):
            node = Node(self.mapping_biolink_nodetype[row['node_type']], row['all_names'], row['description'], row['knowledge_source'], row['link'], row['synonyms'], row['is_pathogen'])
            node.id = row['node_id']
            self.add_node(node)

        for _, row in tqdm(edge_df.iterrows(), desc="Add edges to the graph", total=edge_df.shape[0]):
            edge = Edge(row['source_node'], row['target_node'], row['predicate'], row['knowledge_source'])
            self.add_edge(edge)

        self.logger.info("Load graph successfully!")
        
        
def change_prefix(synonym):
    prefix = synonym.split(':')[0]
    value = synonym.split(':')[1]
    
    if prefix not in ['MONDO', 'OMIM', 'LOINC', 'RXNORM', 'DOID', 'ORPHANET', 'ICD-9', 'ICD-10', 'MeSH', 'UMLS']:
        return None
    
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


def extract_disease_synonyms(name, nodesynonymizer, umls_class, umls_api_key=None):

    
    # Use Node Synonymizer to get canonical curies
    normliazer = nodesynonymizer.get_canonical_curies(names=name)[name]
    if normliazer:
        canonical_curie = normliazer['preferred_curie']
        # Get all synonyms
        synonyms = nodesynonymizer.get_equivalent_nodes(canonical_curie)[canonical_curie]
        if synonyms:
            synonyms = list(set([change_prefix(synonym) for synonym in synonyms if change_prefix(synonym)]))
            return synonyms
        elif change_prefix(canonical_curie):
            return [change_prefix(canonical_curie)]
        else:
            return []
    else:
        if umls_api_key:
            # Use UMLS API to get canonical curies
            temp_curie = umls_class.get_mapping(name)
            if not temp_curie:
                return []
            else:
                normliazer = nodesynonymizer.get_canonical_curies(temp_curie)[temp_curie]
                if normliazer:
                    canonical_curie = normliazer['preferred_curie']
                    # Get all synonyms
                    synonyms = nodesynonymizer.get_equivalent_nodes(canonical_curie)[canonical_curie]
                    if synonyms:
                        synonyms = list(set([change_prefix(synonym) for synonym in synonyms if change_prefix(synonym)]))
                        return synonyms
                    elif change_prefix(canonical_curie):
                        return [change_prefix(canonical_curie)]
                    else:
                        return []
                else:
                    return [change_prefix(temp_curie)]
        else:
            return []


class UMLSMapping:
    
    def __init__(self, apikey):
        self.umls_api = "https://uts-ws.nlm.nih.gov/rest/search/current"
        self.query = {'string': None,'apiKey':apikey, 'pageNumber':0}

    def _parse_response(self, response_json):
        items = response_json['result']['results']
        if len(items) == 0:
            return None
        else:
            return f"UMLS:{items[0]['ui']}"
        
    def get_mapping(self, name = None):
        
        if not isinstance(name, str):
            print("Valid input should be a str or list", flush=True)
            return None
        
        self.query['string'] = name
        response = requests.get(self.umls_api, params=self.query)
        if response.status_code == 200:
            return self._parse_response(response.json())
        else:
            return None

class OxOMapping:

    def __init__(self):
        self.oxo_api = 'https://www.ebi.ac.uk/spot/oxo/api/search'
        self.data = {"ids": None, "distance": 1}
        self.headers = {'content-type': 'application/json'}
        
    def _parse_response(self, response_json, output_only = None):
        mapping_results = response_json['_embedded']['searchResults']
        res_dict = {}
        for result in mapping_results:
            if not output_only:
                res_dict[result['queryId']] = [x['curie'] for x in result['mappingResponseList']]
            else:
                res_dict[result['queryId']] = [x['curie'] for x in result['mappingResponseList'] if x['distance'] == output_only]
                
        return res_dict
    
    def get_mappings(self, ids = None, mappingTarget = None, distance = 1, output_only = None):
        
        if not ids:
            return None
        elif isinstance(ids, str):
            ids = [ids]
        elif isinstance(ids, list):
            ids = ids
        else:
            print("Valid input should be a str or list", flush=True)
            return None
        
        if output_only:
            if not (output_only == 1 or output_only == 2 or output_only == 3):
                print("The parameter 'output_only' should be 1 or 2 or 3", flush=True)
        
        self.data['ids'] = ids
        self.data['distance'] = distance
        if isinstance(mappingTarget, list):
            self.data['mappingTarget'] = mappingTarget
        response = requests.post(self.oxo_api, data=json.dumps(self.data), headers=self.headers)
        if response.status_code == 200:
            return self._parse_response(response.json(), output_only)
        else:
            return None
        
def MappingLink(synonym):
    prefix = synonym.split(':')[0]
    if prefix in ['MONDO', 'OMIM', 'LOINC', 'RXNORM', 'DOID', 'ORPHANET', 'ICD-9', 'ICD-10', 'MeSH', 'UMLS', 'HP', 'NBO', 'SYMP', 'PSY', 'DRUGBANK', 'KEGG', 'DrugCentral', 'VANDF', 'PathWhiz.Compound', 'HMDB', 'CHEMBL.COMPOUND', 'PubChem', 'ChEBI', 'PathWhiz.ProteinComplex', 'UniProtKB', 'GO']:
        if prefix == 'MONDO':
            return f"http://purl.obolibrary.org/obo/{synonym.replace(':', '_')}"
        elif prefix == 'OMIM':
            return f"https://www.omim.org/entry/{synonym.replace('OMIM:', '_')}"
        elif prefix == 'LOINC':
            return None
        elif prefix == 'RXNORM':
            return f"https://bioportal.bioontology.org/ontologies/{prefix}?p=classes&conceptid={synonym.split(':')[1]}"
        elif prefix == 'DOID':
            return f"https://www.ebi.ac.uk/ols/ontologies/doid/terms?obo_id={synonym}/"
        elif prefix == 'ORPHANET':
            return f"https://identifiers.org/{synonym}"
        elif prefix == 'ICD-9':
            return None
        elif prefix == 'ICD-10':
            return None        
        elif prefix == 'MeSH':
            return f"https://id.nlm.nih.gov/mesh/{synonym.split(':')[1]}.html"
        elif prefix == 'HP':
            return f"https://hpo.jax.org/app/browse/term/{synonym}"
        elif prefix == 'NBO':
            return None
        elif prefix == 'SYMP':
            return None
        elif prefix == 'PSY':
            return None
        elif prefix == 'UMLS':
            return f"http://linkedlifedata.com/resource/umls/id/{synonym.split(':')[1]}"
        elif prefix == 'DRUGBANK':
            return f"https://go.drugbank.com/drugs/{synonym.split(':')[1]}"
        elif prefix == 'KEGG':
            return f"https://www.kegg.jp/entry/{synonym.split(':')[1].split('_')[1]}"
        elif prefix == 'DrugCentral':
            return f"https://drugcentral.org/drugcard/{synonym.split(':')[1]}"
        elif prefix == 'VANDF':
            return f"https://bioportal.bioontology.org/ontologies/VANDF?p=classes&conceptid={synonym.split(':')[1]}"
        elif prefix == 'PathWhiz.Compound':
            return None
        elif prefix == 'HMDB':
            return f"https://hmdb.ca/metabolites/{synonym.split(':')[1]}"
        elif prefix == 'CHEMBL.COMPOUND':
            return f"https://www.ebi.ac.uk/chembl/compound_report_card/{synonym.split(':')[1]}/"
        elif prefix == 'PubChem':
            return f"https://pubchem.ncbi.nlm.nih.gov/compound/{synonym.split(':')[1]}"
        elif prefix == 'ChEBI':
            return f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={synonym}"
        elif prefix == 'PathWhiz.ProteinComplex':
            return None
        elif prefix == 'UniProtKB':
            return f"https://www.uniprot.org/uniprotkb/{synonym.split(':')[1]}/entry"
        elif prefix == 'GO':
            return f"https://www.salivaryproteome.org/public/index.php/Special:Ontology_Term/{synonym}"
    else:
        return None