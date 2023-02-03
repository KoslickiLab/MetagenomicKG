import os, sys
project_name = os.environ['project_name']
path_list = os.path.abspath(__file__).split('/')
project_index = path_list.index(project_name)
path_to_scripts = '/'.join(path_list[:(project_index+1)] + ['scripts'])
sys.path.append(path_to_scripts)
from utility_scripts import utils, node_synonymizer
import pandas as pd
import argparse
from tqdm import tqdm
import pickle
import csv as tsv
tsv.field_size_limit(sys.maxsize)
from glob import glob

def get_kegg_iri(kegg_id):
    kegg_entry_link = 'https://www.genome.jp/entry/'

    if kegg_id.split(':')[0] != 'tax':
        return f"{kegg_entry_link}{kegg_id}"
    else:
        if temp_taxaid_mapping.get(str(kegg_id.split(':')[1])) is not None:
            return f"{kegg_entry_link}gn:{temp_taxaid_mapping.get(str(kegg_id.split(':')[1]))}"
        else:
            return f"{kegg_entry_link}gn:{str(kegg_id.split(':')[1])}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--kg", type=str, required=True, help="The full path of original kg")
    parser.add_argument("--kg_nodes", type=str, required=True, help="The full path of the node info in original kg")
    parser.add_argument("--kg2_nodes", type=str, required=True, help="The full path of KG2 'nodes.txt' file")
    parser.add_argument("--kg2_edges", type=str, required=True, help="The full path of KG2 'edges.txt' file")
    parser.add_argument("--organism_tb", type=str, required=True, help="The full path of 'organism_table.txt' file")
    parser.add_argument("--kegg_dir", type=str, required=True, help="The full path to a directory which contains KEGG information")
    parser.add_argument("--log_dir", type=str, required=False, default=str('/'.join(path_list[:(project_index+1)] + ['log_folder'])), help="the path of logfile folder")
    parser.add_argument("--log_name", type=str, required=False, default="add_KG2_info.log", help="log file name")
    parser.add_argument("--outdir", type=str, help="The output dir")
    args = parser.parse_args()

    ## set up logger
    logger = utils.get_logger(os.path.join(args.log_dir, args.log_name))

    ## set mapping
    mapping = {'cpd': 'compound', 'dg': 'dgroup', 'dr': 'drug', 'ds': 'disease', 'ec': 'enzyme', 'gl': 'glycan', 'ko': 'orthology', 'md': 'module', 'ne': 'network', 'path': 'pathway', 'rc': 'rclass', 'rn': 'reaction', 'tax': 'taxonomy'}

    if os.path.exists(args.kg):
        kg = pd.read_csv(args.kg, sep='\t', header=0)
    else:
        logger.error(f"Not Found {args.kg}")
        exit()

    if os.path.exists(args.kg_nodes):
        kg_nodes = pd.read_csv(args.kg_nodes, sep='\t', header=0)
    else:
        logger.error(f"Not Found {args.kg_nodes}")
        exit()

    if os.path.exists(args.kg2_nodes):
        kg2_nodes = pd.read_csv(args.kg2_nodes, sep='\t', header=0)
    else:
        logger.error(f"Not Found {args.kg2_nodes}")
        exit()

    if os.path.exists(args.kg2_edges):
        kg2_edges = pd.read_csv(args.kg2_edges, sep='\t', header=0)
    else:
        logger.error(f"Not Found {args.kg2_edges}")
        exit()

    if os.path.exists(args.organism_tb):
        organism_tb = pd.read_csv(args.organism_tb, sep='\t', header=0)
        temp_taxaid_mapping = {str(row[4]):row[0] for row in organism_tb.to_numpy() if row[4] is not None}
    else:
        logger.error(f"Not Found {args.organism_tb}")
        exit()

    ## read taxonomy information
    logger.info('read taxonomy information')
    kg_nodes = pd.concat([kg_nodes,pd.DataFrame([(mapping[x.split(':')[0]],get_kegg_iri(x)) for x in kg_nodes['node_KG_ids']])], axis=1)
    kg_nodes.columns = list(kg_nodes.columns)[:-2] + ['node_type'] + ['iri']
    kg2_to_kegg = {}
    for kegg_id, kg2_ids in kg_nodes.loc[~kg_nodes['preferred_KG2_ids'].isna(),['node_KG_ids','preferred_KG2_ids']].to_numpy():
        kg2_id_list = kg2_ids.split('##')
        for kg2_id in kg2_id_list:
            if kg2_id in kg2_to_kegg:
                kg2_to_kegg[kg2_id].update([kegg_id])
            else:
                kg2_to_kegg[kg2_id] = set([kegg_id])

    ## process kg2 data
    logger.info('process kg2 data')
    nodesynonymizer = node_synonymizer.NodeSynonymizer()
    if not os.path.exists(os.path.join(args.outdir,'kg2node_mapping.pkl')):
        kg2node_mapping = {}
        kg2_nodes_colname = list(kg2_nodes.columns)
        for row in tqdm(kg2_nodes.to_numpy()):
            kg_node_id = row[kg2_nodes_colname.index('id')]
            normalizer_res = nodesynonymizer.get_canonical_curies(kg_node_id)[kg_node_id]
            if normalizer_res is None:
                preferred_category = None
            else:
                preferred_category = normalizer_res['preferred_category']
            kg2node_mapping[kg_node_id] = list(row[:4]) + list(row[5:]) + [preferred_category]
        kg2_nodes_colname = kg2_nodes_colname[:4] + kg2_nodes_colname[5:] + ['preferred_category']
        with open(os.path.join(args.outdir,'kg2node_mapping.pkl'),'wb') as file_out:
            pickle.dump([kg2node_mapping, kg2_nodes_colname], file_out)
    else:
        with open(os.path.join(args.outdir,'kg2node_mapping.pkl'),'rb') as file_in:
            kg2node_mapping, kg2_nodes_colname = pickle.load(file_in)    

    ## add new columns to kg
    logger.info('add new columns to kg')
    kg['predicate'] = 'KEGG_associated_with'
    kg['knowledge_source'] = 'KEGG_api'
    kg['publications'] = None
    kg['publications_info'] = None
    kg['kg2_ids'] = None
    check_edges = {}
    for row in tqdm(kg.to_numpy()):
        if (row[0], row[1]) in check_edges:
            check_edges[(row[0], row[1])].update(set([row[2]]))
        else:
            check_edges[(row[0], row[1])] = set([row[2]])

    ## add kg2 edges and kg2 phenotype features
    logger.info('add kg2 edges and kg2 phenotype features')
    existing_ids = list(kg2_to_kegg.keys())
    kg2_edges_colname = list(kg2_edges.columns)
    new_nodes = set()
    new_edges = []
    for row1 in tqdm(kg2_edges.to_numpy()):
        ## ignore the existing kegg-based edges in kg2
        if 'infores:kegg' in row1[2]:
            continue
        try:
            source, predicate, target = row1[kg2_edges_colname.index('triple')].split('--')
        except:
            continue
        ## ignore the 'biolink:close_match' edges in kg2
        if predicate == 'biolink:close_match':
            continue
        if len(set([source, target]).intersection(set(existing_ids))) == 0:
            continue
        else:
            if (source in existing_ids) and (target in existing_ids):
                for row2 in [tuple([x,y]) for x in kg2_to_kegg[source] for y in kg2_to_kegg[target]]:
                    if (row2[0].split(':')[0] != row2[1].split(':')[0]) and (predicate == 'biolink:subclass_of'):
                        continue
                    if (row2[0], row2[1]) in check_edges:
                        if predicate in check_edges[(row2[0], row2[1])]:
                            continue
                        else:
                            check_edges[(row2[0], row2[1])].update(set([predicate]))
                            new_edges += [tuple(list(row2)+[predicate]+list(row1[2:]))]
                    else:
                        check_edges[(row2[0], row2[1])] = set([predicate])
                        new_edges += [tuple(list(row2)+[predicate]+list(row1[2:]))]
            elif source in existing_ids:
                if kg2node_mapping[target][kg2_nodes_colname.index('preferred_category')] == 'biolink:PhenotypicFeature':
                    iri_link_list = list(kg2_nodes.loc[kg2_nodes['id'] == target,'iri'])
                    iri_link = iri_link_list[0] if len(iri_link_list) != 0 else None
                    for x in kg2_to_kegg[source]:
                        if (x, target) in check_edges:
                            if predicate in check_edges[(x, target)]:
                                continue
                            else:
                                check_edges[(x, target)].update(set([predicate]))
                                normalizer_res = nodesynonymizer.get_canonical_curies(target)[target]
                                if normalizer_res is not None:
                                    new_nodes.update([(target,target,normalizer_res['preferred_curie'],normalizer_res['preferred_name'],normalizer_res['preferred_category'],'PhenotypicFeature'.lower(),iri_link)])
                                else:
                                    new_nodes.update([(target,target,None,None,None,'PhenotypicFeature'.lower(),iri_link)])
                                new_edges += [tuple([x]+[target]+[predicate]+list(row1[2:]))]
                        else:
                            check_edges[(x, target)] = set([predicate])
                            normalizer_res = nodesynonymizer.get_canonical_curies(target)[target]
                            if normalizer_res is not None:
                                new_nodes.update([(target,target,normalizer_res['preferred_curie'],normalizer_res['preferred_name'],normalizer_res['preferred_category'],'PhenotypicFeature'.lower(),iri_link)])
                            else:
                                new_nodes.update([(target,target,None,None,None,'PhenotypicFeature'.lower(),iri_link)])
                            new_edges += [tuple([x]+[target]+[predicate]+list(row1[2:]))]
            elif target in existing_ids:
                if kg2node_mapping[source][kg2_nodes_colname.index('preferred_category')] == 'biolink:PhenotypicFeature':
                    iri_link_list = list(kg2_nodes.loc[kg2_nodes['id'] == source,'iri'])
                    iri_link = iri_link_list[0] if len(iri_link_list) != 0 else None
                    for x in kg2_to_kegg[target]:
                        if (source, x) in check_edges:
                            if predicate in check_edges[(source, x)]:
                                continue
                            else:
                                check_edges[(source, x)].update(set([predicate]))
                                normalizer_res = nodesynonymizer.get_canonical_curies(source)[source]
                                if normalizer_res is not None:
                                    new_nodes.update([(source,source,normalizer_res['preferred_curie'],normalizer_res['preferred_name'],normalizer_res['preferred_category'],'PhenotypicFeature'.lower(),iri_link)])
                                else:
                                    new_nodes.update([(source,source,None,None,None,'PhenotypicFeature'.lower(),iri_link)])
                                new_edges += [tuple([source]+[x]+[predicate]+list(row1[2:]))]
                        else:
                            check_edges[(source, x)] = set([predicate])
                            normalizer_res = nodesynonymizer.get_canonical_curies(source)[source]
                            if normalizer_res is not None:
                                new_nodes.update([(source,source,normalizer_res['preferred_curie'],normalizer_res['preferred_name'],normalizer_res['preferred_category'],'PhenotypicFeature'.lower(),iri_link)])
                            else:
                                new_nodes.update([(source,source,None,None,None,'PhenotypicFeature'.lower(),iri_link)])
                            new_edges += [tuple([source]+[x]+[predicate]+list(row1[2:]))]

    ## merge KG2 info and the existing KEGG-based KG
    new_nodes = pd.DataFrame(new_nodes)
    new_nodes.columns = kg_nodes.columns
    kg_nodes.loc[kg_nodes['node_KG_ids'].isna(),'node_KG_ids'] = None
    kg_nodes.loc[kg_nodes['KG2_ids'].isna(),'KG2_ids'] = None
    kg_nodes.loc[kg_nodes['preferred_KG2_ids'].isna(),'preferred_KG2_ids'] = None
    kg_nodes.loc[kg_nodes['preferred_KG2_names'].isna(),'preferred_KG2_names'] = None
    kg_nodes.loc[kg_nodes['preferred_KG2_node_types'].isna(),'preferred_KG2_node_types'] = None
    kg_nodes.loc[kg_nodes['node_type'].isna(),'node_type'] = None
    kg_nodes.loc[kg_nodes['iri'].isna(),'iri'] = None 
    new_nodes = pd.concat([kg_nodes,new_nodes]).reset_index(drop=True)

    ## add node descriptions
    id_to_description = {}
    id_to_description.update({row[list(kg2_nodes.columns).index('id')]:row[list(kg2_nodes.columns).index('description')] for row in kg2_nodes.to_numpy()})
    id_to_description.update({row[0]:row[1] for path in glob(os.path.join(kegg_dir, 'kegg_*.txt')) for row in pd.read_csv(path, sep='\t', header=0).to_numpy()})

    for index in range(len(new_nodes)):
        new_nodes.loc[index, 'description'] = id_to_description.get(new_nodes.loc[index, 'id:ID'])

    ## remove all kg2 edges whose source is from semmeddb only 
    new_edges = [row for row in new_edges if row[list(kg.columns).index('knowledge_source:string[]')]!='infores:semmeddb']
    ## remove self-edges from KG2
    new_edges = [row for row in new_edges if row[list(kg.columns).index('source')]!=row[list(kg.columns).index('target')]]


    new_edges = pd.DataFrame(new_edges)
    new_edges.columns = kg.columns
    new_edges.loc[new_edges['source'].isna(),'source'] = None
    new_edges.loc[new_edges['target'].isna(),'target'] = None
    new_edges.loc[new_edges['predicate'].isna(),'predicate'] = None
    new_edges.loc[new_edges['knowledge_source'].isna(),'knowledge_source'] = None
    new_edges.loc[new_edges['publications'].isna(),'publications'] = None
    new_edges.loc[new_edges['publications_info'].isna(),'publications_info'] = None
    new_edges.loc[new_edges['kg2_ids'].isna(),'kg2_ids'] = None  
    new_edges = pd.concat([kg,new_edges]).reset_index(drop=True)

    ## output nodes to tsv format
    new_nodes.columns = ['node_id:ID','KG2_ids','preferred_KG2_ids','preferred_KG2_names','preferred_KG2_node_types','node_type','iri','description']
    node_col_tsvfile_out = open(os.path.join(args.outdir, 'nodes_header.tsv'), 'w+')
    node_col_tsvwrite = tsv.writer(node_col_tsvfile_out, delimiter="\t",  quoting=tsv.QUOTE_MINIMAL)
    node_col_tsvwrite.writerows([list(new_nodes.columns)])
    node_col_tsvfile_out.close()
    node_tsvfile_out = open(os.path.join(args.outdir, 'nodes.tsv'), 'w+')
    node_tsvwrite = tsv.writer(node_tsvfile_out, delimiter="\t",  quoting=tsv.QUOTE_MINIMAL)
    rows = [list(row) for row in new_nodes.to_numpy()]
    node_tsvwrite.writerows(rows)
    node_tsvfile_out.close()

    ## output edges to tsv format
    new_edges.columns = ['source','target','predicate','knowledge_source:string[]','publications:string[]','publications_info','kg2_ids:string[]']
    new_edges['id'] = list(range(1,len(new_edges)+1))
    new_edges[':TYPE'] = new_edges['predicate']
    new_edges[':START_ID'] = new_edges['source']
    new_edges[':END_ID'] = new_edges['target']
    edge_col_tsvfile_out = open(os.path.join(args.outdir, 'edges_header.tsv'), 'w+')
    edge_col_tsvwrite = tsv.writer(edge_col_tsvfile_out, delimiter="\t",  quoting=tsv.QUOTE_MINIMAL)
    edge_col_tsvwrite.writerows([list(new_edges.columns)])
    edge_col_tsvfile_out.close()
    edge_tsvfile_out = open(os.path.join(args.outdir, 'edges.tsv'), 'w+')
    edge_tsvwrite = tsv.writer(edge_tsvfile_out, delimiter="\t",  quoting=tsv.QUOTE_MINIMAL)
    rows = [list(row) for row in new_edges.to_numpy()]
    edge_tsvwrite.writerows(rows)
    edge_tsvfile_out.close()