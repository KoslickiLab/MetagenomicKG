import os, sys
import argparse
project_name = os.environ['project_name']
path_list = os.path.abspath(__file__).split('/')
project_index = path_list.index(project_name)
path_to_scripts = '/'.join(path_list[:(project_index+1)] + ['scripts'])
sys.path.append(path_to_scripts)
from utility_scripts import utils
import pandas as pd
from glob import glob
from tqdm import tqdm
import pickle
import csv as tsv
tsv.field_size_limit(sys.maxsize)
from Bio import SeqIO
import requests

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--kg_nodes_header", type=str, required=True, help="The full path of the node header")
    parser.add_argument("--kg_nodes", type=str, required=True, help="The full path of the node info in kg")
    parser.add_argument("--kg_edges_header", type=str, required=True, help="The full path of the edge header")
    parser.add_argument("--kg_edges", type=str, required=True, help="The full path of the edge info in kg")
    parser.add_argument("--ko_hierarchy", type=str, required=True, help="The full path of ko hierarchy information")
    parser.add_argument("--ko_gene_list", type=str, required=True, help="The full path of ko gene list")
    parser.add_argument("--PATRIC_pathogen", type=str, required=True, help="The full path of PATRIC pathogen taxa list")
    parser.add_argument("--KEGG_pathogen", type=str, required=True, help="The full path of KEGG pathogen taxa list")
    parser.add_argument("--CARD_faa", type=str, required=True, help="The full path of CARD protein fasta file")
    parser.add_argument("--log_dir", type=str, required=False, default=str('/'.join(path_list[:(project_index+1)] + ['log_folder'])), help="the path of logfile folder")
    parser.add_argument("--log_name", type=str, required=False, default="add_kg_attributes.log", help="log file name")
    parser.add_argument("--outdir", type=str, help="The output dir")
    args = parser.parse_args()

    KEGG_api_link = 'http://rest.kegg.jp'

    ## set up logger
    logger = utils.get_logger(os.path.join(args.log_dir, args.log_name))

    ## setting for nodes
    logger.info(f"setting for nodes")
    node_tsv_file = open(args.kg_nodes, 'r')
    node_read_tsv = tsv.reader(node_tsv_file, delimiter="\t")
    node_header = pd.read_csv(args.kg_nodes_header, sep='\t', header=None)
    node_column_names = list(node_header.loc[0,:])
    node_KG_ids_index = node_column_names.index('node_id:ID')

    # ## setting for edges
    # logger.info(f"setting for edges")
    # edge_tsv_file = open(args.kg_edges, 'r')
    # edge_read_tsv = tsv.reader(edge_tsv_file, delimiter="\t")
    # edge_header = pd.read_csv(args.kg_edges_header, sep='\t', header=None)
    # edge_column_names = list(edge_header.loc[0,:])

    ## read ko hierarchy information
    logger.info(f"read ko hierarchy information")
    with open(args.ko_hierarchy, 'rb') as file_in:
        ko_hierarchy = pickle.load(file_in)

    ## read ko gene list
    logger.info(f"read ko gene list")
    with open(args.ko_gene_list,'rb') as file_in:
        gene_list_per_ko = pickle.load(file_in)

    ## read PATRIC pathogen list
    logger.info(f"read PATRIC pathogen list")
    PATRIC_pathogen = pd.read_csv(args.PATRIC_pathogen, sep='\t', header=None)
    pathogen_dict = {f'tax:{str(x)}':1 for x in PATRIC_pathogen[0]}

    ## read KEGG pathogen list
    logger.info(f"read KEGG pathogen list")
    KEGG_pathogen = pd.read_csv(args.KEGG_pathogen, sep='\t', header=None)
    pathogen_dict.update({f'tax:{str(x)}':1 for x in KEGG_pathogen[0]})

    ## download KEGG Addendum to KEGG Orthology link
    logger.info(f"Download KEGG Addendum to KEGG Orthology link")
    ko_ag_tb = {}
    link = f"{KEGG_api_link}/link/ag/ko"
    res = requests.get(link)
    if res.status_code == 200:
        logger.info(f"Successuflly link KEGG Addendum to KEGG Orthology from {link}")
        table = pd.DataFrame([x.split('\t') for x in res.text.split('\n') if x.split('\t')[0]])
        table.columns = ['kegg_ko_id','kegg_ag_id']
    else:
        logger.error(f"Fail to link reactions to koids from {link}")
    for row in table.to_numpy():
        if row[0].split(':')[1] in ko_ag_tb:
            ko_ag_tb[row[0].split(':')[1]].update([row[1]])
        else:
            ko_ag_tb[row[0].split(':')[1]] = set([row[1]])

    ## mapping CARD protein to KEGG Orthology
    CARD_KO_set = set([f"ag:{record.id.split('|')[1].split('.')[0]}" for record in SeqIO.parse(args.CARD_faa, "fasta")])
    CARD_KO_dict = {}
    for ko_id in ko_ag_tb:
        if len(ko_ag_tb[ko_id].intersection(CARD_KO_set)) > 0:
            CARD_KO_dict[ko_id] = [len(ko_ag_tb[ko_id].intersection(CARD_KO_set)),'ǂ'.join(list(ko_ag_tb[ko_id].intersection(CARD_KO_set)))]

    ## add ko hierarchy and ko gene list
    new_node_list = []
    for row in node_read_tsv:
        row = [x if x!='nan' else '' for x in row]
        node_id = row[node_KG_ids_index]
        if node_id.split(':')[0] == 'ko':
            kegg_ko_id = node_id.split(':')[1]
            temp_res = ko_hierarchy.get(kegg_ko_id)
            row += [None]
            if temp_res is not None:
                row += ['ǂ'.join(['|'.join(x.split('|')[:-1]) for x in temp_res])]
            else:
                row += [None]
            temp_res = gene_list_per_ko.get(kegg_ko_id)
            if temp_res is not None:
                row += ['ǂ'.join(temp_res)]
            else:
                row += [None]
            temp_res = ko_ag_tb.get(kegg_ko_id)
            if temp_res is not None:
                row += ['ǂ'.join(list(temp_res))]
            else:
                row += [None]
            temp_res = CARD_KO_dict.get(kegg_ko_id)
            if temp_res is not None:
                row += temp_res
            else:
                row += [None, None]
            new_node_list += [row]
        elif node_id.split(':')[0] == 'tax':
            temp_res = pathogen_dict.get(node_id)
            if temp_res is not None:
                row += [True, None, None, None, None, None]
            else:
                row += [False, None, None, None, None, None]
            new_node_list += [row]
        else:
            row += [None, None, None, None, None, None]
            new_node_list += [row]
    node_column_names += ['is_pathogen','ko_hierarchy', 'ko_gene_list', 'Addendum_list', 'CARD_match_num', 'CARD_match_ids']
    new_nodes = pd.DataFrame(new_node_list)
    new_nodes.columns = node_column_names

    ## output nodes to tsv format
    new_nodes.loc[new_nodes['node_id:ID'].isna(),'node_id:ID'] = None
    new_nodes.loc[new_nodes['KG2_ids'].isna(),'KG2_ids'] = None
    new_nodes.loc[new_nodes['preferred_KG2_ids'].isna(),'preferred_KG2_ids'] = None
    new_nodes.loc[new_nodes['preferred_KG2_names'].isna(),'preferred_KG2_names'] = None
    new_nodes.loc[new_nodes['preferred_KG2_node_types'].isna(),'preferred_KG2_node_types'] = None
    new_nodes.loc[new_nodes['node_type'].isna(),'node_type'] = None
    new_nodes.loc[new_nodes['iri'].isna(),'iri'] = None
    new_nodes.loc[new_nodes['description'].isna(),'description'] = None
    new_nodes.loc[new_nodes['is_pathogen'].isna(),'is_pathogen'] = None
    new_nodes.loc[new_nodes['ko_hierarchy'].isna(),'ko_hierarchy'] = None
    new_nodes.loc[new_nodes['ko_gene_list'].isna(),'ko_gene_list'] = None
    new_nodes.loc[new_nodes['Addendum_list'].isna(),'Addendum_list'] = None
    new_nodes.loc[new_nodes['CARD_match_num'].isna(),'CARD_match_num'] = 0
    new_nodes.loc[new_nodes['CARD_match_ids'].isna(),'CARD_match_ids'] = None
    new_nodes[':LABEL'] = new_nodes['node_type']
    new_nodes.columns = ['node_id:ID','KG2_ids','preferred_KG2_ids','preferred_KG2_names','preferred_KG2_node_types','node_type','iri','description','is_pathogen','ko_hierarchy:string[]', 'ko_gene_list:string[]', 'Addendum_list:string[]', 'CARD_match_num:float', 'CARD_match_ids:string[]', ':LABEL']
    node_col_tsvfile_out = open(os.path.join(args.outdir, 'nodes_header.tsv'), 'w+')
    node_col_tsvwrite = tsv.writer(node_col_tsvfile_out, delimiter="\t",  quoting=tsv.QUOTE_MINIMAL)
    node_col_tsvwrite.writerows([list(new_nodes.columns)])
    node_col_tsvfile_out.close()
    node_tsvfile_out = open(os.path.join(args.outdir, 'nodes.tsv'), 'w+')
    node_tsvwrite = tsv.writer(node_tsvfile_out, delimiter="\t",  quoting=tsv.QUOTE_MINIMAL)
    rows = [list(row) for row in new_nodes.to_numpy()]
    node_tsvwrite.writerows(rows)
    node_tsvfile_out.close()

    # ## output edges to tsv format
    # edge_col_tsvfile_out = open(os.path.join(args.outdir, 'edges_header.tsv'), 'w+')
    # edge_col_tsvwrite = tsv.writer(edge_col_tsvfile_out, delimiter="\t",  quoting=tsv.QUOTE_MINIMAL)
    # edge_col_tsvwrite.writerows([list(new_edges.columns)])
    # edge_col_tsvfile_out.close()
    # edge_tsvfile_out = open(os.path.join(args.outdir, 'edges.tsv'), 'w+')
    # edge_tsvwrite = tsv.writer(edge_tsvfile_out, delimiter="\t",  quoting=tsv.QUOTE_MINIMAL)
    # rows = [list(row) for row in new_edges.to_numpy()]
    # edge_tsvwrite.writerows(rows)
    # edge_tsvfile_out.close()