import os, sys
project_name = os.environ['project_name']
path_list = os.path.abspath(__file__).split('/')
project_index = path_list.index(project_name)
path_to_scripts = '/'.join(path_list[:(project_index+1)] + ['scripts'])
sys.path.append(path_to_scripts)
from utility_scripts import utils
import pandas as pd
from tqdm import trange
from multiprocessing import Pool
from itertools import combinations
import pickle
import argparse


def process_fun(gene_node):
    neighbor_nodes = gene_connect[gene_node]
    reconnect = pd.DataFrame([x for x in combinations(neighbor_nodes, 2)])
    if len(reconnect) != 0:
        reconnect.columns = ['source','target']
        ## make edge bidirectional
        reconnect_rev = reconnect[['target','source']]
        reconnect_rev.columns = ['source','target']
        kg_temp = pd.concat([reconnect,reconnect_rev]).reset_index(drop=True)
    else:
        kg_temp = pd.DataFrame(columns=['source','target'])

    return kg_temp

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--kg", type=str, help="The full path of original kg")
    parser.add_argument("--batchsize", type=int, required=False, default=50000, help="batch size for parallel running")
    parser.add_argument("--log_dir", type=str, required=False, default=str('/'.join(path_list[:(project_index+1)] + ['log_folder'])), help="the path of logfile folder")
    parser.add_argument("--log_name", type=str, required=False, default="filter_gene_nodes.log", help="log file name")
    args = parser.parse_args()

    ## set up logger
    logger = utils.get_logger(os.path.join(args.log_dir, args.log_name))

    if os.path.exists(args.kg):
        kg = pd.read_csv(args.kg, sep='\t', header=0)
    else:
        logger.error(f"Not Found {args.kg}")
        exit()
    dest_path = os.path.dirname(args.kg)

    if not os.path.exists(os.path.join(dest_path,'gene_connect.pkl')):
        logger.info(f"Find all nodes connected to 'gene' node and store them into a dictionary and save it as 'gene_connect.pkl'")
        gene_connect = dict()
        for index in trange(len(kg)):
            source = kg.loc[index,'source']
            target = kg.loc[index,'target']
            if source.split(':')[0] == 'gene':
                if source in gene_connect:
                    gene_connect[source].update(set([target]))
                else:
                    gene_connect[source] = set([target])
            elif target.split(':')[0] == 'gene':
                if target in gene_connect:
                    gene_connect[target].update(set([source]))
                else:
                    gene_connect[target] = set([source])
        with open(os.path.join(dest_path,'gene_connect.pkl'),'wb') as outfile:
            pickle.dump(gene_connect,outfile)
    else:
        logger.info(f"File 'gene_connect.pkl' is detected!")
        with open(os.path.join(dest_path,'gene_connect.pkl'),'rb') as infile:
            gene_connect = pickle.load(infile)

    ## remove all edges with 'gene' nodes
    kg = kg.loc[~((kg['source'].str.contains('gene:')) | (kg['target'].str.contains('gene:'))),:].reset_index(drop=True)

    ## get gene nodes
    gene_nodes = list(gene_connect.keys())

    if not os.path.exists(os.path.join(dest_path,'reconnect.pkl')):
        logger.info(f"Reconnect neighbor nodes that are originally connected to the gene nodes")
        reconnect = []
        batch =list(range(0,len(gene_nodes),args.batchsize))
        batch.append(len(gene_nodes))
        for i in trange(1, len(batch)):
            start = batch[i-1]
            end = batch[i]
            with Pool(processes=64) as excutator:
                res = excutator.map(process_fun, gene_nodes[start:end])
                reconnect += [pd.concat(res).reset_index(drop=True)]
        reconnect = pd.concat(reconnect).reset_index(drop=True)
        ## remove duplicate edges
        reconnect = reconnect.drop_duplicates().reset_index(drop=True)
        with open(os.path.join(dest_path,'reconnect.pkl'),'wb') as outfile:
            pickle.dump(reconnect,outfile)
    else:
        logger.info(f"File 'reconnect.pkl' is detected!")
        with open(os.path.join(dest_path,'reconnect.pkl'),'rb') as infile:
            reconnect = pickle.load(infile)

    ## add reconnected neighbor nodes into knowledge graph
    kg = pd.concat([kg,reconnect]).reset_index(drop=True)
    kg = kg.drop_duplicates().reset_index(drop=True)

    kg.to_csv(os.path.join(dest_path,'kg_gene_filtered.txt'), sep='\t', index=None)

