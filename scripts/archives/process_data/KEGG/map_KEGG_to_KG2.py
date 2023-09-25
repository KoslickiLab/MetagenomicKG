import os, sys
import argparse
project_name = os.environ['project_name']
path_list = os.path.abspath(__file__).split('/')
project_index = path_list.index(project_name)
path_to_scripts = '/'.join(path_list[:(project_index+1)] + ['scripts'])
sys.path.append(path_to_scripts)
from utility_scripts import utils, node_synonymizer
import pandas as pd
from glob import glob
from tqdm import tqdm
from bs4 import BeautifulSoup
import time
import requests
import re
import pickle
from multiprocessing import Pool

def IdMapping(KEGG_id):

    tag = KEGG_id.split(':')[0]

    if tag == 'cpd':
        KG2_id = KEGG_id.replace('cpd:','KEGG.COMPOUND:')
    elif tag == 'dr':
        KG2_id = KEGG_id.replace('dr:','KEGG.DRUG:')
    elif tag == 'ec':
        KG2_id = KEGG_id.replace('ec:','KEGG.ENZYME:')
    elif tag == 'gl':
        KG2_id = KEGG_id.replace('gl:','KEGG.GLYCAN:')
    elif tag == 'rn':
        KG2_id = KEGG_id.replace('rn:','KEGG.REACTION:')
    elif tag == 'tax':
        KG2_id = KEGG_id.replace('tax:','NCBITaxon:')
    elif tag == 'ds':
        link = f"https://www.genome.jp/entry/{KEGG_id}"
        while True:
            try:
                r = requests.get(link)
                break
            except:
                print(f"Get something wrong when calling {link}. Try to sleep a while and retry")
                time.sleep(300)
        if r.status_code == 200:
            soup = BeautifulSoup(r.text, 'html.parser')
            KG2_id = []
            ## extract MESH id
            try:
                temp_res = soup.find_all('a',{'href': re.compile('mesh')})
                KG2_id += [f"MESH:{x.text}" for x in temp_res]
            except:
                print(f"Error: Fail to get the MESH id for {KEGG_id}",flush=True)
            ## extract OMIN ids
            try:
                temp_res = soup.find_all('a',{'href': re.compile('omim.org')})
                KG2_id += [f"OMIM:{x.text}" for x in temp_res]
            except:
                print(f"Error: Fail to get the OMIN ids for {KEGG_id}",flush=True)
            if len(KG2_id) == 0:
                KG2_id = None
            else:
                KG2_id = '##'.join(KG2_id)
        else:
            print(f"Error: Fail to get the taxa id for {KEGG_id}",flush=True)
            KG2_id = None
    else:
        KG2_id = None
        
    return [KEGG_id, KG2_id]


def get_preferred_info(KG2_ids):

    if KG2_ids is None:
        return [None, None, None]
    else:
        KG2_id_list = KG2_ids.split('##')
        normalizer = nodesynonymizer.get_canonical_curies(KG2_id_list)
        preferred_KG2_ids = []
        preferred_KG2_names = []
        preferred_KG2_nodetype = []
        for node in normalizer:
            if normalizer[node] is not None:
                preferred_KG2_ids += [normalizer[node]['preferred_curie']]
                preferred_KG2_names += [normalizer[node]['preferred_name']]
                preferred_KG2_nodetype += [normalizer[node]['preferred_category']]
        if len(preferred_KG2_ids) > 1:
            preferred_KG2_ids = '##'.join(list(set(preferred_KG2_ids)))
            preferred_KG2_names = '##'.join(list(set(preferred_KG2_names)))
            preferred_KG2_nodetype = '##'.join(list(set(preferred_KG2_nodetype)))
        else:
            preferred_KG2_ids = '##'.join(preferred_KG2_ids)
            preferred_KG2_names = '##'.join(preferred_KG2_names)
            preferred_KG2_nodetype = '##'.join(preferred_KG2_nodetype)
    
    return [preferred_KG2_ids, preferred_KG2_names, preferred_KG2_nodetype]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--KG_nodes", type=str, required=True, help="Multiple options from Fungi, Archaea, Bacteria", default=str('/'.join(path_list[:(project_index+1)] + ['data', 'KEGG', 'KG_info', 'nodes.txt'])))
    parser.add_argument("--log_dir", type=str, required=False, default=str('/'.join(path_list[:(project_index+1)] + ['log_folder'])), help="the path of logfile folder")
    parser.add_argument("--log_name", type=str, required=False, default="map_KEGG_to_KG2.log", help="log file name")
    parser.add_argument("--processes", type=int, required=False, default=64, help="number of cpu to use")
    parser.add_argument("--outdir", type=str, help="The output dir", default=str('/'.join(path_list[:(project_index+1)] + ['data', 'KEGG', 'KG_info'])))
    args = parser.parse_args()

    KEGG_api_link = 'http://rest.kegg.jp'

    ## set up logger
    logger = utils.get_logger(os.path.join(args.log_dir, args.log_name))

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    ## read node information
    logger.info(f'read node information from {args.KG_nodes}')
    if not os.path.exists(args.KG_nodes):
        logger.error(f'Not found {args.KG_nodes}')
        exit()
    else:
        KG_nodes = pd.read_csv(args.KG_nodes, sep='\t', header=0)
        KG_nodes.columns = ['node_KG_ids']

    ## map KEGG ids to KG2 ids in parallel
    ds_id_list = [node for node in KG_nodes['node_KG_ids'] if node.split(':')[0]=='ds']
    batch =list(range(0,len(ds_id_list),20))
    batch.append(len(ds_id_list))
    get_data = True
    final_res = []
    for i in tqdm(range(1, len(batch))):
        start = batch[i-1]
        end = batch[i]
        while get_data:
            try:
                with Pool(processes=20) as p:
                    res = p.map(IdMapping, ds_id_list[start:end])
                    get_data = False
            except:
                time.sleep(300)
        final_res += [x for x in res]
        time.sleep(5)
        get_data = True
    
    nonds_id_list = [node for node in KG_nodes['node_KG_ids'] if node.split(':')[0]!='ds']
    with Pool(processes=args.processes) as p:
        res = p.map(IdMapping, nonds_id_list)
        final_res += [x for x in res]
    
    KG_nodes = pd.DataFrame(final_res)
    KG_nodes.columns = ['node_KG_ids', 'KG2_ids']

    ## find preferred KG2 ids
    nodesynonymizer = node_synonymizer.NodeSynonymizer()
    KG_nodes.apply(lambda row: get_preferred_info(row[1]), axis=1, result_type='expand')
    temp_df = KG_nodes.apply(lambda row: get_preferred_info(row[1]), axis=1, result_type='expand')
    temp_df.columns = ['preferred_KG2_ids', 'preferred_KG2_names', 'preferred_KG2_node_types']
    KG_nodes = pd.concat([KG_nodes, temp_df], axis=1)
    KG_nodes.to_csv(os.path.join(args.outdir, os.path.basename(args.KG_nodes)), sep='\t', index=None)