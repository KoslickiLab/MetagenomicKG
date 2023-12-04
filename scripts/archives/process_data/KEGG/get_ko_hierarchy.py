import os, sys
project_name = os.environ['project_name']
path_list = os.path.abspath(__file__).split('/')
project_index = path_list.index(project_name)
path_to_scripts = '/'.join(path_list[:(project_index+1)] + ['scripts'])
sys.path.append(path_to_scripts)
from utility_scripts import utils
import argparse
import logging
from tqdm import tqdm, trange
import pickle
from itertools import chain
import pandas as pd
import requests
from bs4 import BeautifulSoup
import time
import re
from glob import glob
import json

def organize_hierarchy(hiearchy_json, regex='^\w?\w?\d{5} ', prefix='', stop_level=None):

    def _iterate_multidimensional(res, hiearchy_json, res_list):
        if isinstance(hiearchy_json,dict):
            for k,v in hiearchy_json.items():
                if k == 'name':
                    temp_name = re.sub(' \[.*\]','',hiearchy_json['name'])
                    res += f"|{temp_name}"
                    if re.search(regex, hiearchy_json['name']) is not None:
                        res_list += [res]
                elif k == 'children':
                    for elem in hiearchy_json['children']:
                        _iterate_multidimensional(res, elem, res_list)
        else:
            self.logger.error(f"{hiearchy_json} is not dictionary")
            raise

    res_list = []
    _iterate_multidimensional('', hiearchy_json, res_list)
    if stop_level is None:
        res_dict = {prefix+string.split('|')[-1].split(' ')[0]:'|'.join(string.split('|')[1:]) for string in res_list}
        return res_dict
    else:
        res_dict = {prefix+string.split('|')[-1].split(' ')[0].split('\t')[0]:'|'.join(string.split('|')[1:(stop_level+1)]) for string in res_list}
        return res_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, help="The output dir")
    parser.add_argument("--log_dir", type=str, required=False, default=str('/'.join(path_list[:(project_index+1)] + ['log_folder'])), help="the path of logfile folder")
    parser.add_argument("--log_name", type=str, required=False, default="get_ko_hierarchy.log", help="log file name")
    args = parser.parse_args()

    KEGG_api_link = 'http://rest.kegg.jp'

    ## set up logger
    logger = utils.get_logger(os.path.join(args.log_dir, args.log_name))

    ## get brite table
    link = f"{KEGG_api_link}/list/brite"
    res = requests.get(link)
    if res.status_code == 200:
        brite_table = pd.DataFrame([x.split('\t') for x in res.text.split('\n') if x.split('\t')[0]])
        brite_table.columns = ['kegg_brite_id','desc']
    else:
        logger.error(f"Fail to download KEGG brite information from {link}")
        exit()

    ## download KEGG KO associated hierarchy and process hierarchy
    logger.info(f"Download KEGG KO associated hierarchy")
    ko_hierarchy_dict = dict()
    for brite_id, desc in brite_table.to_numpy():
        logger.info(f"Processing brite id {brite_id}")
        check_link = f"{KEGG_api_link}/get/{brite_id}"
        res = requests.get(check_link)
        if res.status_code == 200:
            m = re.search('K\d{5}', res.text)
        else:
            logger.error(f"Fail to download KEGG brite information from {check_link}")
            continue
        if m is None:
            logger.warning(f"Brite ID {brite_id} doesn't contain KO ids and thus skip it.")
            continue
        link = f"{KEGG_api_link}/get/{brite_id}/json"
        json_res = requests.get(link)
        if json_res.status_code == 200:
            temp_dict = organize_hierarchy(json_res.json(), regex='^K\d{5} ')
            for key in temp_dict:
                if key in ko_hierarchy_dict:
                    ko_hierarchy_dict[key] += [temp_dict[key]]
                else:
                    ko_hierarchy_dict[key] = [temp_dict[key]]
        else:
            logger.error(f"Fail to download KEGG brite information from {link}")

    with open(os.path.join(args.outdir,f"kegg_ko_brite_info.pkl"),'wb') as outfile:
        pickle.dump(ko_hierarchy_dict,outfile)
