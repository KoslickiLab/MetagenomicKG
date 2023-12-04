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
from bs4 import BeautifulSoup
import time
import requests
import re
import pickle
from multiprocessing import Pool

def extract_taxid(gn_id):
    link = f"https://www.genome.jp/entry/{gn_id}"
    while True:
        try:
            r = requests.get(link)
            break
        except:
            print(f"Get something wrong when calling {link}. Try to sleep a while and retry")
            time.sleep(300)
    if r.status_code == 200:
        soup = BeautifulSoup(r.text, 'html.parser')
        try:
            return soup.find_all('a',{'href': re.compile('Taxonomy')})[0].text
        except:
            print(f"Error: Fail to get the taxa id for {gn_id}",flush=True)
            return None
    else:
        print(f"Error: Fail to get the taxa id for {gn_id}",flush=True)
        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--organisms", type=str, nargs='*', help="Multiple options from Fungi, Archaea, Bacteria", default=['Archaea','Bacteria', 'Fungi'])
    parser.add_argument("--log_dir", type=str, required=False, default=str('/'.join(path_list[:(project_index+1)] + ['log_folder'])), help="the path of logfile folder")
    parser.add_argument("--log_name", type=str, required=False, default="build_knowledge_graph.log", help="log file name")
    parser.add_argument("--outdir", type=str, help="The output dir")
    args = parser.parse_args()

    KEGG_api_link = 'http://rest.kegg.jp'

    ## set up logger
    logger = utils.get_logger(os.path.join(args.log_dir, args.log_name))

    logger.info(f"Organism Taxa only consider {args.organisms}")
    args.organisms = [x.lower() for x in args.organisms]

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    ## collect kg nodes
    logger.info('collect kg nodes')
    nodes = []
    if not os.path.exists(os.path.join(args.outdir, 'gn_to_taxaid.pkl')):
        gn_to_taxaid = {}
        exist_gn_to_taxaid = False
    else:
        logger.info(f'Read {os.path.join(args.outdir, "gn_to_taxaid.pkl")}')
        with open(os.path.join(args.outdir, 'gn_to_taxaid.pkl'),'rb') as file_in:
            gn_to_taxaid = pickle.load(file_in)
        exist_gn_to_taxaid = True
    kegg_file_list = glob(os.path.join(args.outdir,'kegg_*.txt'))
    for path in kegg_file_list:
        logger.info(f'Read {path}')
        temp_table = pd.read_csv(path, sep='\t', header=0)
        if path == os.path.join(args.outdir,'kegg_gns.txt'):
            if not exist_gn_to_taxaid:
                temp_target_gn_ids = list(temp_table['kegg_gn_id'])
                batch =list(range(0,len(temp_target_gn_ids),20))
                batch.append(len(temp_target_gn_ids))
                get_data = True
                final_res = []
                for i in tqdm(range(1, len(batch))):
                    start = batch[i-1]
                    end = batch[i]
                    while get_data:
                        try:
                            with Pool(processes=20) as excutator:
                                res = excutator.map(extract_taxid, temp_target_gn_ids[start:end])
                                get_data = False
                        except:
                            time.sleep(300)
                    final_res += [x for x in res]
                    time.sleep(5)
                    get_data = True
                temp = {key:value for key, value in zip(temp_target_gn_ids,final_res)}
                gn_to_taxaid.update({key:f'tax:{str(value)}' for key, value in temp.items() if value})
            else:
                logger.warning(f'Skip {os.path.join(args.outdir,"kegg_gns.txt")} because "gn_to_taxaid" variable is imported from {os.path.join(args.outdir, "gn_to_taxaid.pkl")}')
        else:
            nodes += list(set(temp_table.iloc[:,0]))

    ## add organism taxa ids
    logger.info(f"Read {os.path.join(args.outdir,'organisms','organism_table.txt')}")
    organism_table = pd.read_csv(os.path.join(args.outdir,'organisms','organism_table.txt'), sep='\t', header=0)
    excluded_taxids = list(organism_table.loc[~organism_table.apply(lambda row: len(set(args.organisms).intersection(set(row[3].lower().split(';')))) > 0, axis=1),'taxaid'].astype(int).astype(str).replace('^','tax:',regex=True))
    organism_table = organism_table.loc[organism_table.apply(lambda row: len(set(args.organisms).intersection(set(row[3].lower().split(';')))) > 0, axis=1),:].reset_index(drop=True)
    organism_table = organism_table.loc[~organism_table.taxaid.isna(),:].reset_index(drop=True)
    organism_table['org_code'] = organism_table['org_code'].astype(str)
    organism_table['taxaid'] = organism_table['taxaid'].astype(int).astype(str)
    organism_table['taxaid'] = organism_table['taxaid'].replace('^','tax:',regex=True)
    orgcode_to_taxaid = dict(zip(list(organism_table['org_code']),list(organism_table['taxaid'])))
    nodes += list(set(organism_table.loc[~organism_table['taxaid'].isna(),'taxaid']))
    if not exist_gn_to_taxaid:
        gn_to_taxaid.update({f'gn:{org_code}':taxaid for org_code, taxaid in organism_table[['org_code','taxaid']].to_numpy() if taxaid})

    ## add viruses taxa ids
    logger.info(f"Read {os.path.join(args.outdir,'viruses','virus_table.txt')}")
    virus_table = pd.read_csv(os.path.join(args.outdir,'viruses','virus_table.txt'), sep='\t', header=0)
    virus_table.loc[virus_table['rs_ncbi_seq_ids'].isnull(),'rs_ncbi_seq_ids'] = 'None'
    virus_table['rs_ncbi_seq_ids'] = virus_table['rs_ncbi_seq_ids'].apply(lambda row: eval(row) if type(row) is str else row)
    virus_table.loc[virus_table['gb_ncbi_seq_ids'].isnull(),'gb_ncbi_seq_ids'] = 'None'
    virus_table['gb_ncbi_seq_ids'] = virus_table['gb_ncbi_seq_ids'].apply(lambda row: eval(row) if type(row) is str else row)
    virus_table['taxaid'] = virus_table['taxaid'].astype(int).astype(str)
    virus_table['taxaid'] = virus_table['taxaid'].replace('^','tax:',regex=True)
    nodes += list(set(virus_table.loc[~virus_table['taxaid'].isna(),'taxaid']))
    if not exist_gn_to_taxaid:
        gn_to_taxaid.update({taxaid.replace('tax:','gn:'):taxaid for taxaid in virus_table['taxaid'] if taxaid})
    if not exist_gn_to_taxaid:
        with open(os.path.join(args.outdir, 'gn_to_taxaid.pkl'),'wb') as file_out:
            pickle.dump(gn_to_taxaid, file_out)

    ## exclude the organisms that don't belong to args.organisms
    gn_to_taxaid = {key:value for key, value in gn_to_taxaid.items() if value not in excluded_taxids}
    nodes += list(set(gn_to_taxaid.values()))
    nodes = list(set(nodes))

    ## integrate all link files
    logger.info('integrate all link files')
    temp_link_file_list = glob(os.path.join(args.outdir,'link_*'))
    temp_kg = pd.DataFrame(columns=['source','target'])
    for file_path in temp_link_file_list:
        logger.info(f'Read {file_path}')
        infile = pd.read_csv(file_path, sep='\t', header=0)
        infile.columns = ['source','target']
        if 'path' in os.path.basename(file_path):
            infile['source'] = infile['source'].replace('path:\D*','path:map',regex=True)
            infile['target'] = infile['target'].replace('path:\D*','path:map',regex=True)
        # if 'enzyme' in os.path.basename(file_path):
        #     infile['source'] = infile['source'].replace('ec:1.5.1.541','ec:1.5.1.-', regex=True).replace('ec:5.6.2.$','ec:5.6.2.-', regex=True)
        #     infile['target'] = infile['target'].replace('ec:1.5.1.541','ec:1.5.1.-', regex=True).replace('ec:5.6.2.$','ec:5.6.2.-', regex=True)
        if ('gn_to' in os.path.basename(file_path)) or ('to_gn' in os.path.basename(file_path)):
            infile['source'] = infile['source'].apply(lambda row: gn_to_taxaid.get(row) if 'gn:' in row else row)
            infile['target'] = infile['target'].apply(lambda row: gn_to_taxaid.get(row) if 'gn:' in row else row)
            infile = infile.dropna().reset_index(drop=True)
        temp_kg = pd.concat([temp_kg,infile])
    ## remove duplicate edges
    temp_kg = temp_kg.drop_duplicates().reset_index(drop=True)
    ## make edge bidirectional
    temp_kg_rev = temp_kg[['target','source']]
    temp_kg_rev.columns = ['source','target']
    kg = pd.concat([temp_kg,temp_kg_rev]).reset_index(drop=True)

    ## connect kg to virus information
    logger.info('connect kg to virus information')
    temp_link_file_list = glob(os.path.join(args.outdir,'viruses','link_*'))
    temp_kg = pd.DataFrame(columns=['source','target'])
    for file_path in temp_link_file_list:
        infile = pd.read_csv(file_path, sep='\t', header=0)
        infile['kegg_gene_id'] = infile['kegg_gene_id'].replace(':','_',regex=True).replace('^','gene:',regex=True)
        infile.columns = ['source','target']
        temp_kg = pd.concat([temp_kg,infile])

    temp_link_file_list = glob(os.path.join(args.outdir,'viruses','kegg_gene_info','*'))
    for file_path in temp_link_file_list:
        infile = pd.read_csv(file_path, sep='\t', header=0)
        infile = infile[['taxaid','kegg_gene_id']]
        infile['kegg_gene_id'] = infile['kegg_gene_id'].replace(':','_',regex=True).replace('^','gene:',regex=True)
        infile.taxaid = infile.taxaid.astype(str)
        infile['taxaid'] = infile['taxaid'].replace('^','tax:',regex=True)
        infile['kegg_gene_id'] = infile['kegg_gene_id'].replace(':','_',regex=True).replace('^','gene:',regex=True)
        infile.columns = ['source','target']
        temp_kg = pd.concat([temp_kg,infile])
    ## remove duplicate edges
    temp_kg = temp_kg.drop_duplicates().reset_index(drop=True)
    ## make edge bidirectional
    temp_kg_rev = temp_kg[['target','source']]
    temp_kg_rev.columns = ['source','target']
    kg = pd.concat([kg,temp_kg,temp_kg_rev]).reset_index(drop=True)
    kg = kg.drop_duplicates().reset_index(drop=True)

    ## connect kg to organism information
    logger.info('connect kg to organism information')
    temp_link_dir_list = glob(os.path.join(args.outdir,'organisms','*_to_*'))
    temp_kg = pd.DataFrame(columns=['source','target'])
    for file_path1 in temp_link_dir_list:
        logger.info(f'processing file {file_path1}')
        if os.path.basename(file_path1) == 'convert_geneid_to_uniprot':
            continue
        if os.path.basename(file_path1) == 'pathway_to_gene':
            temp_path_list = glob(os.path.join(file_path1,'*'))
            for file_path2 in tqdm(temp_path_list):
                orgcode = os.path.basename(file_path2).split('_')[0]
                taxaid = orgcode_to_taxaid.get(orgcode)
                if taxaid is None:
                    continue
                infile = pd.read_csv(file_path2, sep='\t', header=0)
                infile['kegg_gene_id'] = infile['kegg_gene_id'].replace(':','_',regex=True).replace('^','gene:',regex=True)
                infile['kegg_pathway_id'] = infile['kegg_pathway_id'].replace('path:\D*','path:map',regex=True)
                infile.columns = ['source','target']
                temp_kg = pd.concat([temp_kg,infile])
        elif os.path.basename(file_path1) == 'module_to_gene':
            temp_path_list = glob(os.path.join(file_path1,'*'))
            for file_path2 in tqdm(temp_path_list):
                orgcode = os.path.basename(file_path2).split('_')[0]
                taxaid = orgcode_to_taxaid.get(orgcode)
                if taxaid is None:
                    continue
                infile = pd.read_csv(file_path2, sep='\t', header=0)
                infile['kegg_gene_id'] = infile['kegg_gene_id'].replace(':','_',regex=True).replace('^','gene:',regex=True)
                infile['kegg_module_id'] = infile['kegg_module_id'].replace(':\D*_',':',regex=True)
                infile.columns = ['source','target']
                temp_kg = pd.concat([temp_kg,infile])
        else:
            temp_path_list = glob(os.path.join(file_path1,'*'))
            for file_path2 in tqdm(temp_path_list):
                orgcode = os.path.basename(file_path2).split('_')[0]
                taxaid = orgcode_to_taxaid.get(orgcode)
                if taxaid is None:
                    continue
                infile = pd.read_csv(file_path2, sep='\t', header=0)
                infile['kegg_gene_id'] = infile['kegg_gene_id'].replace(':','_',regex=True).replace('^','gene:',regex=True)
                # col_name = [col_name for col_name in list(infile.columns) if col_name != 'kegg_gene_id'][0]
                # infile[col_name] = infile[col_name].replace('ec:1.5.1.541','ec:1.5.1.-', regex=True).replace('ec:5.6.2.$','ec:5.6.2.-', regex=True)
                infile.columns = ['source','target']
                temp_kg = pd.concat([temp_kg,infile])
    ## remove duplicate edges
    temp_kg = temp_kg.drop_duplicates().reset_index(drop=True)
    ## make edge bidirectional
    temp_kg_rev = temp_kg[['target','source']]
    temp_kg_rev.columns = ['source','target']
    kg = pd.concat([kg,temp_kg,temp_kg_rev]).reset_index(drop=True)
    kg = kg.drop_duplicates().reset_index(drop=True)

    ## connect organism genes to organism
    logger.info('connect organism genes to organism')
    temp_link_file_list = glob(os.path.join(args.outdir,'organisms','kegg_gene_info','*'))
    temp_kg = pd.DataFrame(columns=['source','target'])
    for path in tqdm(temp_link_file_list):
        orgcode = os.path.basename(path).split('_')[0]
        taxaid = orgcode_to_taxaid.get(orgcode)
        if taxaid is not None:
            infile = pd.read_csv(path, sep='\t', header=0)
            infile['kegg_gene_id'] = infile['kegg_gene_id'].replace(':','_',regex=True).replace('^','gene:',regex=True)
            temp = pd.DataFrame(list(zip([taxaid]*len(infile['kegg_gene_id']),list(infile['kegg_gene_id']))))
            temp.columns = ['source','target']
            temp_kg = pd.concat([temp_kg,temp])
    ## remove duplicate edges
    temp_kg = temp_kg.drop_duplicates().reset_index(drop=True)
    ## make edge bidirectional
    temp_kg_rev = temp_kg[['target','source']]
    temp_kg_rev.columns = ['source','target']
    kg = pd.concat([kg,temp_kg,temp_kg_rev]).reset_index(drop=True)
    kg = kg.drop_duplicates().reset_index(drop=True)

    ## connect organism pathways to organism
    logger.info('connect organism pathways to organism')
    temp_link_file_list = glob(os.path.join(args.outdir,'organisms','kegg_pathway_info','*'))
    temp_kg = pd.DataFrame(columns=['source','target'])
    for path in tqdm(temp_link_file_list):
        orgcode = os.path.basename(path).split('_')[0]
        taxaid = orgcode_to_taxaid.get(orgcode)
        if taxaid is not None:
            infile = pd.read_csv(path, sep='\t', header=0)
            infile['kegg_pathway_id'] = infile['kegg_pathway_id'].replace('path:\D*','path:map',regex=True)
            temp = pd.DataFrame(list(zip([taxaid]*len(infile['kegg_pathway_id']),list(infile['kegg_pathway_id']))))
            temp.columns = ['source','target']
            temp_kg = pd.concat([temp_kg,temp])
    ## remove duplicate edges
    temp_kg = temp_kg.drop_duplicates().reset_index(drop=True)
    ## make edge bidirectional
    temp_kg_rev = temp_kg[['target','source']]
    temp_kg_rev.columns = ['source','target']
    kg = pd.concat([kg,temp_kg,temp_kg_rev]).reset_index(drop=True)
    kg = kg.drop_duplicates().reset_index(drop=True)

    ## check if there are some nodes in the edges but miss in node list
    edge_nodes = set(list(kg['source']) + list(kg['target']))
    edge_nodes = set([node for node in edge_nodes if node.split(':')[0] != 'gene'])
    if len(edge_nodes.difference(set(nodes))) > 0:
        additional_nodes = list(edge_nodes.difference(set(nodes)))
        logger.warning(f"Found the following nodes {additional_nodes} are not in the KEGG nodes collected from 'list' API. They will be dropped out.")
        kg = kg.loc[~((kg['source'].isin(additional_nodes)) | (kg['target'].isin(additional_nodes))),:].reset_index(drop=True)

    ## exclude edges that contain organisms that doesn't belong to args.organisms
    kg = kg.loc[~((kg['source'].isin(excluded_taxids)) | (kg['target'].isin(excluded_taxids))),:].reset_index(drop=True)

    ## save kg
    if not os.path.exists(os.path.join(args.outdir,'KG_info')):
        os.makedirs(os.path.join(args.outdir,'KG_info'))
    kg.to_csv(os.path.join(args.outdir,'KG_info','kg.txt'), sep='\t', index=None)

    ## save kg nodes
    pd.DataFrame(sorted(nodes), columns=['node_name']).to_csv(os.path.join(args.outdir,'KG_info','nodes.txt'), sep='\t', index=None)
