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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ko", type=str, required=True, help="the path of all koids")
    parser.add_argument("--organism_tb", type=str, required=True, help="the path of organism table")
    parser.add_argument("--organism_ko_dir", type=str, required=True, help="the path to a directory that contains organism ko")
    parser.add_argument("--virus_ko_file", type=str, required=True, help="the path to a file that contains ko to virus genes")
    parser.add_argument("--organisms", type=str, nargs='*', help="Multiple options from Fungi, Archaea, Bacteria", default=['Archaea','Bacteria', 'Fungi'])
    parser.add_argument("--outdir", type=str, help="The output dir")
    parser.add_argument("--log_dir", type=str, required=False, default=str('/'.join(path_list[:(project_index+1)] + ['log_folder'])), help="the path of logfile folder")
    parser.add_argument("--log_name", type=str, required=False, default="gather_genes_per_ko.log", help="log file name")
    args = parser.parse_args()

    ## set up logger
    logger = utils.get_logger(os.path.join(args.log_dir, args.log_name))

    logger.info(f"Organism Taxa only consider {args.organisms}")
    args.organisms = [x.lower() for x in args.organisms]

    ## read all ko ids
    logger.info(f"Read all ko ids")
    all_ko_ids = pd.read_csv(args.ko, sep='\t', header=0)
    ko_dict = {x.split(':')[1]:[] for x in all_ko_ids['kegg_ko_id']}

    ## read organism table and ko to gene
    logger.info(f"Read organism table and ko to gene")
    organism_table = pd.read_csv(args.organism_tb, sep='\t', header=0)
    organism_table = organism_table.loc[organism_table.apply(lambda row: len(set(args.organisms).intersection(set(row[3].lower().split(';')))) > 0, axis=1),:].reset_index(drop=True)
    organism_table = organism_table.loc[~organism_table.taxaid.isna(),:].reset_index(drop=True)
    organism_table['org_code'] = organism_table['org_code'].astype(str)
    organism_table['taxaid'] = organism_table['taxaid'].astype(int).astype(str)
    organism_table['taxaid'] = organism_table['taxaid'].replace('^','tax:',regex=True)
    orgcode_to_taxaid = dict(zip(list(organism_table['org_code']),list(organism_table['taxaid'])))

    ko_to_gene = glob(os.path.join(args.organism_ko_dir,'*'))
    for path in tqdm(ko_to_gene):
        orgcode = os.path.basename(path).split('_')[0]
        if orgcode in orgcode_to_taxaid:
            temp_df = pd.read_csv(path, sep='\t', header=0)
            for row in temp_df.to_numpy():
                gene_id, ko_id = row
                ko_id = ko_id.split(':')[1]
                if ko_id in ko_dict:
                    ko_dict[ko_id] += [gene_id]

    ## read ko to virus genes
    logger.info(f"Read ko to virus genes")
    link_ko_to_virus_gene = pd.read_csv(args.virus_ko_file, sep='\t', header=0)
    for row in tqdm(link_ko_to_virus_gene.to_numpy()):
        gene_id, ko_id = row
        ko_id = ko_id.split(':')[1]
        if ko_id in ko_dict:
            ko_dict[ko_id] += [gene_id]

    with open(os.path.join(args.outdir,f"gene_list_per_ko.pkl"),'wb') as outfile:
        pickle.dump(ko_dict,outfile)
    