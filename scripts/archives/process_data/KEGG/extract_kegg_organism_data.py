
import os, sys
import argparse
project_name = os.environ['project_name']
path_list = os.path.abspath(__file__).split('/')
project_index = path_list.index(project_name)
path_to_scripts = '/'.join(path_list[:(project_index+1)] + ['scripts'])
sys.path.append(path_to_scripts)
from utility_scripts import utils
from tqdm import tqdm, trange
import pickle
from itertools import chain
import pandas as pd
import pytaxonkit
import requests
from bs4 import BeautifulSoup
from multiprocessing import Pool, cpu_count
import time
import re
from glob import glob

# def get_logger():
#     logger = logging.getLogger()
#     logger.setLevel(logging.DEBUG)
#     formatter = logging.Formatter('%(asctime)s  [%(levelname)s]  %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
#     ch = logging.StreamHandler(sys.stdout)
#     ch.setFormatter(formatter)
#     logger.addHandler(ch)
#     return logger

def extract_taxid(org_code):
    link = f"https://www.genome.jp/kegg-bin/show_organism?org={org_code}"
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
            print(f"Error: Fail to get the taxa id for {org_code}",flush=True)
            return None
    else:
        print(f"Error: Fail to get the taxa id for {org_code}",flush=True)
        return None

def extract_ref_seq_id(params):
    org_code, seq_db = params
    link = f"https://www.genome.jp/kegg-bin/show_organism?org={org_code}"
    while True:
        try:
            r = requests.get(link)
            break
        except:
            print(f"Get something wrong when calling {link}. Try to sleep a while and retry")
            time.sleep(300)
    if r.status_code == 200:
        soup = BeautifulSoup(r.text, 'html.parser')
        if seq_db == 'gb':
            try:
                return [x.text for x in soup.find_all('a',{'href': re.compile('www.ncbi.nlm.nih.gov\/nuccore')})]
            except:
                print(f"Error: Fail to get the ref seq id for {org_code}",flush=True)
                return None
        elif seq_db == 'rs':
            try:
                return [x.text for x in soup.find_all('a',{'href': re.compile('www.genome.jp/dbget-bin/www_bget\?refseq')})]
            except:
                print(f"Error: Fail to get the ref seq id for {org_code}",flush=True)
                return None
    else:
        print(f"Error: Fail to get the ref seq id for {org_code}",flush=True)
        return None

def download_kegg_gene(params):
    org_code, out_loc = params
    link1 = f"{KEGG_api_link}/list/{org_code}"
    while True:
        try:
            r = requests.get(link1)
            break
        except:
            print(f"Get something wrong when calling {link1}. Try to sleep a while and retry")
            time.sleep(300)
    if r.status_code == 200:
        print(f"Successuflly download gene/protein information from {link1}")
        table1 = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])
        table1.columns = ['kegg_gene_id','desc']
    else:
        print(f"Error: Fail to download gene/protein information from {link1}")
        return 0
    # link2 = f"http://rest.kegg.jp/conv/{org_code}/uniprot"
    # r = requests.get(link2)
    # if r.status_code == 200:
    #     table2 = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])
    #     table2.columns = ['uniprot_id','kegg_gene_id']
    # else:
    #     print(f"Error: Fail to convert kegg gene ids to uniprot ids via {link2}") 
    #     return 0 
    # combine = table1.merge(table2, on='kegg_gene_id').reset_index(drop=True)
    # combine.uniprot_id = combine.uniprot_id.str.replace('up','UniProt')
    table1.to_csv(os.path.join(out_loc,f"{org_code}_kegg_genes.txt"), sep='\t', index=None)
    return 1

def download_kegg_pathway(params):
    org_code, out_loc = params
    link1 = f"{KEGG_api_link}/list/pathway/{org_code}"
    while True:
        try:
            r = requests.get(link1)
            break
        except:
            print(f"Get something wrong when calling {link1}. Try to sleep a while and retry")
            time.sleep(300)
    if r.status_code == 200:
        print(f"Successuflly download pathway information from {link1}")
        table1 = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])
        table1.columns = ['kegg_pathway_id','desc']
    else:
        print(f"Error: Fail to download pathway information from {link1}")
        return 0
    table1.to_csv(os.path.join(out_loc,f"{org_code}_kegg_pathways.txt"), sep='\t', index=None)
    return 1

def link_pathway_to_gene(params):
    org_code, out_loc = params
    link1 = f"{KEGG_api_link}/link/pathway/{org_code}"
    while True:
        try:
            r = requests.get(link1)
            break
        except:
            print(f"Get something wrong when calling {link1}. Try to sleep a while and retry")
            time.sleep(300)
    if r.status_code == 200:
        print(f"Successuflly link pathways to genes from {link1}")
        table1 = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])
        table1.columns = ['kegg_gene_id','kegg_pathway_id']
    else:
        print(f"Error: Fail to link pathways to genes from {link1}")
        return 0
    table1.to_csv(os.path.join(out_loc,f"{org_code}_link_pathway_to_gene.txt"), sep='\t', index=None)
    return 1

def link_ko_to_gene(params):
    org_code, out_loc = params
    link1 = f"{KEGG_api_link}/link/ko/{org_code}"
    while True:
        try:
            r = requests.get(link1)
            break
        except:
            print(f"Get something wrong when calling {link1}. Try to sleep a while and retry")
            time.sleep(300)
    if r.status_code == 200:
        print(f"Successuflly link koids to genes from {link1}")
        table1 = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])
        table1.columns = ['kegg_gene_id','kegg_ko_id']
    else:
        print(f"Error: Fail to link koids to genes from {link1}")
        return 0
    table1.to_csv(os.path.join(out_loc,f"{org_code}_link_ko_to_gene.txt"), sep='\t', index=None)
    return 1

def link_module_to_gene(params):
    org_code, out_loc = params
    link1 = f"{KEGG_api_link}/link/module/{org_code}"
    while True:
        try:
            r = requests.get(link1)
            break
        except:
            print(f"Get something wrong when calling {link1}. Try to sleep a while and retry")
            time.sleep(300)
    if r.status_code == 200:
        print(f"Successuflly link modules to genes from {link1}")
        table1 = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])
        table1.columns = ['kegg_gene_id','kegg_module_id']
    else:
        print(f"Error: Fail to link modules to genes from {link1}")
        return 0
    table1.to_csv(os.path.join(out_loc,f"{org_code}_link_module_to_gene.txt"), sep='\t', index=None)
    return 1

def link_enzyme_to_gene(params):
    org_code, out_loc = params
    link1 = f"{KEGG_api_link}/link/enzyme/{org_code}"
    while True:
        try:
            r = requests.get(link1)
            break
        except:
            print(f"Get something wrong when calling {link1}. Try to sleep a while and retry")
            time.sleep(300)
    if r.status_code == 200:
        print(f"Successuflly link enzymes to genes from {link1}")
        table1 = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])
        table1.columns = ['kegg_gene_id','kegg_enzyme_id']
    else:
        print(f"Error: Fail to link enzymes to genes from {link1}")
        return 0
    table1.to_csv(os.path.join(out_loc,f"{org_code}_link_enzyme_to_gene.txt"), sep='\t', index=None)
    return 1

def convert_geneid_to_uniprot(params):
    org_code, out_loc = params
    link1 = f"{KEGG_api_link}/conv/uniprot/{org_code}"
    while True:
        try:
            r = requests.get(link1)
            break
        except:
            print(f"Get something wrong when calling {link1}. Try to sleep a while and retry")
            time.sleep(300)
    if r.status_code == 200:
        print(f"Successuflly convert gene ids to uniprot ids from {link1}")
        table1 = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])
        table1.columns = ['kegg_gene_id','uniprot_gene_id']
    else:
        print(f"Error: Fail to convert gene ids to uniprot ids from {link1}")
        return 0
    table1.to_csv(os.path.join(out_loc,f"{org_code}_convert_gene_kegg_to_uniprot.txt"), sep='\t', index=None)
    return 1

def extract_taxaid_seq(inlist):
    temp = '|'.join(inlist)
    if 'ORTHOLOGY' in temp:
        koid = 'ko:'+[re.sub('\s.*','',re.sub('ORTHOLOGY\s*','',line)) for line in inlist if 'ORTHOLOGY' in line][0]
    else:
        koid = None
    if 'AASEQ' in temp:
        aaseq = re.sub('\d*','','|'.join(inlist).split('AASEQ       ')[1].split('|COMMENT     ')[0].split('|NTSEQ     ')[0]).replace('|            ','').replace('|','')
    else:
        aaseq = None
    if 'NTSEQ' in temp:
        ntseq = re.sub('\d*','','|'.join(inlist).split('NTSEQ       ')[1].split('|COMMENT     ')[0]).replace('|            ','').replace('|','')
    else:
        ntseq = None
    return koid, aaseq, ntseq

def process_query(instr):
    link = KEGG_api_link + f"/get/{instr}"
    while True:
        try:
            res = requests.get(link)
            break
        except:
            print(f"Get something wrong when calling {link}. Try to sleep a while and retry")
            time.sleep(300)
    if res.status_code == 200:
        print(f"Successuflly extract info from {link}")
        res = [tuple([a])+b for a, b in zip(instr.split('+'),list(map(extract_taxaid_seq,[x.split('\n') for x in res.text.split('///')])))]
        return res
    else:
        print(f"Error: Fail to extract info from {link}")
        return []


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--organisms", type=str, nargs='*', help="Multiple options from Fungi, Archaea, Bacteria", default=['Archaea','Bacteria', 'Fungi'])
    parser.add_argument("--log_dir", type=str, required=False, default=str('/'.join(path_list[:(project_index+1)] + ['log_folder'])), help="the path of logfile folder")
    parser.add_argument("--log_name", type=str, required=False, default="extract_kegg_organism_data.log", help="log file name")
    parser.add_argument("--outdir", type=str, help="The output dir")
    args = parser.parse_args()

    KEGG_api_link = 'http://rest.kegg.jp'

    ## set up logger
    logger = utils.get_logger(os.path.join(args.log_dir, args.log_name))
    args.organisms = [x.lower() for x in args.organisms]

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    ## download KEGG organism table
    if not os.path.exists(os.path.join(args.outdir,'organism_table.txt')):
        link = KEGG_api_link + '/list/organism'
        res = requests.get(link)
        if res.status_code == 200:
            logger.info('Organism table is sucessfully downloaded!')
            organism_table = pd.DataFrame([x.split('\t') for x in res.text.split('\n') if x.split('\t')[0]])
            organism_table.columns = ['T_number','org_code','name','lineage']

        ## extract taxa ids for each KEGG organism code
        org_code_list = list(organism_table['org_code'])
        org_code_taxaids = list(map(extract_taxid, org_code_list))
        params = zip(org_code_list, ['gb']*len(org_code_list))
        gb_ncbi_seq_ids = list(map(extract_ref_seq_id, params))
        params = zip(org_code_list, ['rs']*len(org_code_list))
        rs_ncbi_seq_ids = list(map(extract_ref_seq_id, params))
        organism_table = pd.concat([organism_table,pd.DataFrame(org_code_taxaids),pd.Series(gb_ncbi_seq_ids),pd.Series(rs_ncbi_seq_ids)], axis=1)
        organism_table.columns = ['T_number','org_code','name','lineage', 'taxaid','gb_ncbi_seq_ids','rs_ncbi_seq_ids']
        organism_table['org_code'] = organism_table['org_code'].astype(str)
        organism_table['taxaid'] = organism_table['taxaid'].astype(int).astype(str)
        organism_table.to_csv(os.path.join(args.outdir,'organism_table.txt'),sep='\t',index=None)
    else:
        organism_table = pd.read_csv(os.path.join(args.outdir,'organism_table.txt'), sep='\t', header=0)
        organism_table = organism_table.loc[~organism_table.taxaid.isna(),:].reset_index(drop=True)
        organism_table['org_code'] = organism_table['org_code'].astype(str)
        organism_table['taxaid'] = organism_table['taxaid'].astype(int).astype(str)

    organism_table = organism_table.loc[organism_table.apply(lambda row: len(set(args.organisms).intersection(set(row[3].lower().split(';')))) > 0, axis=1),:].reset_index(drop=True)

    ## extract gene/protein information from KEGG
    if not os.path.exists(os.path.join(args.outdir,'kegg_gene_info')):
        os.makedirs(os.path.join(args.outdir,'kegg_gene_info'))

    org_code_list = list(organism_table['org_code'])
    params = zip(org_code_list, [os.path.join(args.outdir,'kegg_gene_info')]*len(org_code_list))
    res = list(map(download_kegg_gene, params))

    all_gene_table_list = glob(os.path.join(args.outdir,'kegg_gene_info','*'))
    for index in trange(len(all_gene_table_list)):
        infile = pd.read_csv(all_gene_table_list[index], sep='\t', header=0)
        if 'koid' in infile.columns or 'koid_x' in infile.columns:
            if 'koid_x' in infile.columns:
                infile = infile[['kegg_gene_id','desc','koid_x','aaseq_x','ntseq_x']]
                infile.columns = ['kegg_gene_id','desc','koid','aaseq','ntseq']
                infile.to_csv(all_gene_table_list[index], sep='\t', index=None)
            continue

        kegg_gene_id_list = list(infile['kegg_gene_id'])
        batch =list(range(0,len(kegg_gene_id_list),10))
        batch.append(len(kegg_gene_id_list))
        kegg_gene_id_list = ['+'.join(kegg_gene_id_list[batch[i-1]:batch[i]]) for i in range(1, len(batch))]
        batch =list(range(0,len(kegg_gene_id_list),20))
        batch.append(len(kegg_gene_id_list))

        final_res = []
        get_data = True
        for i in range(1, len(batch)):
            start = batch[i-1]
            end = batch[i]
            while get_data:
                try:
                    with Pool(processes=20) as excutator:
                        res = excutator.map(process_query, kegg_gene_id_list[start:end])
                        get_data = False
                except:
                    time.sleep(300)
            final_res += [y for x in res for y in x]
            time.sleep(5)
            get_data = True
        final_res = pd.DataFrame(final_res)
        final_res.columns = ['kegg_gene_id','koid','aaseq','ntseq']
        outfile = infile.merge(final_res, on='kegg_gene_id', how='left').reset_index(drop=True)
        outfile.to_csv(all_gene_table_list[index], sep='\t', index=None)

    ## extract pathway information from KEGG
    if not os.path.exists(os.path.join(args.outdir,'kegg_pathway_info')):
        os.makedirs(os.path.join(args.outdir,'kegg_pathway_info'))

    org_code_list = list(organism_table['org_code'])
    org_code_list = list(set(org_code_list).difference(set([os.path.basename(path).split('_')[0] for path in glob(os.path.join(args.outdir,'kegg_pathway_info','*'))])))
    params = zip(org_code_list, [os.path.join(args.outdir,'kegg_pathway_info')]*len(org_code_list))
    res = list(map(download_kegg_pathway, params))

    ## link pathway information to gene information
    if not os.path.exists(os.path.join(args.outdir,'pathway_to_gene')):
        os.makedirs(os.path.join(args.outdir,'pathway_to_gene'))

    org_code_list = list(organism_table['org_code'])
    org_code_list = list(set(org_code_list).difference(set([os.path.basename(path).split('_')[0] for path in glob(os.path.join(args.outdir,'pathway_to_gene','*'))])))
    params = zip(org_code_list, [os.path.join(args.outdir,'pathway_to_gene')]*len(org_code_list))
    res = list(map(link_pathway_to_gene, params))

    ## link ko information to gene information
    if not os.path.exists(os.path.join(args.outdir,'ko_to_gene')):
        os.makedirs(os.path.join(args.outdir,'ko_to_gene'))

    org_code_list = list(organism_table['org_code'])
    org_code_list = list(set(org_code_list).difference(set([os.path.basename(path).split('_')[0] for path in glob(os.path.join(args.outdir,'ko_to_gene','*'))])))
    params = zip(org_code_list, [os.path.join(args.outdir,'ko_to_gene')]*len(org_code_list))
    res = list(map(link_ko_to_gene, params))

    ## link module information to gene information
    if not os.path.exists(os.path.join(args.outdir,'module_to_gene')):
        os.makedirs(os.path.join(args.outdir,'module_to_gene'))

    org_code_list = list(organism_table['org_code'])
    org_code_list = list(set(org_code_list).difference(set([os.path.basename(path).split('_')[0] for path in glob(os.path.join(args.outdir,'module_to_gene','*'))])))
    params = zip(org_code_list, [os.path.join(args.outdir,'module_to_gene')]*len(org_code_list))
    res = list(map(link_module_to_gene, params))

    ## link enzyme information to gene information
    if not os.path.exists(os.path.join(args.outdir,'enzyme_to_gene')):
        os.makedirs(os.path.join(args.outdir,'enzyme_to_gene'))

    org_code_list = list(organism_table['org_code'])
    org_code_list = list(set(org_code_list).difference(set([os.path.basename(path).split('_')[0] for path in glob(os.path.join(args.outdir,'enzyme_to_gene','*'))])))
    params = zip(org_code_list, [os.path.join(args.outdir,'enzyme_to_gene')]*len(org_code_list))
    res = list(map(link_enzyme_to_gene, params))

    ## convert kegg gene ids to uniprot ids
    if not os.path.exists(os.path.join(args.outdir,'convert_geneid_to_uniprot')):
        os.makedirs(os.path.join(args.outdir,'convert_geneid_to_uniprot'))

    org_code_list = list(organism_table['org_code'])
    org_code_list = list(set(org_code_list).difference(set([os.path.basename(path).split('_')[0] for path in glob(os.path.join(args.outdir,'convert_geneid_to_uniprot','*'))])))
    params = zip(org_code_list, [os.path.join(args.outdir,'convert_geneid_to_uniprot')]*len(org_code_list))
    res = list(map(convert_geneid_to_uniprot, params))