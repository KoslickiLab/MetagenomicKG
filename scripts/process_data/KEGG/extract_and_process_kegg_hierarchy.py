import os, sys
project_name = os.environ['project_name']
path_list = os.path.abspath(__file__).split('/')
project_index = path_list.index(project_name)
path_to_scripts = '/'.join(path_list[:(project_index+1)] + ['scripts'])
sys.path.append(path_to_scripts)
from utility_scripts import utils
import argparse
import logging
import pickle
import pandas as pd
import requests
import time
import re
from tqdm import tqdm

def get_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s  [%(levelname)s]  %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger

class KeggHierarchy(object):
    """
        This is a class object to extract and process the hierarchy structure of KEGG
    """

    def __init__(self, args, logger=None):
        self.KEGG_api_link = 'http://rest.kegg.jp'
        if logger is None:
            self.logger = get_logger()
        else:
            self.logger = logger
        self.args = args

    def get_request_link(self, operation, obj1, obj2=None):
        if operation == "list":
            if obj1 is None:
                self.logger.error("Please specify the parameter 'obj1'")
                return None
            else:
                return f"{self.KEGG_api_link}/{operation}/{obj1}"
        elif operation == "link":
            if obj1 is None or obj2 is None:
                self.logger.error("Please specify the parameters 'obj1' abd 'obj2'")
                return None
            else:
                return f"{self.KEGG_api_link}/{operation}/{obj1}/{obj2}"
        elif operation == "get":
            if obj1 is None:
                self.logger.error("Please specify the parameter 'obj1'")
                return None
            else:
                return f"{self.KEGG_api_link}/{operation}/{obj1}/json"
        else:
            self.logger.error("Currently this class only supports 'list' and 'link'")
            return None

    def run(self, discarded_brite_id_list=['br:br08902','br:br08904']):
        ## download KEGG brite information
        link = self.get_request_link(operation='list',obj1='brite')
        if link is not None:
            while True:
                try:
                    r = requests.get(link)
                    break
                except:
                    ## avoid rejection due to the frequent rquest
                    time.sleep(20)
            if r.status_code == 200:
                brite_table = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])
                brite_table.columns = ['kegg_brite_id','desc']
            else:
                self.logger.error(f"Fail to download KEGG brite information from {link}")
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

        # discarded_brite_id_list = ['br:br08902','br:br08904','br:br08330','br:br08332','br:br08411','br:br08410','br:br08420','br:br08601','br:br08610',
        #                         'br:br08611','br:br08602','br:br08603','br:br08620','br:br08621','br:br08605','br:br03220']
        ## filter out some kegg hierarchy structure
        hierarchy_dict = dict()
        hierarchy_raw_json = dict()

        brite_list = brite_table.values
        ncbi_brite_id_list = list(brite_table.loc[brite_table['desc'].str.contains('NCBI'),'kegg_brite_id'])

        ## read kg information
        kg = pd.read_csv(self.args.kg, sep='\t', header=0)
        nodes = list(set(list(kg['source'])+list(kg['target'])))
        taxa_nodes = [node for node in nodes if node.split(':')[0] == 'tax']

        for index, (brite_id, br_name) in enumerate(brite_list):
            brite_id, br_name = brite_list[index]
            self.logger.info(f"{index+1}/{len(brite_list)}: Extracting info from {brite_id} {br_name}")
            if brite_id not in discarded_brite_id_list:
                if brite_id not in ncbi_brite_id_list:
                    link = self.get_request_link(operation='get', obj1=brite_id)
                    while True:
                        try:
                            res = requests.get(link)
                            break
                        except:
                            ## avoid rejection due to the frequent rquest
                            time.sleep(20)
                    if res.status_code == 200:
                        if brite_id == 'br:br08901':
                            res_dict = self.organize_hierarchy1(res.json(), prefix='map')
                        elif brite_id == 'br:ko00001':
                            res_dict = self.organize_hierarchy1(res.json(), regex='^K\d{5} ')
                        elif brite_id == 'br:br08907':
                            res_dict = self.organize_hierarchy1(res.json(), regex='^N\d{5} ')
                        elif brite_id == 'br:br08202':
                            res_dict = self.organize_hierarchy1(res.json(), regex='^R\d{5}\t', stop_level=2)
                        # elif brite_id == 'br:br08204':
                        #     res_dict = self.organize_hierarchy1(res.json(), regex='^\d+\.\d+\.\d+\.\d+')
                        elif brite_id == 'br:ko01000':
                            res_dict = self.organize_hierarchy1(res.json(), regex='^([\d+-](\.)?)+ ')
                            res_dict = {ec_id:'|'.join([re.sub('\.$','',re.sub('[ ]+[^\d]*','',x)) for x in path.split('|')]) for ec_id, path in res_dict.items()}
                            kegg_enzymes = pd.read_csv(os.path.join(self.args.outdir,'kegg_enzymes.txt'), sep='\t',header=0)
                            res_dict_temp = {enzyme_id.replace('ec:',''):'|'.join(['ko01000','.'.join(enzyme_id.replace('ec:','').split('.')[:1]),'.'.join(enzyme_id.replace('ec:','').split('.')[:2]),'.'.join(enzyme_id.replace('ec:','').split('.')[:3])]) for enzyme_id in kegg_enzymes['kegg_enzyme_id']}
                            res_dict.update(res_dict_temp)
                        else:
                            res_dict = self.organize_hierarchy1(res.json())
                        hierarchy_raw_json[brite_id] = res.json()
                        for k, v in res_dict.items():
                            if k in hierarchy_dict:
                                hierarchy_dict[k] += [v]
                            else:
                                hierarchy_dict[k] = [v]
                    else:
                        self.logger.error(f"Fail to download hierarchical information from {link}")
                else:
                    # organism_virus_hierarchy_dict = dict()
                    # link = self.get_request_link(operation='get', obj1=brite_id)
                    # while True:
                    #     try:
                    #         res = requests.get(link)
                    #         break
                    #     except:
                    #         ## avoid rejection due to the frequent rquest
                    #         time.sleep(20)
                    # if res.status_code == 200:
                    #     res_dict =  self.organize_hierarchy2(res.json())
                    #     hierarchy_raw_json[brite_id] = res.json()
                    #     for key in res_dict:
                    #         if key in organism_virus_hierarchy_dict:
                    #             organism_virus_hierarchy_dict[key] += res_dict[key]
                    #         else:
                    #             organism_virus_hierarchy_dict[key] = res_dict[key]
                    # else:
                    #     self.logger.error(f"Fail to download hierarchical information from {link}")

                    ## process virus information
                    if brite_id == "br:br08620":
                        check_taxa_hierarchy = dict()
                        # virus_table = pd.read_csv(os.path.join(self.args.outdir,"kegg_viruses","virus_table.txt"), sep='\t', header=0)
                        # virus_table.loc[virus_table['rs_ncbi_seq_ids'].isnull(),'rs_ncbi_seq_ids'] = 'None'
                        # virus_table['rs_ncbi_seq_ids'] = virus_table['rs_ncbi_seq_ids'].apply(lambda row: eval(row) if type(row) is str else row)
                        # virus_table.loc[virus_table['gb_ncbi_seq_id'].isnull(),'gb_ncbi_seq_id'] = 'None'
                        # virus_table['gb_ncbi_seq_id'] = virus_table['gb_ncbi_seq_id'].apply(lambda row: eval(row) if type(row) is str else row)
                        # virus_table['taxaid'] = virus_table['taxaid'].astype(int)
                        # virus_table_with_seq = virus_table.loc[virus_table.apply(lambda row: row[2] is not None or row[3] is not None, axis=1),:].reset_index(drop=True)
                        # temp = set(virus_table_with_seq['taxaid'].astype(str))
                        # self.logger.info(f"Checking viruses' hierarchy")
                        # for node in temp:
                        #     if 'tax:'+node in taxa_nodes:
                        #         temp_path = organism_virus_hierarchy_dict.get(node)
                        #         if temp_path is not None:
                        #             if 'tax:'+node not in check_taxa_hierarchy:
                        #                 # check_taxa_hierarchy['tax:'+node] = ['|'.join(['root'] + ['Viruses'] + x.split('|')[1:]) for x in temp_path]
                        #                 check_taxa_hierarchy['tax:'+node] = temp_path
                        #         else:
                        #             self.logger.warning(f"Taxa ID {node} can't find the corresponding hierarchical info")
                        taxon_table = pd.read_csv(os.path.join(self.args.outdir,"all_taxon_info.txt"), sep='\t', header=None)
                        taxon_table.columns = ['orig_taxids', 'lineage', 'lineage_taxids', 'rank', 'lineage_ranks', 'species_taxids', 'selected_lineage', 'selected_lineage_taxids', 'selected_lineage_ranks']
                        taxon_table = taxon_table.loc[~taxon_table['selected_lineage'].isna(),:].reset_index(drop=True)
                        viruses_taxon_table = taxon_table.loc[taxon_table['selected_lineage'].str.contains('Viruses;'),:].reset_index(drop=True)
                        for index in range(len(viruses_taxon_table)):
                            node = viruses_taxon_table.loc[index,'species_taxids']
                            path_temp = viruses_taxon_table.loc[index,'selected_lineage']
                            check_taxa_hierarchy[f'tax:{node}'] = ['|'.join(['br08620']+ path_temp.split(';')[:-1])]
                        hierarchy_dict.update(check_taxa_hierarchy)

                    ## process organism information
                    if brite_id == "br:br08610":
                        check_taxa_hierarchy = dict()
                        # organism_table = pd.read_csv(os.path.join(self.args.outdir,"kegg_organisms","organism_table.txt"), sep='\t', header=0)
                        # organism_table = organism_table.loc[~organism_table.taxaid.isna(),:].reset_index(drop=True)
                        # organism_table['org_code'] = organism_table['org_code'].astype(str)
                        # organism_table['taxaid'] = organism_table['taxaid'].astype(int)
                        # organism_table['rs_ncbi_seq_ids'] = organism_table['rs_ncbi_seq_ids'].apply(lambda row: eval(row) if type(row) is str else row)
                        # organism_table['gb_ncbi_seq_id'] = organism_table['gb_ncbi_seq_id'].apply(lambda row: eval(row) if type(row) is str else row)
                        # organism_table = organism_table.loc[organism_table.apply(lambda row: len(set([x.lower() for x in self.args.organisms]).intersection(set(row[3].lower().split(';')))) > 0, axis=1),:].reset_index(drop=True)
                        # organism_table_with_seq = organism_table.loc[organism_table.apply(lambda row: len(row[5]) > 0 or len(row[6]) > 0, axis=1),:].reset_index(drop=True)
                        # # organism_table_with_seq = organism_table_with_seq.loc[organism_table_with_seq.apply(lambda row: 'tax:'+str(row[4]) in taxa_nodes, axis=1),:].reset_index(drop=True)
                        # self.logger.info(f"Checking organisms' hierarchy")
                        # for index in range(len(organism_table_with_seq)):
                        #     org_code, taxaid = organism_table_with_seq.loc[index,['org_code','taxaid']]
                        #     if 'tax:'+str(taxaid) in taxa_nodes:
                        #         temp_path = organism_virus_hierarchy_dict.get(org_code)
                        #         if temp_path is not None:
                        #             if 'tax:'+str(taxaid) not in check_taxa_hierarchy:
                        #                 # check_taxa_hierarchy['tax:'+str(taxaid)] = ['|'.join(['root'] + x.split('|')[1:]) for x in temp_path]
                        #                 check_taxa_hierarchy['tax:'+str(taxaid)] = temp_path
                        #         else:
                        #             self.logger.warning(f"Taxa ID {str(taxaid)} can't find the corresponding hierarchical info")
                        taxon_table = pd.read_csv(os.path.join(self.args.outdir,"all_taxon_info.txt"), sep='\t', header=None)
                        taxon_table.columns = ['orig_taxids', 'lineage', 'lineage_taxids', 'rank', 'lineage_ranks', 'species_taxids', 'selected_lineage', 'selected_lineage_taxids', 'selected_lineage_ranks']
                        taxon_table = taxon_table.loc[~taxon_table['selected_lineage'].isna(),:].reset_index(drop=True)
                        organisms_taxon_table = taxon_table.loc[~taxon_table['selected_lineage'].str.contains('Viruses;'),:].reset_index(drop=True)
                        for index in range(len(organisms_taxon_table)):
                            node = organisms_taxon_table.loc[index,'species_taxids']
                            path_temp = organisms_taxon_table.loc[index,'selected_lineage']
                            check_taxa_hierarchy[f'tax:{node}'] = ['|'.join(['br08610']+ path_temp.split(';')[:-1])]
                        hierarchy_dict.update(check_taxa_hierarchy)
            else:
                self.logger.warning(f"Skip {brite_id}\t{br_name}")

        return [brite_table, hierarchy_raw_json, hierarchy_dict]

    def organize_hierarchy1(self, hiearchy_json, regex='^\w?\w?\d{5} ', prefix='', stop_level=None):

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
            res_dict = {prefix+string.split('|')[-1].split(' ')[0]:'|'.join(string.split('|')[1:-1]) for string in res_list}
            return res_dict
        else:
            res_dict = {prefix+string.split('|')[-1].split(' ')[0].split('\t')[0]:'|'.join(string.split('|')[1:(stop_level+1)]) for string in res_list}
            return res_dict

    def organize_hierarchy2(self, hiearchy_json):

        def _iterate_multidimensional(res, hiearchy_json, res_list):
            if isinstance(hiearchy_json,dict):
                for k,v in hiearchy_json.items():
                    if k == 'name':
                        temp_name = re.sub(' \(.*\)','',re.sub(' \[.*\]','',hiearchy_json['name']))
                        res += f"|{temp_name}"
                        if len(temp_name.split('  ')) > 1:
                            res_list += [res]
                    elif k == 'children':
                        for elem in hiearchy_json['children']:
                            _iterate_multidimensional(res, elem, res_list)
            else:
                self.logger.error(f"{hiearchy_json} is not dictionary")
                raise

        res_list = []
        _iterate_multidimensional('', hiearchy_json, res_list)
        res_dict = dict()
        for string in res_list:
            part_list = string.split('|')
            key = part_list[-1].split('  ')[0]
            temp_val = part_list[-1].split('  ')[1]
            if part_list[-2] == temp_val:
                val = '|'.join(part_list[1:-1])
            else:
                val = '|'.join(part_list[1:-1] + [temp_val])
            if key in res_dict:
                res_dict[key] += [val]
            else:
                res_dict[key] = [val]
        return res_dict

    def convert_to_hierarchy_pairs(self, brite_table, hierarchy_dict, weight=1):
        nodes = list(pd.read_csv(os.path.join(self.args.outdir,'KG_info','nodes.txt'), sep='\t', header=0)['node_name'])

        node_mapping = {}
        for node in nodes:
            if node.split(':')[0] == 'tax':
                node_mapping[node] = node
            else:
                node_mapping[node.split(':')[1]] = node
        brite_to_des = {brite_id.split(':')[1]:des for brite_id, des in brite_table.values}
        type_mapping = {'cpd':'kegg cpd', 'ec':'kegg cpd', 'gl':'kegg gl', 'ko': 'kegg ko', 'md': 'kegg md', 'ne':'kegg ne', 'path':'kegg path', 'rn':'kegg rn', 'tax':'kegg tax'}
        hierarchy_dict = {node_mapping[key]:value for key, value in hierarchy_dict.items() if node_mapping.get(key)}

        hierarchy_pairs = []
        hierarchy_pairs += [(node, type_mapping[node.split(':')[0]], weight) for node in nodes]
        self.logger.info(f"generate hierarchical pairs")
        for node in tqdm(hierarchy_dict):
            if node.split(':')[0] != 'rn':
                paths = hierarchy_dict[node]
                for path in paths:
                    branch = path.split('|')
                    if branch[0] in ["br08610", "br08620"]:
                        branch[0] = "Taxonomy lineage"
                    else:
                        branch[0] = brite_to_des[branch[0]]
                    branch += [node]
                    for index in range(1, len(branch)):
                        branch[index-1], branch[index] = self.rename_brach_item(branch[index-1], branch[index])
                    branch.reverse()
                    b_len = len(branch)
                    hierarchy_pairs += [(branch[index1], branch[index2], weight) for index1 in range(b_len) for index2 in range(index1+1, b_len)]
            else:
                paths = hierarchy_dict[node]
                for path in paths:
                    branch = path.split('|')
                    if branch[-1] == '-.-.-.-':
                        continue
                    branch[0] = brite_to_des[branch[0]]
                    branch += [node]
                    for index in range(1, len(branch)):
                        branch[index-1], branch[index] = branch[index-1], re.sub('\.$','',branch[index].split(' ')[0])
                    branch.reverse()
                    b_len = len(branch)
                    hierarchy_pairs += [(branch[index1], branch[index2], weight) for index1 in range(b_len) for index2 in range(index1+1, b_len)]

        return list(set(hierarchy_pairs))

    @staticmethod
    def rename_brach_item(parent, child):
        if child == 'Other':
            child = f'Other sub {parent}'
        elif child == 'Others':
            child = f'Other sub {parent}'
        elif child == 'Unclassified':
            child = f'Unclassified sub {parent}'
        elif child == 'Not specified':
            child = f'Not specified sub {parent}'
        else:
            pass
        return parent, child


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--kg", type=str, help="The full path of original kg")
    parser.add_argument("--organisms", type=str, nargs='*', help="Multiple options from Fungi, Archaea, Bacteria", default=['Archaea','Bacteria', 'Fungi'])
    parser.add_argument("--log_dir", type=str, required=False, default=str('/'.join(path_list[:(project_index+1)] + ['log_folder'])), help="the path of logfile folder")
    parser.add_argument("--log_name", type=str, required=False, default="extract_and_process_kegg_hierarchy.log", help="log file name")
    parser.add_argument("--outdir", type=str, help="The output dir")
    args = parser.parse_args()


    ## set up logger
    logger = utils.get_logger(os.path.join(args.log_dir, args.log_name))

    exclude_list = [
        "br:br08902\tBRITE hierarchy files",
        "br:br08904\tBRITE table files",
        "br:ko00002\tKEGG modules",
        "br:br08303\tAnatomical Therapeutic Chemical (ATC) classification",
        "br:br08302\tUSP drug classification",
        "br:br08301\tTherapeutic category of drugs in Japan",
        "br:br08313\tClassification of Japanese OTC drugs",
        "br:br08312\tRisk category of Japanese OTC drugs",
        "br:br08304\tTraditional Chinese Medicine in Japan",
        "br:br08305\tCrude drugs",
        "br:br08331\tAnimal drugs in Japan",
        "br:br08330\tDrug groups",
        "br:br08332\tDrug classes",
        "br:br08310\tTarget-based classification of drugs",
        "br:br08307\tAntimicrobials",
        "br:br08327\tAntimicrobials abbreviations",
        "br:br08311\tDrugs listed in the Japanese Pharmacopoeia",
        "br:br08402\tHuman diseases",
        "br:br08401\tInfectious diseases",
        "br:br08403\tHuman diseases in ICD-11 classification",
        "br:br08411\tICD-11 International Classification of Diseases",
        "br:br08410\tICD-10 International Classification of Diseases",
        "br:br08420\tICD-O-3 International Classification of Diseases for Oncology",
        "br:br08601\tKEGG organisms",
        "br:br08611\tKEGG organisms in taxonomic groups",
        "br:br08612\tKEGG organisms: animals",
        "br:br08613\tKEGG organisms: plants",
        "br:br08614\tKEGG organisms: fungi",
        "br:br08615\tKEGG organisms: protists",
        "br:br08621\tKEGG viruses in taxonomic groups",
        "br:br08605\tPlant pathogens",
        "br:br03220\tEnveloped virus entry",
        "br:br08319\tNew drug approvals in the USA",
        "br:br08329\tNew drug approvals in Europe",
        "br:br08318\tNew drug approvals in Japan",
        "br:br08328\tNew drug approvals in the USA, Europe and Japan",
        "br:br08309\tDrug metabolizing enzymes and transporters",
        "br:br08341\tPharmacogenomic biomarkers",
        "br:br08324\tProdrugs",
        "br:br08317\tTopical steroids potency",
        "br:br08315\tRx-to-OTC switch list in the USA",
        "br:br08314\tRx-to-OTC switch list in Japan",
        "br:br08442\tTumor markers",
        "br:br08441\tCancer-associated carbohydrates",
        "br:br08431\tCarbohydrates in viral and bacterial infections"
    ]

    kegghierarchy = KeggHierarchy(args, logger)
    brite_table, hierarchy_raw_json, hierarchy_dict = kegghierarchy.run(discarded_brite_id_list=[x.split('\t')[0] for x in exclude_list])
    hierarchy_pairs = kegghierarchy.convert_to_hierarchy_pairs(brite_table, hierarchy_dict)

    with open(os.path.join(args.outdir,"kegg_brite_info.pkl"),'wb') as outfile:
        pickle.dump([brite_table,hierarchy_raw_json,hierarchy_dict,hierarchy_pairs],outfile)

    hierarchical_edge_pairs = pd.DataFrame(hierarchy_pairs, columns=['offspring','ancestor','weight'])
    hierarchical_edge_pairs.to_csv(os.path.join(args.outdir,"hierarchical_edge_pairs.txt"), sep='\t', index=None)

# |Other -> |Other sub-parent
# |Others -> |Other sub-parent
# |Unclassified -> |Unclassified sub-parent
# |Not specified -> |Not specified sub-parent
