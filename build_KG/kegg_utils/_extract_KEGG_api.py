## Import standard libraries
import os
import sys
import pandas as pd
import requests

## Import custom libraries
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/..')
from utils import get_logger

class GetKeggLinkData(object):
    """This class is used to download KEGG link data via KEGG API."""

    def __init__(self):
        self.KEGG_api_link = 'http://rest.kegg.jp'
        self.logger = get_logger()
        self.session = requests.Session()

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
        else:
            self.logger.error("Currently this class only supports 'list' and 'link'")
            return None

    def download_kegg_compound(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='list',obj1='compound')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly download compound information from {link}")
                table = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])[[0,1]]
                table.columns = ['kegg_compound_id','desc']
            else:
                self.logger.error(f"Fail to download compound information from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"kegg_compounds.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def download_kegg_pathway(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='list',obj1='pathway')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly download pathway information from {link}")
                table = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])[[0,1]]
                table.columns = ['kegg_pathway_id','desc']
            else:
                self.logger.error(f"Fail to download pathway information from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"kegg_pathways.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def download_kegg_module(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='list',obj1='module')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly download module information from {link}")
                table = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])[[0,1]]
                table.columns = ['kegg_module_id','desc']
            else:
                self.logger.error(f"Fail to download module information from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"kegg_modules.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def download_kegg_glycan(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='list',obj1='glycan')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly download glycan information from {link}")
                table = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])[[0,1]]
                table.columns = ['kegg_glycan_id','desc']
            else:
                self.logger.error(f"Fail to download glycan information from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"kegg_glycans.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def download_kegg_reaction(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='list',obj1='reaction')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly download reaction information from {link}")
                table = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])[[0,1]]
                table.columns = ['kegg_reaction_id','desc']
            else:
                self.logger.error(f"Fail to download reaction information from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"kegg_reactions.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def download_kegg_enzyme(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='list',obj1='enzyme')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly download enzyme information from {link}")
                table = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])[[0,1]]
                table.columns = ['kegg_enzyme_id','desc']
            else:
                self.logger.error(f"Fail to download enzyme information from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"kegg_enzymes.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def download_kegg_network(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='list',obj1='network')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly download network information from {link}")
                table = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])[[0,1]]
                table.columns = ['kegg_network_id','desc']
            else:
                self.logger.error(f"Fail to download network information from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"kegg_networks.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def download_kegg_ko(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='list',obj1='ko')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly download ko information from {link}")
                table = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])[[0,1]]
                table.columns = ['kegg_ko_id','desc']
            else:
                self.logger.error(f"Fail to download ko information from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"kegg_koids.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def download_kegg_disease(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='list',obj1='disease')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly download disease information from {link}")
                table = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])[[0,1]]
                table.columns = ['kegg_disease_id','desc']
            else:
                self.logger.error(f"Fail to download disease information from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"kegg_diseases.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def download_kegg_drug(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='list',obj1='drug')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly download drug information from {link}")
                table = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])[[0,1]]
                table.columns = ['kegg_drug_id','desc']
            else:
                self.logger.error(f"Fail to download drug information from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"kegg_drugs.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def download_kegg_gn(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='list',obj1='gn')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly download genome information from {link}")
                table = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])[[0,1]]
                table.columns = ['kegg_gn_id','desc']
            else:
                self.logger.error(f"Fail to download genome information from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"kegg_gns.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def download_kegg_rclass(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='list',obj1='rclass')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly download reaction class information from {link}")
                table = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])[[0,1]]
                table.columns = ['kegg_rclass_id','desc']
            else:
                self.logger.error(f"Fail to download reaction class information from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"kegg_rclasses.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def download_kegg_dgroup(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='list',obj1='dgroup')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly download drug class information from {link}")
                table = pd.DataFrame([x.split('\t') for x in r.text.split('\n') if x.split('\t')[0]])[[0,1]]
                table.columns = ['kegg_dgroup_id','desc']
            else:
                self.logger.error(f"Fail to download drug class information from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"kegg_dgroups.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_compound_to_pathway(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='compound',obj2='pathway')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link compounds to pathways from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:has_participant','biolink:participates_in'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_pathway_id','kegg_compound_id','source_to_target','target_to_source']
                table['kegg_pathway_id'] = table['kegg_pathway_id'].str.replace('path:[a-z]*','path:map',regex=True)
                table = table.drop_duplicates().reset_index(drop=True)
            else:
                self.logger.error(f"Fail to link compounds to pathways from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_compound_to_pathway.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_compound_to_module(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='compound',obj2='module')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link compounds to modules from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:has_participant','N/A'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_module_id','kegg_compound_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link compounds to modules from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_compound_to_module.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_compound_to_glycan(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='compound',obj2='glycan')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link compounds to glycans from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:physically_interacts_with','biolink:physically_interacts_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_glycan_id','kegg_compound_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link compounds to glycans from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_compound_to_glycan.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_compound_to_reaction(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='compound',obj2='reaction')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link compounds to reactions from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:has_participant','biolink:participates_in']  for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_reaction_id','kegg_compound_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link compounds to reactions from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_compound_to_reaction.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_compound_to_enzyme(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='compound',obj2='enzyme')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link compounds to enzymes from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:physically_interacts_with','biolink:physically_interacts_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_enzyme_id','kegg_compound_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link compounds to enzymes from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_compound_to_enzyme.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_compound_to_network(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='compound',obj2='network')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link compounds to networks from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:has_participant','biolink:participates_in'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_network_id','kegg_compound_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link compounds to networks from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_compound_to_network.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_compound_to_drug(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='compound',obj2='drug')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link compounds to drugs from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:same_as','biolink:same_as'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_drug_id','kegg_compound_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link compounds to drugs from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_compound_to_drug.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")

    def link_compound_to_gn(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='compound',obj2='gn')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link compounds to genomes from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:produces','N/A'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_gn_id','kegg_compound_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link compounds to genomes from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_compound_to_gn.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")

    def link_pathway_to_module(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='pathway',obj2='module')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link pathways to modules from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:participates_in','biolink:has_participant'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_module_id','kegg_pathway_id','source_to_target','target_to_source']
                table['kegg_pathway_id'] = table['kegg_pathway_id'].str.replace('path:[a-z]*','path:map',regex=True)
                table = table.drop_duplicates().reset_index(drop=True)
            else:
                self.logger.error(f"Fail to link pathways to modules from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_pathway_to_module.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_pathway_to_glycan(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='pathway',obj2='glycan')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link pathways to glycans from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:participates_in','biolink:has_participant']  for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_glycan_id','kegg_pathway_id','source_to_target','target_to_source']
                table['kegg_pathway_id'] = table['kegg_pathway_id'].str.replace('path:[a-z]*','path:map',regex=True)
                table = table.drop_duplicates().reset_index(drop=True)
            else:
                self.logger.error(f"Fail to link pathways to glycans from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_pathway_to_glycan.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_pathway_to_reaction(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='pathway',obj2='reaction')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link pathways to reactions from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:participates_in','biolink:has_participant'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_reaction_id','kegg_pathway_id','source_to_target','target_to_source']
                table['kegg_pathway_id'] = table['kegg_pathway_id'].str.replace('path:[a-z]*','path:map',regex=True)
                table = table.drop_duplicates().reset_index(drop=True)
            else:
                self.logger.error(f"Fail to link pathways to reactions from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_pathway_to_reaction.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_pathway_to_enzyme(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='pathway',obj2='enzyme')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link pathways to enzymes from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:participates_in','biolink:has_participant'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_enzyme_id','kegg_pathway_id','source_to_target','target_to_source']
                table['kegg_pathway_id'] = table['kegg_pathway_id'].str.replace('path:[a-z]*','path:map',regex=True)
                table = table.drop_duplicates().reset_index(drop=True)
            else:
                self.logger.error(f"Fail to link pathways to enzymes from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_pathway_to_enzyme.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_pathway_to_network(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='pathway',obj2='network')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link pathways to networks from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:has_participant', 'biolink:participates_in'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_network_id','kegg_pathway_id','source_to_target','target_to_source']
                table['kegg_pathway_id'] = table['kegg_pathway_id'].str.replace('path:[a-z]*','path:map',regex=True)
                table = table.drop_duplicates().reset_index(drop=True)
            else:
                self.logger.error(f"Fail to link pathways to networks from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_pathway_to_network.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_pathway_to_ko(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='pathway',obj2='ko')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link pathways to koids from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:participates_in','biolink:has_participant'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_ko_id','kegg_pathway_id','source_to_target','target_to_source']
                table['kegg_pathway_id'] = table['kegg_pathway_id'].str.replace('path:[a-z]*','path:map',regex=True)
                table = table.drop_duplicates().reset_index(drop=True)
            else:
                self.logger.error(f"Fail to link pathways to koids from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_pathway_to_ko.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_pathway_to_disease(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='pathway',obj2='disease')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link pathways to diseases from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:associated_with','biolink:associated_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_disease_id','kegg_pathway_id','source_to_target','target_to_source']
                table['kegg_pathway_id'] = table['kegg_pathway_id'].str.replace('path:[a-z]*','path:map',regex=True)
                table = table.drop_duplicates().reset_index(drop=True)
            else:
                self.logger.error(f"Fail to link pathways to diseases from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_pathway_to_disease.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_pathway_to_drug(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='pathway',obj2='drug')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link pathways to drugs from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:associated_with','biolink:associated_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_drug_id','kegg_pathway_id','source_to_target','target_to_source']
                table['kegg_pathway_id'] = table['kegg_pathway_id'].str.replace('path:[a-z]*','path:map',regex=True)
                table = table.drop_duplicates().reset_index(drop=True)
            else:
                self.logger.error(f"Fail to link pathways to drugs from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_pathway_to_drug.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_pathway_to_gn(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='pathway',obj2='gn')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link pathways to genomes from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:associated_with','biolink:associated_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_gn_id','kegg_pathway_id','source_to_target','target_to_source']
                table['kegg_pathway_id'] = table['kegg_pathway_id'].str.replace('path:[a-z]*','path:map',regex=True)
                table = table.drop_duplicates().reset_index(drop=True)
            else:
                self.logger.error(f"Fail to link pathways to genomes from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_pathway_to_gn.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_pathway_to_rclass(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='pathway',obj2='rclass')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link pathways to reaction classes from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:participates_in','biolink:has_participant'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_rclass_id','kegg_pathway_id','source_to_target','target_to_source']
                table['kegg_pathway_id'] = table['kegg_pathway_id'].str.replace('path:[a-z]*','path:map',regex=True)
                table = table.drop_duplicates().reset_index(drop=True)
            else:
                self.logger.error(f"Fail to link pathways to reaction classes from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_pathway_to_rclass.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_module_to_glycan(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='module',obj2='glycan')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link modules to glycans from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:participates_in','biolink:has_participant'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_glycan_id','kegg_module_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link modules to glycans from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_module_to_glycan.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_module_to_reaction(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='module',obj2='reaction')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link modules to reactions from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:participates_in','biolink:has_participant'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_reaction_id','kegg_module_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link modules to reactions from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_module_to_reaction.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_module_to_enzyme(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='module',obj2='enzyme')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link modules to enzymes from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:participates_in','biolink:has_participant'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_enzyme_id','kegg_module_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link modules to enzymes from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_module_to_enzyme.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_module_to_ko(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='module',obj2='ko')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link modules to koids from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:participates_in','biolink:has_participant'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_ko_id','kegg_module_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link modules to koids from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_module_to_ko.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_module_to_disease(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='module',obj2='disease')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link modules to diseases from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:associated_with','biolink:associated_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_disease_id','kegg_module_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link modules to diseases from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_module_to_disease.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_module_to_gn(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='module',obj2='gn')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link modules to genomes from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:associated_with','biolink:associated_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_gn_id','kegg_module_id','source_to_target','target_to_source']
                table['kegg_module_id'] = table['kegg_module_id'].str.replace(':[a-z]*_',':',regex=True)
                table = table.drop_duplicates().reset_index(drop=True)
            else:
                self.logger.error(f"Fail to link modules to genomes from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_module_to_gn.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_glycan_to_reaction(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='glycan',obj2='reaction')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link glycans to reactions from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:has_participant','biolink:participates_in'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_reaction_id','kegg_glycan_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link glycans to reactions from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_glycan_to_reaction.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_glycan_to_enzyme(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='glycan',obj2='enzyme')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link glycans to enzymes from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:physically_interacts_with','biolink:physically_interacts_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_enzyme_id','kegg_glycan_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link glycans to enzymes from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_glycan_to_enzyme.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_reaction_to_enzyme(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='reaction',obj2='enzyme')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link reactions to enzymes from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:participates_in','biolink:has_participant'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_enzyme_id','kegg_reaction_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link reactions to enzymes from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_reaction_to_enzyme.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_reaction_to_ko(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='reaction',obj2='ko')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link reactions to koids from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:has_participant', 'biolink:participates_in'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_ko_id','kegg_reaction_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link reactions to koids from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_reaction_to_ko.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_reaction_to_rclass(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='reaction',obj2='rclass')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link reactions to rclasses from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:superclass_of','biolink:subclass_of'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_rclass_id','kegg_reaction_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link reactions to rclasses from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_reaction_to_rclass.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_enzyme_to_ko(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='enzyme',obj2='ko')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link enzymes to koids from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:associated_with','biolink:associated_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_ko_id','kegg_enzyme_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link enzymes to koids from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_enzyme_to_ko.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_enzyme_to_rclass(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='enzyme',obj2='rclass')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link enzymes to rclasses from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:associated_with','biolink:associated_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_rclass_id','kegg_enzyme_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link enzymes to rclasses from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_enzyme_to_rclass.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_network_to_ko(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='network',obj2='ko')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link networks to koids from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:associated_with','biolink:associated_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_ko_id','kegg_network_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link networks to koids from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_network_to_ko.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_network_to_disease(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='network',obj2='disease')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link networks to diseases from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:associated_with','biolink:associated_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_disease_id','kegg_network_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link networks to diseases from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_network_to_disease.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_ko_to_disease(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='ko',obj2='disease')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link koids to diseases from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:associated_with','biolink:associated_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_disease_id','kegg_ko_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link koids to diseases from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_ko_to_disease.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_ko_to_drug(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='ko',obj2='drug')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link koids to drugs from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:associated_with','biolink:associated_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_drug_id','kegg_ko_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link koids to drugs from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_ko_to_drug.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_ko_to_rclass(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='ko',obj2='rclass')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link koids to rclasses from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:associated_with','biolink:associated_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_rclass_id','kegg_ko_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link koids to rclasses from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_ko_to_rclass.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_disease_to_drug(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='disease',obj2='drug')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link diseases to drugs from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:treats','N/A'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_drug_id','kegg_disease_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link diseases to drugs from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_disease_to_drug.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_disease_to_gn(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='disease',obj2='gn')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link diseases to genomes from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:associated_with','biolink:associated_with'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_gn_id','kegg_disease_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link diseases to genomes from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_disease_to_gn.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1

    def link_drug_to_dgroup(self, out_loc=None):
        if out_loc is None:
            self.logger.error("Please specify the parameters 'out_loc'")
            return 1
        else:
            if not os.path.exists(out_loc):
                os.makedirs(out_loc)
        link = self.get_request_link(operation='link',obj1='drug',obj2='dgroup')
        if link is not None:
            r = self.session.get(link)
            if r.status_code == 200:
                self.logger.info(f"Successuflly link drugs to dgroups from {link}")
                table = pd.DataFrame([x.split('\t')+['biolink:superclass_of','biolink:subclass_of'] for x in r.text.split('\n') if x.split('\t')[0]])
                table.columns = ['kegg_dgroup_id','kegg_drug_id','source_to_target','target_to_source']
            else:
                self.logger.error(f"Fail to link drugs to dgroups from {link}")
                return 0
            table.to_csv(os.path.join(out_loc,f"link_drug_to_dgroup.txt"), sep='\t', index=None)
            return 1
        else:
            self.logger.error(f"There are something wrong to get request link")
            return 1


