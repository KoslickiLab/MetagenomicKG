## Import standard libraries
import os
import sys
from tqdm import tqdm, trange
from glob import glob
import pandas as pd
import yaml

## Import custom libraries
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/../../build_KG/')
from kegg_utils.extract_KEGG_data import KEGGData

## Load configuration file
with open(f'{os.path.dirname(os.path.realpath(__file__))}/../../config.yml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)

## Set up paths
kegg_data_dir = config['BUILD_KG_VARIABLES']['KEGG_FTP_DATA_DIR']
output_dir = './'

# Create KEGGData object
keggdata = KEGGData(kegg_data_dir, output_dir)

# Create a regular expression pattern from the microbe labels
microbe_labels = ['Archaea', 'Bacteria']
pattern = '|'.join(microbe_labels)
microbe_gn_table = keggdata.all_gn_table.loc[keggdata.all_gn_table['kegg_lineage'].str.contains(pattern, regex=True)].reset_index(drop=True)
microbe_gn_table[['gn_id','assembly_id']].to_csv('KEGG_microbe_assembly_mapping.tsv', sep='\t', index=False)