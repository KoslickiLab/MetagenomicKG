import os
import sys
from multiprocessing import Pool, cpu_count
import pandas as pd
from collections import Counter
from tqdm import tqdm, trange
import pickle
import argparse
import logging
from Bio import Entrez
Entrez.email = "test@example.com"

def get_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s  [%(levelname)s]  %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger

def download_seq(seq_id):

    # Downloading...
    try:
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
    except:
        print(f"Error: Fail to donwload nucleotide sequence for {seq_id}", flush=True)
        return None
    seq = handle.read().replace('\n\n','\n')
    handle.close()
    return seq

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", type=str, help="The full path of virus/organism table")
    parser.add_argument("--organisms", type=str, nargs='*', help="Multiple options from Fungi, Archaea, Bacteria", default=None)
    parser.add_argument("--col", type=str, help="Download seqs based on the ids from specific column", default='rs_ncbi_seq_ids')
    parser.add_argument("--outfile", type=str, help="The full path of output file")
    args = parser.parse_args()

    logger = get_logger()

    if args.organisms is not None:
        args.organisms = [x.lower() for x in args.organisms]

    ## read virus/organism table
    logger.info("Read virus/organism table")
    table = pd.read_csv(args.table, sep='\t', header=0)
    table.loc[table['rs_ncbi_seq_ids'].isnull(),'rs_ncbi_seq_ids'] = 'None'
    table['rs_ncbi_seq_ids'] = table['rs_ncbi_seq_ids'].apply(lambda row: eval(row) if type(row) is str else row)
    table.loc[table['gb_ncbi_seq_ids'].isnull(),'gb_ncbi_seq_ids'] = 'None'
    table['gb_ncbi_seq_ids'] = table['gb_ncbi_seq_ids'].apply(lambda row: eval(row) if type(row) is str else row)
    table['taxaid'] = table['taxaid'].astype(int)

    if args.organisms is not None:
        table = table.loc[table.apply(lambda row: len(set(args.organisms).intersection(set(row[3].lower().split(';')))) > 0, axis=1),:].reset_index(drop=True)

    ## start to download sequences
    logger.info("Start to download sequences")
    seq_id_list = [y for x in table[args.col] if x is not None for y in x ]
    seqs = list(map(download_seq, tqdm(seq_id_list)))

    ## save sequences
    with open(args.outfile,'w') as out_handle:
        for seq in seqs:
            if seq is not None:
                out_handle.write(seq)


