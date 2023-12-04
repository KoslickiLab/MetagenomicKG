## Import standard libraries
import os
import sys
import gzip
from tqdm import tqdm, trange
from glob import glob
import pandas as pd
import tarfile
import io
import argparse
import re
import requests
import pytaxonkit
import gzip
from typing import List, Dict, Tuple, Union, Any, Optional
from Bio import Entrez
import logging
Entrez.email = "test@example.com"

## Import custom libraries
from utils import get_logger
from _extract_KEGG_api import GetKeggLinkData

class KEGGData:
    def __init__(self, kegg_data_dir: str, output_dir: str):
        
        self.logger = get_logger()
        self.kegg_data_dir = kegg_data_dir
        self.logger.info("generating organism table")
        self.organism_table = self._read_KEGG_organism_taxonomy()
        self.logger.info("generating virus table")
        self.virus_table = self._read_KEGG_virus_taxonomy()
        self.logger.info("generating all genome table")
        self.all_gn_table = self._get_all_KEGG_genome_info()
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    @staticmethod
    def _extract_info(org_code: str, entry: str):
        """
        Extract information from a KEGG entry
        """
        entry_lines = entry.split("\n")
        entry_id = ""
        symbol = ""
        desc = ""
        aaseq = ""
        ntseq = ""

        index = 0
        max_index = len(entry_lines)
        while index < max_index:
            line = entry_lines[index]
            if line.startswith("ENTRY"):
                entry_id = org_code+':'+line.split()[1]
            elif line.startswith("SYMBOL"):
                symbol = line.split(" ", 1)[1].strip()
            elif line.startswith("NAME"):
                desc = line.split(" ", 1)[1].strip()
            elif line.startswith("AASEQ"):
                index += 1
                aaseq = ""
                while not entry_lines[index].startswith("NTSEQ"):
                    aaseq += entry_lines[index].strip()
                    index += 1
                index -= 1
            elif line.startswith("NTSEQ"):
                index += 1
                ntseq = ""
                while index < max_index and len(entry_lines[index]) != 0:
                    ntseq += entry_lines[index].strip()
                    index += 1
                index -= 1 
            index += 1

        return [entry_id, symbol, desc, aaseq, ntseq]

    @staticmethod
    def extract_dblinks(entry: str, prefix: str):

        def _parse_ids(line: str):
            line = line.split(' ')
            return [f"{line[0]}{line[index]}" for index in range(1, len(line))]

        entry_lines = entry.split("\n")
        # Initialize the variables
        entry_id = ""
        dblink_ids = []

        # Iterate through the lines
        index = 0
        max_index = len(entry_lines)
        while index < max_index:
            line = entry_lines[index]
            if line.startswith("ENTRY"):
                entry_id = prefix+':'+line.split()[1]
            elif line.startswith("DBLINKS"):
                dblink_ids += _parse_ids(line.replace("DBLINKS", "").strip())
                index += 1
                while entry_lines[index].startswith(" "):
                    dblink_ids += _parse_ids(entry_lines[index].strip())
                    index += 1
                break
            index += 1

        return entry_id, dblink_ids

    def _read_KEGG_organism_taxonomy(self):
        """
        Read KEGG FTP organism taxonomy files
        """
        taxonomy_file_path = os.path.join(self.kegg_data_dir, "genes", "misc", "taxonomy")
        if not os.path.exists(taxonomy_file_path):
            self.logger.warning(f"File {taxonomy_file_path} does not exist!")
            return None
        
        taxonomy_rank_file_path = os.path.join(self.kegg_data_dir, "genes", "misc", "taxonomic_rank")
        if not os.path.exists(taxonomy_rank_file_path):
            self.logger.warning(f"File {taxonomy_rank_file_path} does not exist!")
            return None

        # Open the file for reading
        with open(taxonomy_file_path, "r") as file_in:
            # Initialize an empty list to store the output table rows
            output_table = []
            # Initialize a list to store the current hierarchy levels
            hierarchy = []
            # Iterate through each line in the input text
            for line in file_in:
                line = line.strip()
                # Skip empty lines
                if not line:
                    continue
                # If the line starts with '#', it's a hierarchy level
                if line.startswith('#'):
                    # Count the number of '#' characters to determine the hierarchy level
                    level = line.count('#') - 1
                    # Update the hierarchy list with the new level
                    hierarchy = hierarchy[:level] + [line.replace('#', '').strip()]
                else:
                    fields = line.split('\t')
                    org_id = fields[2]
                    org_code = fields[1]
                    org_desc = fields[3]
                    org_hierarchy = ';'.join(hierarchy)
                    output_table.append([org_id, org_code, org_hierarchy, org_desc])

        # Convert the output table to a pandas dataframe
        table1 = pd.DataFrame(output_table, columns=['gn_id', 'org_code', 'org_kegg_lineage', 'org_desc'])

        # Open the file for reading
        with open(taxonomy_rank_file_path, "r") as file_in:
            # Initialize an empty list to store the output table rows
            output_table = []
            # Iterate through each line in the input text
            for line in file_in:
                line = line.strip()
                # Skip empty lines
                if not line:
                    continue
                if not line.startswith('#'):
                    fields = line.split('\t')
                    org_code = fields[0]
                    org_rank = fields[1]
                    output_table.append([org_code, org_rank])

        # Convert the output table to a pandas dataframe
        table2 = pd.DataFrame(output_table, columns=['org_code', 'taxon_id'])

        # Merge the two dataframes
        table = table1.merge(table2, on='org_code')
        ## re-arrange columns
        table =  table[['gn_id', 'org_code', 'taxon_id', 'org_desc', 'org_kegg_lineage']].reset_index(drop=True)

        ## Return taxonomy dataframe
        return table


    def _read_KEGG_virus_taxonomy(self):
        """
        Read KEGG FTP virus taxonomy files
        """
        taxonomy_file_path = glob(os.path.join(self.kegg_data_dir, "genes", "viruses", "*.tax.gz"))[0]
        if not os.path.exists(taxonomy_file_path):
            self.logger.warning(f"File {taxonomy_file_path} does not exist!")
            return None
        else:
            # Read data
            with gzip.open(taxonomy_file_path, 'rb') as f:
                text = f.read()
                text = text.decode('utf-8').split("\n")

        # Convert the output table to a pandas dataframe
        table1 = pd.DataFrame([[row.split('\t')[0], row.split('\t')[1].replace('tax:',''), row.split('\t')[2]] for row in text if len(row) > 0], columns=['gene_id', 'taxon_id', 'virus_kegg_lineage'])

        gene_path = os.path.join(self.kegg_data_dir, "genes")
        in_path = os.path.join(gene_path, "viruses", "viruses_link.tar.gz")
        if not os.path.exists(in_path):
            self.logger.warning(f"File {in_path} does not exist!")
            return None
        else:
            # Open the compressed file for reading as binary data
            with open(in_path, 'rb') as f:
                # Decompress the file
                decompressed_data = gzip.decompress(f.read())
                # Wrap the decompressed data in an io.BytesIO object
                fileobj = io.BytesIO(decompressed_data)
                # Open the tar archive from the wrapped data
                with tarfile.open(fileobj=fileobj, mode='r') as tar:
                    gn_contents = tar.extractfile('vg_gn.list').read().decode('utf-8').split('\n')

        # Convert the output table to a pandas dataframe
        table2 = pd.DataFrame([list(row.split('\t')) for row in gn_contents if len(row) > 0], columns=['gene_id', 'gn_id'])

        # Merge the two dataframes
        table = table1.merge(table2, on='gene_id')
        ## re-arrange columns
        table =  table[['gene_id','gn_id','taxon_id','virus_kegg_lineage']].reset_index(drop=True)

        ## Return taxonomy dataframe
        return table


    def _get_all_KEGG_genome_info(self):
        """
        Get all KEGG genome information
        """
        def extract_info(entry: str):
            entry_lines = entry.split("\n")
            gn_id = ''
            org_code = ''
            desc = ''
            taxon_id = ''
            key_words = ''
            lineage =''
            sequence_ids = []
            assembly_id = ''

            # Iterate through the lines
            for line in entry_lines:
                if line.startswith("ENTRY"):
                    gn_id = line.split()[1]
                elif line.startswith("ORG_CODE"):
                    org_code = line.split()[1]
                elif line.startswith("NAME"):
                    desc = line.split("NAME")[1].strip()
                elif line.startswith("TAXONOMY"):
                    taxon_id = line.split("TAX:")[1].split()[0]
                elif line.startswith("KEYWORDS"):
                    key_words = line.split("KEYWORDS")[1].strip()
                elif line.startswith("  LINEAGE"):
                    lineage = line.split("LINEAGE")[1].strip()
                elif line.startswith("  SEQUENCE"):
                    temp_line = line.split("SEQUENCE")[1].strip().replace('RS:','').replace('GB:','')
                    temp_line = re.sub('\s*\([^)]*\)', '', temp_line)
                    sequence_ids += temp_line.split(' ')
                elif line.startswith("DATA_SOURCE"):
                    match = re.search('(GCF|GCA)_\d+\.\d+', line)
                    if match:
                        assembly_id = match.group()
                
            return [gn_id, org_code, desc, taxon_id, key_words, lineage, assembly_id, sequence_ids]


        # Get the input file path
        gene_path = os.path.join(self.kegg_data_dir, "genes")
        in_path = os.path.join(gene_path, "genome.tar.gz")
        if not os.path.exists(in_path):
            self.logger.warning(f"File {in_path} does not exist!")
        else:
            # Open the compressed file for reading as binary data
            with open(in_path, 'rb') as f:
                # Decompress the file
                decompressed_data = gzip.decompress(f.read())
                # Wrap the decompressed data in an io.BytesIO object
                fileobj = io.BytesIO(decompressed_data)
                # Open the tar archive from the wrapped data
                with tarfile.open(fileobj=fileobj, mode='r') as tar:
                    text = tar.extractfile('genome/genome').read().decode('utf-8').split("///\n")

            # Convert the output table to a pandas dataframe
            table = pd.DataFrame([extract_info(entry) for entry in text if len(entry) > 0], columns=['gn_id', 'org_code', 'desc', 'taxon_id', 'key_words', 'kegg_lineage', 'assembly_id', 'sequence_ids'])
            # find NCBI lineage and their taxon ids
            result = pytaxonkit.lineage(list(table['taxon_id']))
            result['TaxID'] = result['TaxID'].astype('str')
            result = result[['TaxID','FullLineageTaxIDs','FullLineage','Rank']]
            result.columns = ['taxon_id','ncbi_full_lineage_taxids','ncbi_full_lineage','ncbi_rank']
            table = table.merge(result, on='taxon_id', how='left')

            # remove duplicates
            nonduplicate_indexes = table[['gn_id','org_code']].drop_duplicates().index
            table = table.loc[nonduplicate_indexes,:].reset_index(drop=True)

            ## Return dataframe
            return table

    @staticmethod
    def _download_genome_accession(accession_list: List, output_dir: str, file_name: str):
        """
        Download genome based on accession from NCBI
        """

        # Downloading...
        seqs = []
        for accession in accession_list:
            try:
                handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            except:
                print(f"Error: Fail to donwload nucleotide sequence for accession:{accession}", flush=True)
                continue
            seqs += [handle.read().replace('\n\n','\n')]
            handle.close()

        # save sequences
        with open(os.path.join(output_dir, file_name), 'w') as f:
            for seq in seqs:
                if seq is not None:
                    f.write(seq)


    @staticmethod
    def _download_genome_assembly(assembly_accession: str, output_dir: str, file_name: str):
        """
        Download genome based on assembly accession from NCBI
        """
        # Downloading...
        try:
            # Fetch the assembly summary for the given assembly accession
            esearch_handle = Entrez.esearch(db="assembly", term=assembly_accession, retmax=1)
            esearch_record = Entrez.read(esearch_handle)
            assembly_id = esearch_record["IdList"][0]
            esearch_handle.close()

            # Fetch assembly details
            esummary_handle = Entrez.esummary(db="assembly", id=assembly_id)
            esummary_record = Entrez.read(esummary_handle)
            ftp_path = esummary_record["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
            esummary_handle.close()

            # Download and save the compressed FASTA file
            assembly_name = ftp_path.split("/")[-1]
            fasta_file_name = f"{assembly_name}_genomic.fna.gz"
            fasta_url = f"{ftp_path}/{fasta_file_name}".replace("ftp://", "https://") 
            compressed_output_file = os.path.join(output_dir, fasta_file_name)
            response = requests.get(fasta_url)
            if response.status_code == 200:
                with open(compressed_output_file, "wb") as out_file:
                    out_file.write(response.content)

                # Uncompress the FASTA file
                uncompressed_output_file = os.path.join(output_dir, file_name)
                with gzip.open(compressed_output_file, "rb") as compressed_file:
                    with open(uncompressed_output_file, "wb") as uncompressed_file:
                        uncompressed_file.write(compressed_file.read())

                # Remove the compressed file
                os.remove(compressed_output_file)

                return True
            else:
                print(f"Error: Fail to donwload assembly for {assembly_accession}", flush=True)
                return False

        except:
            print(f"Error: Fail to donwload assembly for {assembly_accession}", flush=True)
            return False


    def extract_organism_gene_seq_info(self, org_code: str):
        """
        Extract gene sequence information from the organism files
        """
        # Get the input file path
        gene_path = os.path.join(self.kegg_data_dir, "genes")
        in_path = glob(os.path.join(gene_path, "organisms", org_code, "*.ent.gz"))[0]
        if not os.path.exists(in_path):
            self.logger.warning(f"File {in_path} does not exist!")
        else:
            # Read data
            with gzip.open(in_path, 'rb') as f:
                text = f.read()
                text = text.decode('utf-8').split("///\n")

            out_res = ['\t'.join(self._extract_info(org_code, entry)) for entry in text]

            # Save data
            out_path = os.path.join(self.output_dir,  "organisms", "kegg_gene_info", f"{org_code}_genes.txt")
            if not os.path.exists(os.path.dirname(out_path)):
                os.makedirs(os.path.dirname(out_path))
            header = "gene_id\tsymbol\tdesc\taaseq\tntseq\n"
            with open(out_path, 'w') as f:
                f.write(header)
                f.write('\n'.join(out_res))

    def extract_organism_gene_link_info(self, org_code: str):
        """
        Extract organism gene link information from the organism-related files
        """

        gene_path = os.path.join(self.kegg_data_dir, "genes")
        in_path = os.path.join(gene_path, "organisms", org_code, f"{org_code}_link.tar.gz")
        if not os.path.exists(in_path):
            self.logger.warning(f"File {in_path} does not exist!")
        else:
            # Open the compressed file for reading as binary data
            with open(in_path, 'rb') as f:
                # Decompress the file
                decompressed_data = gzip.decompress(f.read())
                # Wrap the decompressed data in an io.BytesIO object
                fileobj = io.BytesIO(decompressed_data)
                # Open the tar archive from the wrapped data
                with tarfile.open(fileobj=fileobj, mode='r') as tar:
                    for content in ['brite_contents', 'enzyme_contents', 'ko_contents', 'module_contents', 'ncbi_proteinid_contents', 'pathway_contents', 'pdb_contents', 'pfam_contents', 'uniprot_contents']:
                        try:
                            file_name = f"{org_code}_{content.replace('_contents','').replace('_','-')}.list"
                            temp = tar.extractfile(file_name).read().decode('utf-8').split('\n')
                            exec(f"{content}=temp")
                        except:
                            pass

            # group multiple records to a single record
            brite_dict = {}
            enzyme_dict = {}
            ko_dict = {}
            module_dict = {}
            ncbi_proteinid_dict= {}
            pathway_dict = {}
            pdb_dict = {}
            pfam_dict = {}
            uniprot_dict = {}
            for content in ['brite_contents', 'enzyme_contents', 'ko_contents', 'module_contents', 'ncbi_proteinid_contents', 'pathway_contents', 'pdb_contents', 'pfam_contents', 'uniprot_contents']:
                if content not in locals():
                    continue
                file_contents = eval(content)
                file_dict = eval(content.replace('contents', 'dict'))
                for line in file_contents:
                    if len(line) == 0:
                        continue
                    gene_id, this_content = line.split('\t')
                    if gene_id in file_dict:
                        file_dict[gene_id] += f";{this_content}"
                    else:
                        file_dict[gene_id] = this_content
                exec(f"del {content}")

            # Save dictionary to a text file
            for x in ['brite', 'enzyme', 'ko', 'module', 'ncbi_proteinid', 'pathway', 'pdb', 'pfam', 'uniprot']:
                this_dic = eval(f'{x}_dict')
                if len(this_dic) == 0:
                    continue
                out_path = os.path.join(self.output_dir, "organisms", f"link_{x}_to_gene", f"{org_code}_link_{x}_to_gene.txt")
                if not os.path.exists(os.path.dirname(out_path)):
                    os.makedirs(os.path.dirname(out_path))
                with open(out_path, 'w') as f:
                    for gene_id, content in this_dic.items():
                        f.write(f"{gene_id}\t{content}\n")

    def extract_virus_seq_info(self):
        """
        Extract sequence information from the virus files
        """

        gene_path = os.path.join(self.kegg_data_dir, "genes")
        in_path = glob(os.path.join(gene_path, "viruses", "*.ent.gz"))[0]
        if not os.path.exists(in_path):
            self.logger.warning(f"File {in_path} does not exist!")
        else:
            # Read data
            with gzip.open(in_path, 'rb') as f:
                text = f.read()
                text = text.decode('utf-8').split("///\n")

            out_res = ['\t'.join(self._extract_info('vg', entry)) for entry in text]

            # Save data
            out_path = os.path.join(self.output_dir, "viruses", "kegg_gene_info", f"viruses_genes.txt")
            if not os.path.exists(os.path.dirname(out_path)):
                os.makedirs(os.path.dirname(out_path))
            header = "gene_id\tsymbol\tdesc\taaseq\tntseq\n"
            with open(out_path, 'w') as f:
                f.write(header)
                f.write('\n'.join(out_res))

    def extract_virus_gene_link_info(self):
        """
        Extract virus gene link information from the virus-related files
        """

        gene_path = os.path.join(self.kegg_data_dir, "genes")
        in_path = os.path.join(gene_path, "viruses", "viruses_link.tar.gz")
        if not os.path.exists(in_path):
            self.logger.warning(f"File {in_path} does not exist!")
        else:
            # Open the compressed file for reading as binary data
            with open(in_path, 'rb') as f:
                # Decompress the file
                decompressed_data = gzip.decompress(f.read())
                # Wrap the decompressed data in an io.BytesIO object
                fileobj = io.BytesIO(decompressed_data)
                # Open the tar archive from the wrapped data
                with tarfile.open(fileobj=fileobj, mode='r') as tar:
                    for content in ['enzyme_contents', 'ko_contents', 'ncbi_geneid_contents', 'ncbi_proteinid_contents', 'pdb_contents', 'pfam_contents', 'rs_contents', 'uniprot_contents']:
                        try:
                            file_name = f"vg_{content.replace('_contents','').replace('_','-')}.list"
                            temp = tar.extractfile(file_name).read().decode('utf-8').split('\n')
                            exec(f"{content}=temp")
                        except:
                            pass

            # group multiple records to a single record
            enzyme_dict = {}
            ko_dict = {}
            ncbi_geneid_dict = {}
            ncbi_proteinid_dict = {}
            pdb_dict = {}
            pfam_dict = {}
            rs_dict = {}
            uniprot_dict = {}
            for content in ['enzyme_contents', 'ko_contents', 'ncbi_geneid_contents', 'ncbi_proteinid_contents', 'pdb_contents', 'pfam_contents', 'rs_contents', 'uniprot_contents']:
                if content not in locals():
                    continue
                file_contents = eval(content)
                file_dict = eval(content.replace('contents', 'dict'))
                for line in file_contents:
                    if len(line) == 0:
                        continue
                    gene_id, this_content = line.split('\t')
                    if gene_id in file_dict:
                        file_dict[gene_id] += f";{this_content}"
                    else:
                        file_dict[gene_id] = this_content
                exec(f"del {content}")

            # Save dictionary to a text file
            for x in ['enzyme', 'ko', 'ncbi_geneid', 'ncbi_proteinid', 'pdb', 'pfam', 'rs', 'uniprot']:
                this_dic = eval(f'{x}_dict')
                if len(this_dic) == 0:
                    continue
                out_path = os.path.join(self.output_dir, "viruses", f"vg_link_{x}_to_gene.txt")
                if not os.path.exists(os.path.dirname(out_path)):
                    os.makedirs(os.path.dirname(out_path))
                with open(out_path, 'w') as f:
                    for gene_id, content in this_dic.items():
                        f.write(f"{gene_id}\t{content}\n")

    def extract_link_info_via_api(self):
        """
        Extract link information from the KEGG API
        """

        get_kegg_data = GetKeggLinkData()
        ## download kegg compound information
        get_kegg_data.download_kegg_compound(self.output_dir)
        ## download kegg pathway information
        get_kegg_data.download_kegg_pathway(self.output_dir)
        ## download kegg module information
        get_kegg_data.download_kegg_module(self.output_dir)
        ## download kegg glycan information
        get_kegg_data.download_kegg_glycan(self.output_dir)
        ## download kegg reaction information
        get_kegg_data.download_kegg_reaction(self.output_dir)
        ## download kegg enzyme information
        get_kegg_data.download_kegg_enzyme(self.output_dir)
        ## download kegg network information
        get_kegg_data.download_kegg_network(self.output_dir)
        ## download kegg ko information
        get_kegg_data.download_kegg_ko(self.output_dir)
        ## download kegg disease information
        get_kegg_data.download_kegg_disease(self.output_dir)
        ## download kegg drug information
        get_kegg_data.download_kegg_drug(self.output_dir)
        ## download kegg reaction class information
        get_kegg_data.download_kegg_rclass(self.output_dir)
        ## download kegg drug class information
        get_kegg_data.download_kegg_dgroup(self.output_dir)

        ## link compound information to pathway information
        get_kegg_data.link_compound_to_pathway(self.output_dir)
        ## link compound information to module information
        get_kegg_data.link_compound_to_module(self.output_dir)
        ## link compound information to glycan information
        get_kegg_data.link_compound_to_glycan(self.output_dir)
        ## link compound information to reaction information
        get_kegg_data.link_compound_to_reaction(self.output_dir)
        ## link compound information to enzyme information
        get_kegg_data.link_compound_to_enzyme(self.output_dir)
        ## link compound information to network information
        get_kegg_data.link_compound_to_network(self.output_dir)
        ## link compound information to drug information
        get_kegg_data.link_compound_to_drug(self.output_dir)
        ## link compound information to genome information
        get_kegg_data.link_compound_to_gn(self.output_dir)
        ## link pathway information to module information
        get_kegg_data.link_pathway_to_module(self.output_dir)
        ## link pathway information to glycan information
        get_kegg_data.link_pathway_to_glycan(self.output_dir)
        ## link pathway information to reaction information
        get_kegg_data.link_pathway_to_reaction(self.output_dir)
        ## link pathway information to enzyme information
        get_kegg_data.link_pathway_to_enzyme(self.output_dir)
        ## link pathway information to network information
        get_kegg_data.link_pathway_to_network(self.output_dir)
        ## link pathway information to ko information
        get_kegg_data.link_pathway_to_ko(self.output_dir)
        ## link pathway information to disease information
        get_kegg_data.link_pathway_to_disease(self.output_dir)
        ## link pathway information to drug information
        get_kegg_data.link_pathway_to_drug(self.output_dir)
        ## link pathway information to genome information
        get_kegg_data.link_pathway_to_gn(self.output_dir)
        ## link pathway information to reaction class information
        get_kegg_data.link_pathway_to_rclass(self.output_dir)
        ## link module information to glycan information
        get_kegg_data.link_module_to_glycan(self.output_dir)
        ## link module information to reaction information
        get_kegg_data.link_module_to_reaction(self.output_dir)
        ## link module information to enzyme information
        get_kegg_data.link_module_to_enzyme(self.output_dir)
        ## link module information to ko information
        get_kegg_data.link_module_to_ko(self.output_dir)
        ## link module information to disease information
        get_kegg_data.link_module_to_disease(self.output_dir)
        ## link module information to genome information
        get_kegg_data.link_module_to_gn(self.output_dir)
        ## link glycan information to reaction information
        get_kegg_data.link_glycan_to_reaction(self.output_dir)
        ## link glycan information to reaction information
        get_kegg_data.link_glycan_to_enzyme(self.output_dir)
        ## link reaction information to enzyme information
        get_kegg_data.link_reaction_to_enzyme(self.output_dir)
        ## link reaction information to ko information
        get_kegg_data.link_reaction_to_ko(self.output_dir)
        ## link reaction information to rclass information
        get_kegg_data.link_reaction_to_rclass(self.output_dir)
        ## link enzyme information to ko information
        get_kegg_data.link_enzyme_to_ko(self.output_dir)
        ## link enzyme information to rclass information
        get_kegg_data.link_enzyme_to_rclass(self.output_dir)
        ## link network information to ko information
        get_kegg_data.link_network_to_ko(self.output_dir)
        ## link network information to disease information
        get_kegg_data.link_network_to_disease(self.output_dir)
        ## link ko information to disease information
        get_kegg_data.link_ko_to_disease(self.output_dir)
        ## link ko information to drug information
        get_kegg_data.link_ko_to_drug(self.output_dir)
        ## link ko information to reaction class information
        get_kegg_data.link_ko_to_rclass(self.output_dir)
        ## link disease information to drug information
        get_kegg_data.link_disease_to_drug(self.output_dir)
        ## link disease information to genome information
        get_kegg_data.link_disease_to_gn(self.output_dir)
        ## link drug information to drug group information
        get_kegg_data.link_drug_to_dgroup(self.output_dir)


    def get_ko_hierarchy(self):
        """
        Organize the hierarchy information
        """
        def _iterate_json(res_dict, hiearchy_json, res_list):
            if isinstance(hiearchy_json, dict):
                if len(hiearchy_json.keys()) == 1:
                    name = hiearchy_json['name']
                    if name == 'ko00001':
                        name = 'KEGG Orthology (KO) [BR:ko00001]'
                    koid = name.split(' ',1)[0]
                    if koid not in res_dict:
                        res_dict[koid] = []
                        res_dict[koid].append('|'.join(res_list + [name]))
                    else:
                        res_dict[koid].append('|'.join(res_list + [name]))
                else:
                    name = hiearchy_json['name']
                    if name == 'ko00001':
                        name = 'KEGG Orthology (KO) [BR:ko00001]'
                    for elem in hiearchy_json['children']:
                        _iterate_json(res_dict, elem, res_list + [name])

        link = 'http://rest.kegg.jp/get/br:ko00001/json'
        res = requests.get(link)
        if res.status_code == 200:
            res_dict = {}
            _iterate_json(res_dict, res.json(), [])
            return res_dict
        else:
            self.error(f"Fail to download KEGG hierarchy from {link}")
            return None

    def get_kegg_genes(self, microb_only=True):
        """
        Extract KEGG gene information
        """
        
        if microb_only:
            microbe_labels = ['Archaea', 'Viruses', 'Bacteria', 'Fungi']
            # Create a regular expression pattern from the microbe labels
            pattern = '|'.join(microbe_labels)
            # Extract non-microbial data
            nonmicrobe_gn_table = self.all_gn_table.loc[~self.all_gn_table['kegg_lineage'].str.contains(pattern, regex=True)].reset_index(drop=True)
            # Extract non-microbial organism codes
            removed_organism_list = nonmicrobe_gn_table.loc[nonmicrobe_gn_table['org_code']!='','org_code'].to_list()
        else:
            removed_organism_list = []

        # Get the input file path
        gene_path = os.path.join(self.kegg_data_dir, "genes")
        in_path = os.path.join(gene_path, "ko.tar.gz")
        if not os.path.exists(in_path):
            self.logger.warning(f"File {in_path} does not exist!")
            return None
        else:
            # Open the compressed file for reading as binary data
            with open(in_path, 'rb') as f:
                # Decompress the file
                decompressed_data = gzip.decompress(f.read())
                # Wrap the decompressed data in an io.BytesIO object
                fileobj = io.BytesIO(decompressed_data)
                # Open the tar archive from the wrapped data
                with tarfile.open(fileobj=fileobj, mode='r') as tar:
                    text = tar.extractfile('ko/ko_genes.list').read().decode('utf-8').split("\n")

            # Create a dictionary to store the extracted data
            ko_gene_dict = {}
            for line in tqdm(text):
                if line != '':
                    koid, gene_id = line.split("\t")
                    if koid.split(':')[1] not in ko_gene_dict:
                        if gene_id.split(':')[0] not in removed_organism_list:
                            ko_gene_dict[koid.split(':')[1]] = [gene_id]
                    else:
                        if gene_id.split(':')[0] not in removed_organism_list:
                            ko_gene_dict[koid.split(':')[1]].append(gene_id)
            return ko_gene_dict

    def get_kegg_synonyms(self):
        """
        Extract KEGG synonym ids
        """
        # Get the input file paths
        compound_path = os.path.join(self.kegg_data_dir, "ligand", "compound.tar.gz")
        pathway_path = os.path.join(self.kegg_data_dir, "pathway", "pathway.gz")
        module_path = os.path.join(self.kegg_data_dir, "module", "module.gz")
        glycan_path = os.path.join(self.kegg_data_dir, "ligand", "glycan.tar.gz")
        reaction_path = os.path.join(self.kegg_data_dir, "ligand", "reaction.tar.gz")
        enzyme_path = os.path.join(self.kegg_data_dir, "ligand", "enzyme.tar.gz")
        network_path = os.path.join(self.kegg_data_dir, "medicus", "network.tar.gz")
        ko_path = os.path.join(self.kegg_data_dir, "genes", "ko.tar.gz")
        disease_path = os.path.join(self.kegg_data_dir, "medicus", "disease.tar.gz")
        drug_path = os.path.join(self.kegg_data_dir, "medicus", "drug.tar.gz")
        rclass_path = os.path.join(self.kegg_data_dir, "ligand", "rclass.tar.gz")
        dgroup_path = os.path.join(self.kegg_data_dir, "medicus", "dgroup.tar.gz")
        gz_files = [pathway_path, module_path]
        tar_gz_files = [compound_path, glycan_path, reaction_path, enzyme_path, network_path, ko_path, disease_path, drug_path, rclass_path, dgroup_path]
        mapping = {'pathway': "path", 'module': "md", "compound": "cpd", "glycan": "gl", "reaction": 'rn', "enzyme": "ec", "network": "ne", "ko": "ko", "disease": "ds", "drug": "dr", "rclass": "rc", "dgroup": "dg"}

        # Extract dblink information
        dblink_list = []
        for gz_file in tqdm(gz_files, desc=f"Extracting dblink information"):
            if not os.path.exists(gz_file):
                self.logger.warning(f"File {gz_file} does not exist!")
            else:
                # get prefix
                temp_prefix = os.path.basename(gz_file).replace(".gz", "")
                prefix = mapping[temp_prefix]
                with gzip.open(gz_file, 'rt') as f:
                    entries = f.read().split("///\n")

                dblink_list += [self.extract_dblinks(entry, prefix) for entry in entries if len(entry) !=0]

        for tar_gz_file in tqdm(tar_gz_files, desc=f"Extracting dblink information"):
            if not os.path.exists(tar_gz_file):
                self.logger.warning(f"File {tar_gz_file} does not exist!")
            else:
                # Open the compressed file for reading as binary data
                with open(tar_gz_file, 'rb') as f:
                    # get prefix
                    temp_prefix = os.path.basename(tar_gz_file).replace(".tar.gz", "")
                    prefix = mapping[temp_prefix]
                    if prefix == 'cpd':
                        filename = 'compound/compound'
                    elif prefix == 'gl':
                        filename = 'glycan/glycan'
                    elif prefix == 'rn':
                        filename ='reaction/reaction'
                    elif prefix == 'ec':
                        filename = 'enzyme/enzyme'
                    elif prefix == 'ne':
                        filename = 'network/network'
                    elif prefix == 'ko':
                        filename = 'ko/ko'
                    elif prefix == 'ds':
                        filename = 'disease/disease'
                    elif prefix == 'dr':
                        filename = 'drug/drug'
                    elif prefix == 'rc':
                        filename = 'rclass/rclass'
                    elif prefix == 'dg':
                        filename = 'dgroup/dgroup'
                    # Decompress the file
                    decompressed_data = gzip.decompress(f.read())
                    # Wrap the decompressed data in an io.BytesIO object
                    fileobj = io.BytesIO(decompressed_data)
                    # Open the tar archive from the wrapped data
                    with tarfile.open(fileobj=fileobj, mode='r') as tar:
                        try:
                            entries = tar.extractfile(filename).read().decode('utf-8').split("///\n")
                        except:
                            entries = tar.extractfile(filename).read().decode('windows-1252').split("///\n")

                    dblink_list += [self.extract_dblinks(entry, prefix) for entry in entries if len(entry) !=0]

        # convert dblink list to a dictionary
        dblink_dict = {f"KEGG:{kegg_id}":dblink_list for kegg_id, dblink_list in tqdm(dblink_list)}

        return dblink_dict

    def download_gn_seq(self, gn_ids: Union[List, str]):
        """
        Download KEGG genome sequences
        """
        if type(gn_ids) is str:
            gn_ids = [gn_ids]
        
        # Check if the given gn_ids are valid
        diff_gn_ids = set(gn_ids).difference(set(self.all_gn_table['gn_id']))
        if len(diff_gn_ids) > 0:
            self.logger.warning(f"Invalid gn_ids: {diff_gn_ids}")

        valid_gn_ids = list(set(gn_ids).intersection(set(self.all_gn_table['gn_id'])))
        if len(valid_gn_ids) > 0:
            # Download genome sequences
            self.logger.info(f"Downloading genome sequences for {len(valid_gn_ids)} genomes")
            for gn_id in valid_gn_ids:
                self.logger.info(f"Downloading genome sequence for {gn_id}")

                if not os.path.exists(os.path.join(self.output_dir, "fasta")):
                    os.makedirs(os.path.join(self.output_dir, "fasta"))

                # First check if the genome has assembly accession
                assembly_accession = self.all_gn_table.loc[self.all_gn_table['gn_id']==gn_id,'assembly_id'].to_numpy()[0]
                if len(assembly_accession) != 0:
                    
                    res = self._download_genome_assembly(assembly_accession, os.path.join(self.output_dir, "fasta"), f"{gn_id}.fna")
                    if not res:
                        # get RS_id or GB_id
                        seq_ids = self.all_gn_table.loc[self.all_gn_table['gn_id']==gn_id,'sequence_ids'].to_numpy()[0]
                        if len(seq_ids) == 0:
                            self.logger.warning(f"No sequence id for {gn_id}")
                            continue
                        self._download_genome_accession(seq_ids, os.path.join(self.output_dir, "fasta"), f"{gn_id}.fna")

                else:
                    self.logger.warning(f"No assembly accession for {gn_id}, try to use genome accession instead")
                    # get RS_id or GB_id
                    seq_ids = self.all_gn_table.loc[self.all_gn_table['gn_id']==gn_id,'sequence_ids'].to_numpy()[0]
                    if len(seq_ids) == 0:
                        self.logger.warning(f"No sequence id for {gn_id}")
                        continue
                    self._download_genome_accession(seq_ids, os.path.join(self.output_dir, "fasta"), f"{gn_id}.fna")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract KEGG data')
    parser.add_argument('--kegg_data_dir', type=str, help='path of the KEGG FTP data directory')
    parser.add_argument('--microb_only', action='store_true', help="only extract microbial data (e.g. 'Archaea', 'viruses', 'Bacteria', 'Fungi')", default=False)
    parser.add_argument('--output_dir', type=str, help='path of the output directory')
    args = parser.parse_args()

    # Create a logger object
    logger = get_logger()
    logger.setLevel(logging.INFO)

    # Create KEGGData object
    keggdata = KEGGData(args.kegg_data_dir, args.output_dir)
    
    # Create output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    if args.microb_only:
        logger.info("Extracting microbial data (e.g. 'Archaea', 'viruses', 'Bacteria', 'Fungi') only")
        microbe_labels = ['Archaea', 'Viruses', 'Bacteria', 'Fungi']
        # Create a regular expression pattern from the microbe labels
        pattern = '|'.join(microbe_labels)
        # Extract microbial data
        microbe_gn_table = keggdata.all_gn_table.loc[keggdata.all_gn_table['kegg_lineage'].str.contains(pattern, regex=True)].reset_index(drop=True)
        # Extract microbial organism codes
        organism_list = list(set(microbe_gn_table.loc[microbe_gn_table['org_code']!='','org_code'].to_list()))
        logger.info(f"Starting to extract {len(organism_list)} microbial organisms")
        for org_code in tqdm(organism_list):
            keggdata.extract_organism_gene_seq_info(org_code)
            keggdata.extract_organism_gene_link_info(org_code)
        logger.info(f"Finished extracting {len(organism_list)} microbial organisms")
    else:
        organism_list = keggdata.all_gn_table['org_code'].to_list()
        logger.info(f"Starting to extract {len(organism_list)} microbial organisms")
        for org_code in tqdm(organism_list):
            keggdata.extract_organism_gene_seq_info(org_code)
            keggdata.extract_organism_gene_link_info(org_code)
        logger.info(f"Finished extracting {len(organism_list)} microbial organisms")

    logger.info("Starting to extract viruses data")
    keggdata.extract_virus_seq_info()
    keggdata.extract_virus_gene_link_info()
    logger.info("Finished extracting viruses data")
    logger.info("Satrting to extract KEGG link info via APIs")
    keggdata.extract_link_info_via_api()
    logger.info("Finished extracting KEGG link info via APIs")

    logger.info("Starting to download genome sequences")
    archaea_bacteria_labels = ['Archaea', 'Bacteria']
    pattern = '|'.join(archaea_bacteria_labels)
    archaea_bacteria_gn_table = keggdata.all_gn_table.loc[keggdata.all_gn_table['kegg_lineage'].str.contains(pattern, regex=True)].reset_index(drop=True)
    keggdata.download_gn_seq(list(archaea_bacteria_gn_table['gn_id']))
    logger.info("Finished downloading genome sequences")


