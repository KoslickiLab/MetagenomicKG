"""
This script is used to integrate all KEGG data into a knolwedge graph.
"""

## Import standard libraries
import os
import sys
from tqdm import tqdm, trange
from glob import glob
import pandas as pd
import argparse
import re
import logging

## Import custom libraries
from utils import get_logger, read_tsv_file, Node, Edge, KnowledgeGraph
from extract_KEGG_data import KEGGData


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Integrate all KEGG data into a knolwedge graph based on the above graph structure')
    parser.add_argument('--existing_KG_nodes', type=str, help='path of the existing knowledge graph nodes')
    parser.add_argument('--existing_KG_edges', type=str, help='path of the existing knowledge graph edges')
    parser.add_argument('--kegg_data_dir', type=str, help='path of the KEGG FTP data directory')
    parser.add_argument('--kegg_processed_data_dir', type=str, help='path of the KEGG processed data directory')
    parser.add_argument('--gtdb_assignment', type=str, help='path of the gtdb assignment using GTDB-Tk')
    parser.add_argument('--microb_only', action='store_true', help="only extract microbial data (e.g. 'Archaea', 'viruses', 'Bacteria', 'Fungi')", default=False)
    parser.add_argument('--ANI_threshold', type=float, help='ANI threshold to identify the same strain (default 0 for no filtering)', default=0.0)
    parser.add_argument('--AF_threshold', type=float, help='AF threshold to dentify the same strain (default 0 for no filtering)', default=0.0)
    parser.add_argument('--output_dir', type=str, help='path of the output directory')
    args = parser.parse_args()

    # Create a logger object
    logger = get_logger()
    logger.setLevel(logging.INFO)

    # Create KEGGData object
    keggdata = KEGGData(args.kegg_data_dir, args.output_dir)

    ## Get KO hierarchy
    logger.info('Getting KO hierarchy')
    ko_hierarchy = keggdata.get_ko_hierarchy()

    ## Get KO related genes
    logger.info('Getting KO related genes')
    ko_related_genes = keggdata.get_kegg_genes(args.microb_only)

    ## Get KEGG synonyms
    logger.info('Getting KEGG synonyms')
    kegg_synonyms = keggdata.get_kegg_synonyms()
    # merge kegg_synonyms
    temp_kegg_synonyms = {}
    for key in tqdm(kegg_synonyms, desc='Merging KEGG synonyms'):
        temp = key.split(':')
        new_key = temp[0]+':'+temp[1]+'_'+temp[2]
        if temp[1] == 'path':
            temp[2] = re.findall(r'\d+', temp[2])[0]
            new_key = temp[0]+':'+temp[1]+'_map'+temp[2]
        if new_key not in temp_kegg_synonyms:
            temp_kegg_synonyms[new_key] = kegg_synonyms[key]
        else:
            temp_kegg_synonyms[new_key] = list(set(temp_kegg_synonyms[new_key] + kegg_synonyms[key]))
    kegg_synonyms = temp_kegg_synonyms
    
    # remove duplicate synonyms
    synonym_count = {}
    for key in tqdm(kegg_synonyms, desc='Finding duplicate synonyms'):
        for synonym in kegg_synonyms[key]:
            if synonym not in synonym_count:
                synonym_count[synonym] = 1
            else:
                synonym_count[synonym] += 1
    duplicate_synonym = [synonym for synonym in synonym_count if synonym_count[synonym] != 1]
    for key in tqdm(kegg_synonyms, desc='Removing duplicate synonyms'):
        kegg_synonyms[key] = list(set(kegg_synonyms[key]).difference(set(duplicate_synonym)))

    # Create a knowledge graph object
    logger.info("Creating a knowledge graph object...")
    kg = KnowledgeGraph(logger)
    # Load existing knowledge graph nodes and edges
    logger.info("Loading existing knowledge graph nodes and edges...")
    load_dir = args.output_dir
    node_filename = args.existing_KG_nodes.split('/')[-1]
    edge_filename = args.existing_KG_edges.split('/')[-1]
    kg.load_graph(load_dir=load_dir, node_filename=node_filename, edge_filename=edge_filename)

    ## read KEGG archeaa and bacteria assignment
    logger.info('Reading KEGG archeaa and bacteria assignment')
    kegg_archaea_bacteria_assignment = read_tsv_file(args.gtdb_assignment)
    temp_header = kegg_archaea_bacteria_assignment[0]
    mapping_gn_to_microbe_id = {}
    for row in tqdm(kegg_archaea_bacteria_assignment[1:]):
        if (row[temp_header.index('fastani_ani')] != 'N/A' and float(row[temp_header.index('fastani_ani')]) >= args.ANI_threshold) and (row[temp_header.index('fastani_af')] != 'N/A' and float(row[temp_header.index('fastani_af')]) >= args.AF_threshold):
            mapping_gn_to_microbe_id[f"KEGG:gn_{row[temp_header.index('user_genome')]}"] = [f"GTDB:{row[temp_header.index('fastani_reference')]}"]
        else:
            mapping_gn_to_microbe_id[f"KEGG:gn_{row[temp_header.index('user_genome')]}"] = [f"KEGG:gn_{row[temp_header.index('user_genome')]}"]
        mapping_gn_to_microbe_id[f"KEGG:gn_{row[temp_header.index('user_genome')]}"] += [row[temp_header.index('classification')], f"ANI_reference_radius:{row[temp_header.index('fastani_reference_radius')]}; ANI:{row[temp_header.index('fastani_ani')]}; AF:{row[temp_header.index('fastani_af')]}; classification_method:{row[temp_header.index('classification_method')]}; MSA_percent:{row[temp_header.index('msa_percent')]}"]

    if args.microb_only:
        logger.info("Only microbial organism data (e.g. 'Archaea', 'viruses', 'Bacteria', 'Fungi') are used in KG construction.")
        microbe_labels = ['Archaea', 'Viruses', 'Bacteria', 'Fungi']
        # Create a regular expression pattern from the microbe labels
        pattern = '|'.join(microbe_labels)
        # Extract microbial data
        microbe_gn_table = keggdata.all_gn_table.loc[keggdata.all_gn_table['kegg_lineage'].str.contains(pattern, regex=True)].reset_index(drop=True)
    else:
        microbe_gn_table = keggdata.all_gn_table

    # Create output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    ## Get all KEGG-baesd data
    logger.info('Get all KEGG-baesd data')
    parent_dict = {}
    orgcode_to_gnid = {}
    # KEGG genomes
    abb_mapping = {'d': 'superkingdom', 'p': 'phylum', 'c': 'class', 'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species'}
    next_rank = {'superkingdom':'phylum', 'phylum':'class', 'class':'order', 'order':'family', 'family':'genus', 'genus':'species', 'species':'strain'}
    for row in tqdm(microbe_gn_table.to_numpy(), desc='integrating info of KEGG genomes'):
        gn_id, org_code, desc, taxon_id, keywords, kegg_lineage, assembly_id, sequence_ids, ncbi_full_lineage_taxids, ncbi_full_lineage, ncbi_rank  = row
        
        if org_code != '':
            orgcode_to_gnid[org_code] = f"KEGG:gn_{gn_id}"
            
        if assembly_id != '':
            gtdb_assigned_id, gtdb_classification, gtdb_assignment_info = f'GTDB:{assembly_id}', None, None
        elif f"KEGG:gn_{gn_id}" in mapping_gn_to_microbe_id:
            gtdb_assigned_id, gtdb_classification, gtdb_assignment_info = mapping_gn_to_microbe_id[f"KEGG:gn_{gn_id}"]
        else:
            gtdb_assigned_id, gtdb_classification, gtdb_assignment_info = None, None, None 

        if gtdb_assigned_id:
            node_id = kg.find_node_by_synonym(gtdb_assigned_id)
            if node_id:
                ## kg has this node
                existing_node = kg.get_node_by_id(node_id)
                temp_all_names = [desc, ncbi_full_lineage.split(';')[-1]]
                existing_node.all_names = list(set(existing_node.all_names + temp_all_names))
                description_dict = dict(existing_node.description)
                if 'taxid' in description_dict and description_dict['taxid'] == '':
                    description_dict['taxid'] = taxon_id
                    description_dict['rank'] = ncbi_rank
                elif 'taxid' not in description_dict:
                    description_dict['taxid'] = taxon_id
                    description_dict['rank'] = ncbi_rank
                if gtdb_assignment_info and gtdb_assignment_info.split(';')[1].split(':')[1] != 'N/A':
                    temp_description_dict = dict([tuple(item.strip().split(':')) for item in gtdb_assignment_info.split(';')])
                    description_dict.update(temp_description_dict)
                existing_node.description = list(description_dict.items())
                existing_node.knowledge_source = list(set(existing_node.knowledge_source + ['KEGG']))
                temp_synonyms = [f"KEGG:gn_{gn_id}"]
                for synonym in temp_synonyms:
                    kg.map_synonym_to_node_id[synonym] = node_id
                existing_node.synonyms = list(set(existing_node.synonyms + temp_synonyms))
                existing_node.link = list(set(existing_node.link + [f"https://www.genome.jp/entry/gn:{gn_id}"]))
                if assembly_id != '':
                    existing_node.link = list(set(existing_node.link + [f"https://www.ncbi.nlm.nih.gov/assembly/{assembly_id}"]))
                existing_node.link = list(set(existing_node.link + [f"https://www.ncbi.nlm.nih.gov/nuccore/{sequence_id}" for sequence_id in sequence_ids]))
                existing_node.is_pathogen = True if keywords == 'Human pathogen' else existing_node.is_pathogen
            else:
                if f"KEGG:gn_{gn_id}" in mapping_gn_to_microbe_id:
                    gtdb_assigned_id, gtdb_classification, gtdb_assignment_info = mapping_gn_to_microbe_id[f"KEGG:gn_{gn_id}"]
                
                ## kg does not have this node
                temp_all_names = [desc, ncbi_full_lineage.split(';')[-1]]
                temp_description_dict = {'taxid': taxon_id, 'rank': ncbi_rank}
                temp_knowledge_source = ['KEGG']
                temp_synonyms = [f"KEGG:gn_{gn_id}"]
                temp_link = [f"https://www.genome.jp/entry/gn:{gn_id}"]
                if assembly_id != '':
                    temp_link += [f"https://www.ncbi.nlm.nih.gov/assembly/{assembly_id}"]
                temp_link += [f"https://www.ncbi.nlm.nih.gov/nuccore/{sequence_id}" for sequence_id in sequence_ids]
                temp_is_pathogen = True if keywords == 'Human pathogen' else False
                temp_node = Node(node_type="Microbe", all_names=temp_all_names, description=list(temp_description_dict.items()), knowledge_source=temp_knowledge_source, synonyms=temp_synonyms, link=temp_link, is_pathogen=temp_is_pathogen)
                kg.add_node(temp_node)
                if gtdb_classification:                    
                    parent_dict[f"KEGG:gn_{gn_id}"] = f"GTDB:{[x for x in gtdb_classification.split(';')[::-1] if x.split('__')[1] != ''][0].split('__')[1]}"
        else:
            ## This genome might be a vrirus or a fungi or without a assembly_id
            ncbi_synonym = f"NCBI:{ncbi_full_lineage.split(';')[-1]}"
            node_id = kg.find_node_by_synonym(ncbi_synonym)
            if node_id:
                ## kg has this node
                existing_node = kg.get_node_by_id(node_id)
                temp_all_names = [desc]
                existing_node.all_names = list(set(existing_node.all_names + temp_all_names))
                description_dict = dict(existing_node.description)
                if 'taxid' in description_dict and description_dict['taxid'] == '':
                    description_dict['taxid'] = taxon_id
                    description_dict['rank'] = ncbi_rank
                elif 'taxid' not in description_dict:
                    description_dict['taxid'] = taxon_id
                    description_dict['rank'] = ncbi_rank
                existing_node.description = list(description_dict.items())
                existing_node.knowledge_source = list(set(existing_node.knowledge_source + ['KEGG']))
                temp_synonyms = [f"KEGG:gn_{gn_id}", ncbi_synonym]
                for synonym in temp_synonyms:
                    kg.map_synonym_to_node_id[synonym] = node_id
                existing_node.synonyms = list(set(existing_node.synonyms + temp_synonyms))
                existing_node.link = list(set(existing_node.link + [f"https://www.genome.jp/entry/gn:{gn_id}"]))
                if assembly_id != '':
                    existing_node.link = list(set(existing_node.link + [f"https://www.ncbi.nlm.nih.gov/assembly/{assembly_id}"]))
                existing_node.link = list(set(existing_node.link + [f"https://www.ncbi.nlm.nih.gov/nuccore/{sequence_id}" for sequence_id in sequence_ids]))
                existing_node.is_pathogen = True if keywords == 'Human pathogen' else existing_node.is_pathogen
            else:
                ## kg does not have this node
                temp_all_names = [desc]
                temp_description_dict = {'taxid': taxon_id, 'rank': ncbi_rank}
                temp_knowledge_source = ['KEGG']
                temp_synonyms = [f"KEGG:gn_{gn_id}", ncbi_synonym]
                temp_link = [f"https://www.genome.jp/entry/gn:{gn_id}"]
                if assembly_id != '':
                    temp_link += [f"https://www.ncbi.nlm.nih.gov/assembly/{assembly_id}"]
                temp_link += [f"https://www.ncbi.nlm.nih.gov/nuccore/{sequence_id}" for sequence_id in sequence_ids]
                temp_is_pathogen = True if keywords == 'Human pathogen' else False
                temp_node = Node(node_type="Microbe", all_names=temp_all_names, description=list(temp_description_dict.items()), knowledge_source=temp_knowledge_source, synonyms=temp_synonyms, link=temp_link, is_pathogen=temp_is_pathogen)
                kg.add_node(temp_node)
                parent_dict[f"KEGG:gn_{gn_id}"] = f"NCBI:{ncbi_full_lineage.split(';')[-1]}"

    # KEGG compounds
    file_path = os.path.join(args.kegg_processed_data_dir, 'kegg_compounds.txt')
    compound_table = pd.read_csv(file_path, sep='\t', header=0)
    for row in tqdm(compound_table.to_numpy(), desc='integrating info of KEGG compounds'):
        compound_id, desc = row
        if f"KEGG:cpd_{compound_id}" in kegg_synonyms:
            temp_synonyms = kegg_synonyms[f"KEGG:cpd_{compound_id}"]
            
        temp_node = Node(node_type="Compound", all_names=[desc], description=[], knowledge_source=['KEGG'], synonyms=[f"KEGG:cpd_{compound_id}"] + temp_synonyms, link=[f"https://www.genome.jp/entry/cpd:{compound_id}"])
        kg.add_node(temp_node)
        
    # KEGG pathways
    file_path = os.path.join(args.kegg_processed_data_dir, 'kegg_pathways.txt')
    pathway_table = pd.read_csv(file_path, sep='\t', header=0)
    for row in tqdm(pathway_table.to_numpy(), desc='integrating info of KEGG pathways'):
        pathway_id, desc = row
        if f"KEGG:path_{pathway_id}" in kegg_synonyms:
            temp_synonyms = kegg_synonyms[f"KEGG:path_{pathway_id}"]
            
        temp_node = Node(node_type="Pathway", all_names=[desc], description=[], knowledge_source=['KEGG'], synonyms=[f"KEGG:path_{pathway_id}"] + temp_synonyms, link=[f"https://www.genome.jp/entry/path:{pathway_id}"])
        kg.add_node(temp_node)

    # KEGG modules
    file_path = os.path.join(args.kegg_processed_data_dir, 'kegg_modules.txt')
    module_table = pd.read_csv(file_path, sep='\t', header=0)
    for row in tqdm(module_table.to_numpy(), desc='integrating info of KEGG modules'):
        module_id, desc = row
        if f"KEGG:md_{module_id}" in kegg_synonyms:
            temp_synonyms = kegg_synonyms[f"KEGG:md_{module_id}"]
            
        temp_node = Node(node_type="Module", all_names=[desc], description=[], knowledge_source=['KEGG'], synonyms=[f"KEGG:md_{module_id}"] + temp_synonyms, link=[f"https://www.genome.jp/entry/md:{module_id}"])
        kg.add_node(temp_node)

    # KEGG KOs
    file_path = os.path.join(args.kegg_processed_data_dir, 'kegg_koids.txt')
    kos_table = pd.read_csv(file_path, sep='\t', header=0)
    for row in tqdm(kos_table.to_numpy(), desc='integrating info of KEGG KOs'):
        kos_id, desc = row
        temp_description = []
        if kos_id in ko_hierarchy:
            temp_description += [('KO_hierarchy', '#####'.join(ko_hierarchy[kos_id]))]
        if kos_id in ko_related_genes:
            temp_description += [('KO_related_genes', '#####'.join(ko_related_genes[kos_id]))]
        if f"KEGG:ko_{kos_id}" in kegg_synonyms:
            temp_synonyms = kegg_synonyms[f"KEGG:ko_{kos_id}"]
        
        temp_node = Node(node_type="KO", all_names=[desc], description=temp_description, knowledge_source=['KEGG'], synonyms=[f"KEGG:ko_{kos_id}"] + temp_synonyms, link=[f"https://www.genome.jp/entry/ko:{kos_id}"])
        kg.add_node(temp_node)

    # KEGG diseases
    file_path = os.path.join(args.kegg_processed_data_dir, 'kegg_diseases.txt')
    disease_table = pd.read_csv(file_path, sep='\t', header=0)
    for row in tqdm(disease_table.to_numpy(), desc='integrating info of KEGG diseases'):
        disease_id, desc = row
        if f"KEGG:ds_{disease_id}" in kegg_synonyms:
            temp_synonyms = kegg_synonyms[f"KEGG:ds_{disease_id}"]
        
        temp_node = Node(node_type="Disease", all_names=[desc], description=[], knowledge_source=['KEGG'], synonyms=[f"KEGG:ds_{disease_id}"] + temp_synonyms, link=[f"https://www.genome.jp/entry/ds:{disease_id}"])
        kg.add_node(temp_node)

    # KEGG drugs
    file_path = os.path.join(args.kegg_processed_data_dir, 'kegg_drugs.txt')
    drug_table = pd.read_csv(file_path, sep='\t', header=0)
    for row in tqdm(drug_table.to_numpy(), desc='integrating info of KEGG drugs'):
        drug_id, desc = row
        if f"KEGG:dr_{drug_id}" in kegg_synonyms:
            temp_synonyms = kegg_synonyms[f"KEGG:dr_{drug_id}"]
        
        temp_node = Node(node_type="Drug", all_names=[desc], description=[], knowledge_source=['KEGG'], synonyms=[f"KEGG:dr_{drug_id}"] + temp_synonyms, link=[f"https://www.genome.jp/entry/dr:{drug_id}"])
        kg.add_node(temp_node)

    # KEGG reactions
    file_path = os.path.join(args.kegg_processed_data_dir, 'kegg_reactions.txt')
    reaction_table = pd.read_csv(file_path, sep='\t', header=0)
    for row in tqdm(reaction_table.to_numpy(), desc='integrating info of KEGG reactions'):
        reaction_id, desc = row
        if f"KEGG:rn_{reaction_id}" in kegg_synonyms:
            temp_synonyms = kegg_synonyms[f"KEGG:rn_{reaction_id}"]
        
        temp_node = Node(node_type="Reaction", all_names=[desc], description=[], knowledge_source=['KEGG'], synonyms=[f"KEGG:rn_{reaction_id}"] + temp_synonyms, link=[f"https://www.genome.jp/entry/rn:{reaction_id}"])
        kg.add_node(temp_node)

    # KEGG enzymes
    file_path = os.path.join(args.kegg_processed_data_dir, 'kegg_enzymes.txt')
    enzyme_table = pd.read_csv(file_path, sep='\t', header=0)
    for row in tqdm(enzyme_table.to_numpy(), desc='integrating info of KEGG enzymes'):
        enzyme_id, desc = row
        if f"KEGG:ec_{enzyme_id}" in kegg_synonyms:
            temp_synonyms = kegg_synonyms[f"KEGG:ec_{enzyme_id}"]
        
        temp_node = Node(node_type="Enzyme", all_names=[desc], description=[], knowledge_source=['KEGG'], synonyms=[f"KEGG:ec_{enzyme_id}"] + temp_synonyms, link=[f"https://www.genome.jp/entry/ec:{enzyme_id}"])
        kg.add_node(temp_node)

    # KEGG glycans
    file_path = os.path.join(args.kegg_processed_data_dir, 'kegg_glycans.txt')
    glycan_table = pd.read_csv(file_path, sep='\t', header=0)
    for row in tqdm(glycan_table.to_numpy(), desc='integrating info of KEGG glycans'):
        glycan_id, desc = row
        if f"KEGG:gl_{glycan_id}" in kegg_synonyms:
            temp_synonyms = kegg_synonyms[f"KEGG:gl_{glycan_id}"]
        
        temp_node = Node(node_type="Glycan", all_names=[desc], description=[], knowledge_source=['KEGG'], synonyms=[f"KEGG:gl_{glycan_id}"] + temp_synonyms, link=[f"https://www.genome.jp/entry/gl:{glycan_id}"])
        kg.add_node(temp_node)

    # KEGG networks
    file_path = os.path.join(args.kegg_processed_data_dir, 'kegg_networks.txt')
    network_table = pd.read_csv(file_path, sep='\t', header=0)
    for row in tqdm(network_table.to_numpy(), desc='integrating info of KEGG networks'):
        network_id, desc = row
        # if f"KEGG:ne_{network_id}" in kegg_synonyms:
        #     temp_synonyms = kegg_synonyms[f"KEGG:ne_{network_id}"]
        
        temp_node = Node(node_type="Network", all_names=[desc], description=[], knowledge_source=['KEGG'], synonyms=[f"KEGG:ne_{network_id}"] + temp_synonyms, link=[f"https://www.genome.jp/entry/ne:{network_id}"])
        kg.add_node(temp_node)

    # KEGG drug groups
    file_path = os.path.join(args.kegg_processed_data_dir, 'kegg_dgroups.txt')
    dgroup_table = pd.read_csv(file_path, sep='\t', header=0)
    for row in tqdm(dgroup_table.to_numpy(), desc='integrating info of KEGG drug groups'):
        dgroup_id, desc = row
        if f"KEGG:dg_{dgroup_id}" in kegg_synonyms:
            temp_synonyms = kegg_synonyms[f"KEGG:dg_{dgroup_id}"]
        
        temp_node = Node(node_type="Drug_Group", all_names=[desc], description=[], knowledge_source=['KEGG'], synonyms=[f"KEGG:dg_{dgroup_id}"] + temp_synonyms, link=[f"https://www.genome.jp/entry/dg:{dgroup_id}"])
        kg.add_node(temp_node)

    # KEGG reaction classes
    file_path = os.path.join(args.kegg_processed_data_dir, 'kegg_rclasses.txt')
    rclass_table = pd.read_csv(file_path, sep='\t', header=0)
    for row in tqdm(rclass_table.to_numpy(), desc='integrating info of KEGG reaction classes'):
        rclass_id, desc = row
        if f"KEGG:rc_{rclass_id}" in kegg_synonyms:
            temp_synonyms = kegg_synonyms[f"KEGG:rc_{rclass_id}"]
        
        temp_node = Node(node_type="Reaction", all_names=[desc], description=[], knowledge_source=['KEGG'], synonyms=[f"KEGG:rc_{rclass_id}"] + temp_synonyms, link=[f"https://www.genome.jp/entry/rc:{rclass_id}"])
        kg.add_node(temp_node)

    ## Connect genome hierarchy
    logger.info('Connecting genome hierarchy')
    for genome_child in parent_dict:
        genome_parent = parent_dict[genome_child]
        if genome_parent.split(':')[0] == 'NCBI':
            kg_source = 'NCBI'
        elif genome_parent.split(':')[0] == 'GTDB':
            kg_source = 'GTDB'
        
        # Add edge to the knowledge graph
        kg.add_edge(Edge(source_node=genome_parent, target_node=genome_child, predicate='biolink:has_part', knowledge_source=[kg_source]))
        kg.add_edge(Edge(source_node=genome_child, target_node=genome_parent, predicate='biolink:part_of', knowledge_source=[kg_source]))

    ## Get all KEGG-baesd connections
    logger.info('Getting KEGG-baesd connections')
    link_file_list = glob(os.path.join(args.kegg_processed_data_dir,'link_*'))
    for file_path in tqdm(link_file_list, desc='integrating KEGG connections'):
        logger.info(f'Read {file_path}')
        infile = pd.read_csv(file_path, sep='\t', header=0)
        for row in infile.to_numpy():
            source_node, target_node, source_to_target, target_to_source = row
            if os.path.basename(file_path) in ['link_compound_to_gn.txt', 'link_module_to_gn.txt', 'link_disease_to_gn.txt','link_pathway_to_gn.txt']:
                if not source_node.split(':')[1].startswith('T') and source_node.split(':')[1] not in orgcode_to_gnid:
                    continue
                if source_node.split(':')[1].startswith('T'):
                    source_node = 'KEGG:' + source_node.replace(':','_')
                else:
                    source_node = orgcode_to_gnid[source_node.split(':')[1]]
                target_node = 'KEGG:' + target_node.replace(':','_')
            else:
                source_node = 'KEGG:' + source_node.replace(':','_')
                target_node = 'KEGG:' + target_node.replace(':','_')
            
            ## All genomes that connect to disease nodes are pathogens
            if os.path.basename(file_path) == 'link_disease_to_gn.txt':
                existing_node = kg.get_node_by_id(source_node)
                if existing_node:
                    existing_node.is_pathogen = True
            
            # Add edge to the knowledge graph
            if not isinstance(source_to_target, float):
                kg.add_edge(Edge(source_node=source_node, target_node=target_node, predicate=source_to_target, knowledge_source=['KEGG']))
            if not isinstance(target_to_source, float):
                kg.add_edge(Edge(source_node=target_node, target_node=source_node, predicate=target_to_source, knowledge_source=['KEGG']))


    # connect genome to other node types based on genes
    logger.info('connect genome to other node types based on genes')
    # viruses
    vg_id_to_gnid = {row[0]:f"KEGG:{row[1].replace(':','_')}" for row in keggdata.virus_table[['gene_id','gn_id']].to_numpy()}
    link_file_list = glob(os.path.join(args.kegg_processed_data_dir,'viruses','vg_link_*'))
    for file_path in tqdm(link_file_list, desc='integrating viruses gene-based connections'):
        if os.path.basename(file_path).replace('vg_link_','').replace('_to_gene.txt','') in ['ncbi_geneid','uniprot','pfam','rs','pdb','ncbi_proteinid']:
            continue
        infile = pd.read_csv(file_path, sep='\t', header=None)
        for row in infile.to_numpy():
            virus_gene_id, node_ids = row

            if virus_gene_id in vg_id_to_gnid:
                for node_id in node_ids.split(';'):
                    source_node = vg_id_to_gnid[virus_gene_id]
                    target_node = 'KEGG:' + node_id.replace(':','_')
                    # Add edge to the knowledge graph
                    kg.add_edge(Edge(source_node=source_node, target_node=target_node, predicate='biolink:genetically_associated_with', knowledge_source=['KEGG']))
                    kg.add_edge(Edge(source_node=target_node, target_node=source_node, predicate='biolink:genetically_associated_with', knowledge_source=['KEGG']))

    #FIXME: Consider if we need to connect non-virus genome to other node types based on genes because they are too many
    # # bacteria/fungi/archaea
    # flag = False
    # link_file_list = glob(os.path.join(args.kegg_processed_data_dir,'organisms','link_*'))
    # for dir_path in tqdm(link_file_list, desc='integrating bacteria/fungi/archaea gene-based connections'):
    #     if os.path.basename(dir_path).replace('link_','').replace('_to_gene','') in ['brite','uniprot','pfam','rs','pdb','ncbi_proteinid']:
    #         continue
    #     file_type = os.path.basename(dir_path).replace('link_','').replace('_to_gene','')
    #     for file_path in [os.path.join(dir_path,x) for x in os.listdir(dir_path)]:
    #         org_code = os.path.basename(file_path).split('_')[0]
    #         if org_code in orgcode_to_gnid:
    #             infile = pd.read_csv(file_path, sep='\t', header=None)
    #             for row in infile.to_numpy():
    #                 _, node_ids = row
    #                 for node_id in node_ids.split(';'):
    #                     if file_type == 'pathway':
    #                         node_id = node_id.replace(f":{org_code}", ':map').replace(':','_')
    #                     if file_type == 'module':
    #                         node_id = node_id.replace(f"{org_code}_", "").replace(':','_')

    #                     source_node = orgcode_to_gnid[org_code]
    #                     target_node = f"KEGG:{node_id.replace(':','_')}"
                        
    #                     # Add edge to the knowledge graph
    #                     kg.add_edge(Edge(source_node=source_node, target_node=target_node, predicate=['biolink:genetically_associated_with'], knowledge_source=['KEGG']))
    #                     kg.add_edge(Edge(source_node=target_node, target_node=source_node, predicate=['biolink:genetically_associated_with'], knowledge_source=['KEGG']))
    

    # Save the knowledge graph
    logger.info("Saving the knowledge graph...")
    kg.save_graph(save_dir = args.output_dir, node_filename = 'KG_nodes_v2.tsv', edge_filename = 'KG_edges_v2.tsv')
    logger.info("KG node is saved to {}".format(os.path.join(args.output_dir, 'KG_nodes_v2.tsv')))
    logger.info("KG edge is saved to {}".format(os.path.join(args.output_dir, 'KG_edges_v2.tsv')))
    
    logger.info(f'Done!')