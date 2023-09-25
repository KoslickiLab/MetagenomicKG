#!/bin/bash

## set up working path
work_folder=$(pwd)

## set up conda environment
# conda create -n pathon_detection python=3.10

## set up a project environment variable
project_name="pathon_detection"
export project_name

## create required folders
if [ ! -d "${work_folder}/data" ]
then
    mkdir ${work_folder}/data
fi
if [ ! -d "${work_folder}/log_folder" ]
then
    mkdir ${work_folder}/log_folder
fi

## download data
# bash ${work_folder}/scripts/data_download.sh ${work_folder}

## process KEGG
# extract organism data from kegg: organism table, organism associated genes, organism associated pathways, organism gene-pathway connections
# python ${work_folder}/scripts/process_data/KEGG/extract_kegg_organism_data.py --organisms 'Archaea' 'Bacteria' 'Fungi' --outdir ${work_folder}/data/KEGG/organisms
# # extract virus data from kegg: virus table, virus associated genes, virus associated pathways, virus gene-pathway connections
# python ${work_folder}/scripts/process_data/KEGG/extract_kegg_virus_data.py --outdir ${work_folder}/data/KEGG/viruses
# # extract other kegg information: different database links and hierarchical relations
# python ${work_folder}/scripts/process_data/KEGG/extract_other_kegg_data.py --outdir ${work_folder}/data/KEGG
# # build a knowledge graph
# python ${work_folder}/scripts/process_data/KEGG/build_KG.py --organisms 'Archaea' 'Bacteria' 'Fungi' --outdir ${work_folder}/data/KEGG
# map KEGG to KG2
# python ${work_folder}/scripts/process_data/KEGG/map_KEGG_to_KG2.py --KG_nodes ${work_folder}/data/KEGG/KG_info/nodes.txt --processes 150 --outdir ${work_folder}/data/KEGG/KG_info/
# # filter gene category in knowledge graph
# python ${work_folder}/scripts/process_data/KEGG/filter_genes.py --kg ${work_folder}/data/KEGG/KG_info/kg.txt --batchsize 500000
# add KG2 information into our KEGG-based knowledge graph
# python  ${work_folder}/scripts/process_data/KEGG/add_KG2_info.py --kg ${work_folder}/data/KEGG/KG_info/kg_gene_filtered.txt \
#                                                                  --kg_nodes ${work_folder}/data/KEGG/KG_info/nodes.txt \
#                                                                  --kegg_dir ${work_folder}/data/KEGG \
#                                                                  --organism_tb ${work_folder}/data/KEGG/organisms/organism_table.txt \
#                                                                  --kg2_nodes ${work_folder}/data/KG2/nodes.txt \
#                                                                  --kg2_edges ${work_folder}/data/KG2/edges.txt \
#                                                                  --outdir ${work_folder}/data/KEGG/KG_info/
# # get ko hierarchy property
# python ${work_folder}/scripts/process_data/KEGG/get_ko_hierarchy.py --outdir ${work_folder}/data/KEGG/KG_info
# # gather genes to ko property
# python ${work_folder}/scripts/process_data/KEGG/gather_genes_per_ko.py --ko ${work_folder}/data/KEGG/kegg_koids.txt \
#                                                                        --organism_tb ${work_folder}/data/KEGG/organisms/organism_table.txt \
#                                                                        --organism_ko_dir ${work_folder}/data/KEGG/organisms/ko_to_gene \
#                                                                        --virus_ko_file ${work_folder}/data/KEGG/viruses/link_ko_to_virus_gene.txt \
#                                                                        --outdir ${work_folder}/data/KEGG/KG_info/
# # add knowledge graph attribute
# python ${work_folder}/scripts/process_data/KEGG/add_kg_attributes.py --kg_nodes_header ${work_folder}/data/KEGG/KG_info/nodes_header.tsv \
#                                                                      --kg_nodes ${work_folder}/data/KEGG/KG_info/nodes.tsv \
#                                                                      --kg_edges_header ${work_folder}/data/KEGG/KG_info/edges_header.tsv \
#                                                                      --kg_edges ${work_folder}/data/KEGG/KG_info/edges.tsv \
#                                                                      --ko_hierarchy ${work_folder}/data/KEGG/KG_info/kegg_ko_brite_info.pkl \
#                                                                      --ko_gene_list ${work_folder}/data/KEGG/KG_info/gene_list_per_ko.pkl \
#                                                                      --PATRIC_pathogen ${work_folder}/data/Pathogen_data/PATRIC/PATRIC_pathogen.txt \
#                                                                      --KEGG_pathogen ${work_folder}/data/KEGG/KG_info/KEGG_pathogen.txt \
#                                                                      --CARD_faa ${work_folder}/data/Pathogen_data/CARD/CARD.faa \
#                                                                      --outdir ${work_folder}/data/KEGG/KG_info/
#import data into neo4j
bash ${work_folder}/scripts/tsv_to_neo4j.sh kegg_kg2_combination.db /data/shared_data/neo4j-community-3.5.26/conf/neo4j.conf ${work_folder}/data/KEGG/KG_info/neo4j_tsv/




# # download sequences
# python ${work_folder}/scripts/process_data/KEGG/download_seq_fasta.py --table ${work_folder}/data/KEGG/viruses/virus_table.txt --col 'rs_ncbi_seq_ids' --outfile ${work_folder}/data/KEGG/viruses/rs_ncbi_virus.fasta
# python ${work_folder}/scripts/process_data/KEGG/download_seq_fasta.py --table ${work_folder}/data/KEGG/viruses/virus_table.txt --col 'gb_ncbi_seq_ids' --outfile ${work_folder}/data/KEGG/viruses/gb_ncbi_virus.fasta
# python ${work_folder}/scripts/process_data/KEGG/download_seq_fasta.py --table ${work_folder}/data/KEGG/organisms/organism_table.txt --col 'rs_ncbi_seq_ids' --organisms 'Archaea' 'Bacteria' 'Fungi' --outfile ${work_folder}/data/KEGG/organisms/rs_ncbi_organism.fasta
# python ${work_folder}/scripts/process_data/KEGG/download_seq_fasta.py --table ${work_folder}/data/KEGG/organisms/organism_table.txt --col 'gb_ncbi_seq_ids' --organisms 'Archaea' 'Bacteria' 'Fungi' --outfile ${work_folder}/data/KEGG/organisms/gb_ncbi_organism.fasta


# ## process data
# python ${work_folder}/scripts/process_data/process_PATRIC.py --in_dir ${work_folder}/data/Pathogen_data/PATRIC --out_dir ${work_folder}/data/Pathogen_data/PATRIC
