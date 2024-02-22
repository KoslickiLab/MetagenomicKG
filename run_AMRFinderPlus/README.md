# How to run AMRFinderPlus to predict AMR data
This folder contains descriptions and scripts for running [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) to predict antimicrobial resistance (AMR) data.

## Virtual Environment Installation
We recommend using a virtual environment to ensure a clean and isolated workspace for reproducibility. This can be accomplished using either [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Mamba](https://github.com/mamba-org/mamba) (a faster alternative to Conda).

### Using Conda
To create your Conda environment, follow these steps:

```bash
# Clone the YACHT repository
git clone https://github.com/KoslickiLab/MetagenomicKG.git
cd MetagenomicKG

# Create a new virtual environment named 'amrfinderplus_env'
conda env create -f ../envs/amrfinderplus_env.yml

# Activate the newly created environment
conda activate amrfinderplus_env

# Install AMRFinderPlus database
amrfinder -u
```

## Predict AMR of Microbe Assemblies used in KEGG
KEGG database uses the genome assemblies from GenBank and Refseq. To predict AMR of those genome assemblies, the commands below first extract the mapping between KEGG genome ids and GenBank/Refseq ids. According to GenBank/Refseq ids, we utilize [NCBI datasets CLI tool](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/) to download the corresponding genome sequences. Then, we ran AMRFinderPlus to predict AMR based on the downloaded genome sequences.

```bash
## extract mapping between KEGG genome ids and GenBank/Refseq ids.
python ./get_KEGG_microbe_assemby_mapping.py

## download genomes and run AMRFinderPlus
nohup bash run.sh 32 &

## combine all prediction results
python combine.py
```

## Predict AMR of Genome Assemblies used in GTDB
```bash
## download bacteria and archaea genome list from GTDB website
wget https://data.gtdb.ecogenomic.org/releases/release214/214.0/bac120_taxonomy_r214.tsv
wget https://data.gtdb.ecogenomic.org/releases/release214/214.0/ar53_taxonomy_r214.tsv
less ar53_taxonomy_r214.tsv | cut -f 1 | sed 's/RS_//' | sed 's/GB_//' > GTDB_genome_list
less bac120_taxonomy_r214.tsv | cut -f 1 | sed 's/RS_//' | sed 's/GB_//' > GTDB_genome_list

## download genomes and run AMRFinderPlus
nohup bash run.sh 32 &

## combine all prediction results
python combine.py
```

## Predict AMR of Genome Assemblies used in BVBRC