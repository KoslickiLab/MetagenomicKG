# How to run AMRFinderPlus to predict AMR data
This folder contains descriptions and scripts for running [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) to predict antimicrobial resistance (AMR) data.

## Virtual Environment Installation
We recommend using a virtual environment to ensure a clean and isolated workspace for reproducibility. This can be accomplished using either [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Mamba](https://github.com/mamba-org/mamba) (a faster alternative to Conda).

### Using Conda
To create your Conda environment, follow these steps:

```bash
# Clone the YACHT repository
git clone https://github.com/KoslickiLab/MetagenomicKG.git
cd MetagenomicKG/run_AMRFinderPlus/

# Create a new virtual environment named 'amrfinderplus_env'
conda env create -f ../envs/amrfinderplus_env.yml

# Activate the newly created environment
conda activate amrfinderplus_env

# Install AMRFinderPlus database
amrfinder -u
```

## Predict AMR of Microbe Assemblies used in KEGG
KEGG database uses the genome assemblies from GenBank and Refseq. To predict AMR of those genome assemblies, the commands below first extract the mapping between KEGG genome ids and GenBank/Refseq ids. According to GenBank/Refseq ids, we utilized [NCBI datasets CLI tool](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/) to download the corresponding genome sequences. Then, we ran AMRFinderPlus to predict AMR based on the downloaded genome sequences.

```bash
cd ./KEGG

## extract mapping between KEGG genome ids and GenBank/Refseq ids.
python ./get_KEGG_microbe_assemby_mapping.py

## download genomes and run AMRFinderPlus
nohup bash run.sh 32 &

## combine all prediction results
python combine.py
```

## Predict AMR of Genome Assemblies used in GTDB
We first downloaded a list of 715,230 bacteria genome IDs and 17,245 archaea genome IDs based on the latest version of GTDB taxonomy (rs226) from its database FTP server. Then according to their GenBank/Refseq ids, we also utilized NCBI datasets CLI tool to obtrain the corresponding genome sequences, and ran AMRFinderPlus to predict AMR. The commands below are used to implement these procedures. 
```bash
cd ./GTDB

## download bacteria and archaea genome list from GTDB website
wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/bac120_taxonomy_r226.tsv
wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/ar53_taxonomy_r226.tsv
less ar53_taxonomy_r226.tsv | cut -f 1 | sed 's/RS_//' | sed 's/GB_//' > GTDB_genome_list
less bac120_taxonomy_r226.tsv | cut -f 1 | sed 's/RS_//' | sed 's/GB_//' >> GTDB_genome_list

## download genomes and run AMRFinderPlus
nohup bash run.sh 32 &

## combine all prediction results
python combine.py
```

## Predict AMR of Genome Assemblies used in BVBRC
We first downloaded all genome content (including .fna, .gff, .faa) from BVBRC database FTP server. Then we ran AMRFinderPlus to predict AMR. The commands below are used to implement these procedures. 
```bash
cd ./BVBRC

## Download genomes from BVBRC
wget ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_summary
cat genome_summary | sed '1d' | cut -f 1 > genome_id_list.tsv
cat genome_id_list.tsv | parallel -j 32 wget -r -nH -np -N ftp://ftp.bvbrc.org/genomes/{}/

## download genomes and run AMRFinderPlus
nohup bash run.sh 32 &

## combine all prediction results
python combine.py
```

## Predict AMR of Viruses and fungi in NCBI
We first extracted the NCBI taxonomy of viruses and fungi which has been implemented by step0 in the [build-KG pipeline](https://github.com/KoslickiLab/MetagenomicKG/blob/master/run_buildKG_pipeline.smk).After that we downloaded the corresponding genomes from NCBI and ran AMRFinderPlus to predict AMR. The commands below are used to implement these procedures. 
```bash
cd ./NCBI

## predict AMR for viruses
cd viruses
ln -s ../../../data/Micobial_hierarchy/viruses_hierarchy.tsv
## download genomes and run AMRFinderPlus
nohup bash run.sh 32 &

## predict AMR for fungi
cd fungi
ln -s ../../../data/Micobial_hierarchy/fungi_hierarchy.tsv
## download genomes and run AMRFinderPlus
nohup bash run.sh 32 &
```

## Download Pathogen Detection Reference Gene Catalog
We used the command below to download the reference gene catalog list used by AMRFinderPlus.
```bash
wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/Data/latest/ReferenceGeneCatalog.txt
```
