# MetagenomicKG
MetagenomicKG is a novel metagenomics knowledge graph which integrates the commonly-used taxonomic information and other biomedical knowledge.

## Table of contents

- [Metagenomic](#metagenomickg)
  * [About Graph](#about-graph)
  * [Virtual Environment Installation](#virtual-environment-installation)
    + [Using Conda](#using-conda)
    + [Using Mamba](#using-mamba)
  * [Modify `config.yaml` File If Needed ](#modify-configyaml-file-if-needed)
    + [UMLS API key](#umls-api-key)
    + [KEGG FTP Data](#kegg-ftp-data)
  * [Build MetagenomicKG](#build-metagenomickg)
  * [Use Case1 Pathogen Identification](#use-case1-pathogen-identification)
  
## About Graph
MetagenomicKG integrates knowledge from 7 relevant data sources: GTDB taxonomy, NCBI taxonomy, KEGG, RTX-KG2, BV-BRC, MicroPhenoDB, and NCBI AMRFinderPlus Prediction. It consists of 13 node types and 35 edge types (see statistics in the table below).

We have provided the pre-built version of MetagenomicKG, you can download it from [here](). If you would like to rebuild it or reproduce the use case reulsts reported in our paper, you can follow the instruction below.

__Node Statistics__
| **Node Type**             | **Node Count** |
|---------------------------|----------------|
|          Compound         |      19,217    |
|           Disease         |     103,719    |
|            Drug           |      12,004    |
|         Drug Group        |      2,429     |
|           Enzyme          |      8,056     |
|           Glycan          |      11,218    |
|             KO            |      25,745    |
|           Microbe         |     828,660    |
|           Module          |       475      |
|           Network         |      1,417     |
|           Pathway         |       562      |
|     Phenotypic Feature    |     106,754    |
|          Reaction         |      15,130    |
|                           | 1,135,386      |


__Edge Statistics__
| **Edge Type**                  | **Edge Count** |
|--------------------------------|----------------|
| associated with                |      4,494,926 |
| has part                       |       867,991  |
| part of                        |       854,022  |
| physically interacts with      |       443,666  |
| has participant                |       190,447  |
| participates in                |     187,536    |
| subclass of                    |       169,125  |
| related to                     | 96,210         |
| genetically associated with    | 62,096         |
| chemically similar to          | 53,800         |
| biomarker for                  | 44,143         |
| has phenotype                  | 28,875         |
| treats                         | 20,369         |
| interacts with                 | 19,768         |
| close match                    | 17,239         |
| catalyzes                      | 8,752          |
| contraindicated for            | 5,173          |
| same as                        | 3,891          |
| is sequence variant of         | 1,918          |
| gene product of                | 1,918          |
| correlated with                | 1,590          |
| contributes to                 | 238            |
| temporally related to          | 229            |
| causes                         | 202            |
| affects                        | 119            |
| has input                      | 81             |
| regulates                      | 40             |
| produces                       | 27             |
| actively involved in           | 25             |
| has metabolite                 | 23             |
| gene associated with condition | 12             |
| located in                     | 1              |
| derives from                   | 1              |
| disease has basis in           | 1              |
| disrupts                       | 1              |
|                                |                |

## Virtual Environment Installation
We recommend using a virtual environment to ensure a clean and isolated workspace for reproducibility. This can be accomplished using either [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Mamba](https://github.com/mamba-org/mamba) (a faster alternative to Conda).

### Using Conda
To create your Conda environment, follow these steps:

```bash
# Clone the MetagenomicKG repository
git clone https://github.com/KoslickiLab/MetagenomicKG.git
cd MetagenomicKG

# Please note that the versioin of pytorch we used might not be compatible with your nvidia cuda version. So, please first check your version and change it in metagenomickg_env.yml if needed.
conda env create -f envs/metagenomickg_env.yml

# Download required files for the package 'pytaxonkit'
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz

mkdir -p $HOME/.taxonkit
cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/.taxonkit

# Activate the newly created environment
conda activate metagenomickg_env
```

### Using Mamba
If you prefer using Mamba instead of Conda, just simply repalce `conda` with `mamba` in the above commands.

## Modify `config.yaml` File If Needed 
Before rebuilding MetagenomicKG and replicating the results of use cases, you can should modify some global variables in the `config.yaml` file. We listed some required variables below: 

### UMLS API key
Since we utilize the Unified Medical Language System (UMLS) search function via UMLS APIs for identifier mapping, you should first get a UMLS API key (please follow [this instruction](https://documentation.uts.nlm.nih.gov/rest/authentication.html) to get one). After that, replace `UMLS_API_KEY` with your API key in the `config.yaml` file.

### KEGG FTP Data
MetagenomicKG includes KEGG data downloaded from KEGG FTP. According to KEGG policy, we can provide this dataset. To obtain this dataset, you can follow [this instruction](https://www.kegg.jp/kegg/download/). Once you download data, you should replace the path of `KEGG_FTP_DATA_DIR` key in the `config.yaml` file.

## Build MetagenomicKG
We constructed an automatic pipeline to rebuild MetagenomicKG via [Snakemake](https://snakemake.readthedocs.io/en/stable). Since MetagenomicKG uses [RTX-KG2](https://github.com/RTXteam/RTX-KG2), which includes UMLS data, you need to contact authors to demonstrate that you have accepted [the license terms](https://www.nlm.nih.gov/databases/umls.html) in order to get access to download KG2. Once you have access, please download and put the `kg2c-tsv.tar.gz` file to `./data/RTX_KG2` folder.

After downloading the RTX-KG2 TSV files is done, you can run the pipeline via:
```bash
snakemake --cores 16 -s run_buildKG_pipeline.smk targets
``` 

Once it is completed, you can find merged node and edge TSV files (`KG_nodes_v6.tsv` and `KG_edges_v6.tsv`) from `./data/merged_KG` folder.

## Use Case1 Pathogen Identification
To replicate the results of use case 1, you can simply run the Snakemake pipelie via:
```bash
snakemake --cores 16 -s run_usecase1_pipeline.smk targets
``` 

## Contact
If you have any questions or need help, please contact @chunyuma or @ShaopengLiu1 or  @dkoslicki.
