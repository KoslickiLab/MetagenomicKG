# MetagenomicKG
MetagenomicKG is a novel metagenomics knowledge graph which integrates the commonly-used taxonomic information and other biomedical knowledge.

![figure1](https://github.com/KoslickiLab/MetagenomicKG/assets/20031672/22c538a3-7d27-4e79-b07a-35b721314afa)


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
  * [Replicate Use Case2 Sample-specific Graph Embeddings]()
  
## About Graph
MetagenomicKG integrates knowledge from 7 relevant data sources: GTDB taxonomy, NCBI taxonomy, KEGG, RTX-KG2, BV-BRC, MicroPhenoDB, and NCBI AMRFinderPlus Prediction. It consists of 14 node types and 33 edge types (see statistics in the table below).

We have provided the pre-built version of MetagenomicKG, you can download it from [here](https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/cqm5886_psu_edu/EcJBYxxMpvpFsYcWY1ApUaoB2tM11wRHQnZVx1cQpjYY9w?e=qfcYNy). If you would like to rebuild it or reproduce the use case reulsts reported in our paper, you can follow the instruction below.

__Node Statistics__
| **Node Type**      | **Node Count** |
|--------------------|:--------------:|
| Microbe            |     95,2191    |
| Disease            |     10,3129    |
| Phenotypic Feature |     88,708     |
| KO                 |     26,588     |
| Compound           |     19,265     |
| Reaction           |     15,208     |
| Drug               |     12,347     |
| Glycan             |     11,225     |
| Enzyme             |      8,109     |
| AMR                |      5,147     |
| Drug Group         |      2,458     |
| Network            |      1,525     |
| Pathway            |       569      |
| Module             |       481      |
|                    |    1,246,950   |


__Edge Statistics__
| **Edge Type**                            | **Edge Count** |
|------------------------------------------|:--------------:|
| genetically associated with              |   42,134,426   |
| associated with                          |   10,802,017   |
| subclass of                              |    1,141,279   |
| superclass of                            |     981,435    |
| physically interacts with                |     458,541    |
| has participant                          |     192,830    |
| participates in                          |     189,781    |
| related to                               |     95,286     |
| chemically similar to                    |     56,055     |
| biomarker for                            |     47,698     |
| has phenotype                            |     30,291     |
| treats                                   |     21,049     |
| close match                              |     17,389     |
| has part                                 |     13,981     |
| catalyzes                                |      8,819     |
| contraindicated for                      |      5,158     |
| same as                                  |      3,887     |
| is sequence variant of                   |      2,026     |
| gene product of                          |      2,026     |
| correlated with                          |      1,619     |
| contributes to                           |       286      |
| temporally related to                    |       252      |
| affects                                  |       176      |
| causes                                   |       171      |
| has input                                |       82       |
| produces                                 |       27       |
| has metabolite                           |       23       |
| actively involved in                     |       21       |
| gene associated with condition           |        9       |
| disrupts                                 |        4       |
| derives from                             |        1       |
| disease has basis in                     |        1       |
| located in                               |        1       |
|                                          |   56,206,647   |

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

## Replicate Use Case2 Sample-specific Graph Embeddings
To replicate the results of use case 2, you can simply run the Snakemake pipelie via:
```bash
snakemake --cores 16 -s run_usecase2_pipeline.smk targets
```

## Replicate Use Case3 Pathogen Prediction
To replicate the results of use case 3, you can simply run the Snakemake pipelie via:
```bash
snakemake --cores 16 -s run_usecase3_pipeline.smk targets
``` 

## Contact
If you have any questions or need help, please contact @chunyuma or @ShaopengLiu1 or  @dkoslicki.
