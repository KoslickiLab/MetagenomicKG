# MetagenomicKG
MetagenomicKG is a novel metagenomics knowledge graph which integrates the commonly-used taxonomic information and other biomedical knowledge from 7 relevant data sources: GTDB taxonomy, NCBI taxonomy, KEGG, RTX-KG2, BV-BRC, MicroPhenoDB, and NCBI AMRFinderPlus Prediction. 

## Environment Installation
We recommend using a virtual environment to ensure a clean and isolated workspace for reproducibility. This can be accomplished using either [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Mamba](https://github.com/mamba-org/mamba) (a faster alternative to Conda).

#### Using Conda
To create your Conda environment, follow these steps:

```bash
# Clone the YACHT repository
git clone https://github.com/KoslickiLab/MetagenomicKG.git
cd MetagenomicKG

# Create a new virtual environment named 'metagenomickg_env'
conda env create -f env/metagenomickg_env.yml

# Activate the newly created environment
conda activate metagenomickg_env
```

#### Using Mamba
If you prefer using Mamba instead of Conda, just simply repalce `conda` with `mamba` in the above commands.



## Table of contents

- [Goal](#goal-of-this-project)
- [MKG](#about-microbial-knowledge-graph)
- [Current Idea](#current-idea)

## Goal of This Project
The goal of this project is to develop a algorithemic method to identify the causative miroorganism for diseases via the taxonomic profile generated from some metagenoimcs tools (e.g. [metalign](https://github.com/nlapier2/Metalign) and etc.) and biomedical knowledge graph (e.g. kg2) developed by [ARAX](https://github.com/RTXteam/RTX) project.  

## About Microbial Knowledge Graph
The Microbial knowledge graph is a directed graph that currently includes knowledge from five different sources (e.g., [KEGG](https://www.genome.jp/kegg/), [RTX-KG2](https://github.com/RTXteam/RTX-KG2), [BV-BRC](https://www.bv-brc.org/), [MicroPhenoDB](http://lilab2.sysu.edu.cn/microphenodb/#/home). It consists of 13 node types and 35 edge types (see statistics in the table below). For mirobe nodes, we integrate the taxonomy hierarchy of all bacteria and archaea from GTDB, as well as the hierarchy of fungi and viruses from NCBI. 

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


### How to Build
We built an automatic pipeline to build this microbial knowledge graph via [Snakemake](https://snakemake.readthedocs.io/en/stable).

Please follow the steps below to run this pipeline:

1. install conda and then run the following commands to setup enironment
```bash
conda env create -f envs/mkg_envs.yml

## activiate the 'mkg_envs' conda environment
conda activate mkg_envs.yml
```

2. run pipeline via the following command
```bash
snakemake --cores 16 -s run_pipeline.smk targets
``` 


## Current Idea:
1. Novel microbe-disease/phenotypic feature prediction
  * Use GraphSage to convert the microbe and disease/phenotypic feature nodes in KG2 to embedding vectors
  * Use PCA/t-SNE to reduce dimension of embedding vectors and cluster some similar nodes
  * Use the known microbe-disease/phenotypic feature relatinship to find novel association

2. Finding known pathogens in people with/witout disease
* __input__: EHR, meta tax profile
* __challenges__:
  * Define some metric to find causative pathogen
  * How to obtain input data
    * collaborate with clinicians in Hersey 
    * case studies + publications 
      * manually take a look
      * find pubd id in semdedb that contains both microbe + sysmpytoms 
      * Sequence databases, Qitta [16s rRNA], SRA at NCBI WGS metadata)
