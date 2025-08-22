# MetagenomicKG
MetagenomicKG is a novel metagenomics knowledge graph that integrates commonly used taxonomic information (GTDB), functional annotations (KEGG), pathogenicity resources (BV-BRC), and other biomedical knowledge. 

The associated preprint can be found at: https://www.biorxiv.org/content/10.1101/2024.03.14.585056v1. Please cite via:
>Ma, C., Liu, S., & Koslicki, D. (2024). MetagenomicKG: A Knowledge Graph for Metagenomic Applications. bioRxiv, 2024-03.

</br> 

![figure1](https://github.com/KoslickiLab/MetagenomicKG/assets/20031672/22c538a3-7d27-4e79-b07a-35b721314afa)


## Table of contents

- [Metagenomic](#metagenomickg)
  * [About Graph](#about-graph)
    + [Download pre-built MKG](#download-pre-built-mkg)
    + [Neo4j instance](#neo4j-instance)
  * [Virtual Environment Installation](#virtual-environment-installation)
    + [Using Conda](#using-conda)
    + [Using Mamba](#using-mamba)
  * [Configuration Setup](#configuration-setup)
    + [UMLS API key](#umls-api-key)
    + [KEGG FTP Data](#kegg-ftp-data)
  * [Build MetagenomicKG](#build-metagenomickg)
  * [Replicate Use Case1 Hypothesis Generation and Exploration](#replicate-use-case1-hypothesis-generation-and-exploration)
  * [Replicate Use Case2 Sample-specific Graph Embeddings](#replicate-use-case2-sample-specific-graph-embeddings)
  * [Replicate Use Case3 Pathogen Predictions](#replicate-use-case3-pathogen-predictions)
  
## About Graph
MetagenomicKG integrates knowledge from 7 relevant data sources: GTDB taxonomy, NCBI taxonomy, KEGG, RTX-KG2, BV-BRC, MicroPhenoDB, and NCBI AMRFinderPlus Prediction. It consists of 14 node types and 33 edge types (see statistics in the table below).

### Download pre-built MKG
We have provided the pre-built version of MetagenomicKG, you can download the `MetagenomicKG.zip` file from [Zenodo](https://zenodo.org/records/10819216). 

### Neo4j instance
We also host a neo4j instance for MetagenomicKG: [http://mkg.cse.psu.edu:7474/](http://mkg.cse.psu.edu:7474/) (Username:`neo4j`, Password:`klabneo4j`). If you would like to rebuild it or reproduce the use case reulsts reported in our paper, you can follow the instruction below.

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

## Configuration Setup
MetagenomicKG uses a centralized configuration system through the `config.yml` file. Before rebuilding MetagenomicKG and replicating use case results, you should configure the system for your environment.

### Configuration File Structure
The `config.yml` file contains several configuration sections:

```yaml
BUILD_KG_VARIABLES:
  # Required: API keys and data directories
  UMLS_API_KEY: 'your_umls_api_key_here'
  KEGG_FTP_DATA_DIR: '/path/to/your/kegg/data'
  
  # Database configuration
  NODE_SYNONYMIZER_DBNAME: 'node_synonymizer_v1.0_KG2.10.0.sqlite'
  NEO4J_DBNAME: 'MetagenomicsKG'
  
  # Neo4j connection (can be overridden by environment variables)
  NEO4J_BOLT: 'bolt://localhost:7687'
  NEO4J_USERNAME: 'neo4j'
  NEO4J_PASSWORD: 'your_neo4j_password'
  
  # Processing thresholds (can be customized based on your needs)
  ANI_THRESHOLD: 99.5
  AF_THRESHOLD: 0.0
  COVERAGE_THRESHOLD: 80
  IDENTITY_THRESHOLD: 90
  
  # KG file names (customizable for different versions)
  KG_FILES:
    NODES_V1: 'KG_nodes_v1.tsv'
    EDGES_V1: 'KG_edges_v1.tsv'
    # ... additional file versions
```

### Required Configuration Changes

#### UMLS API Key
Since we utilize the Unified Medical Language System (UMLS) search function via UMLS APIs for identifier mapping, you must first obtain a UMLS API key:
1. Follow [this instruction](https://documentation.uts.nlm.nih.gov/rest/authentication.html) to get an API key
2. Replace `UMLS_API_KEY` with your actual API key in `config.yml`

#### KEGG FTP Data
MetagenomicKG includes KEGG data downloaded from KEGG FTP. According to KEGG policy, we cannot provide this dataset:
1. Follow [this instruction](https://www.kegg.jp/kegg/download/) to obtain the dataset
2. Replace `KEGG_FTP_DATA_DIR` with the path to your KEGG data directory in `config.yml`

#### Neo4j Connection (Optional)
Configure Neo4j connection parameters in `config.yml` or set as environment variables:
- **Config file approach**: Update `NEO4J_BOLT`, `NEO4J_USERNAME`, and `NEO4J_PASSWORD` in `config.yml`
- **Environment variables approach** (takes precedence over config file):
  ```bash
  export neo4j_bolt="bolt://your-neo4j-server:7687"
  export neo4j_username="your_username"
  export neo4j_password="your_password"
  ```

### Optional Configuration Customization

#### Processing Thresholds
You can adjust processing thresholds based on your requirements:
- `ANI_THRESHOLD`: Average Nucleotide Identity threshold for strain identification (default: 99.5)
- `AF_THRESHOLD`: Alignment Fraction threshold for strain identification (default: 0.0)
- `COVERAGE_THRESHOLD`: Coverage threshold for AMR gene selection (default: 80)
- `IDENTITY_THRESHOLD`: Identity threshold for AMR gene selection (default: 90)

#### Database and File Names
The configuration system allows you to customize:
- Database names (`NODE_SYNONYMIZER_DBNAME`, `NEO4J_DBNAME`)
- KG file names for different processing versions
- Output file names for final Neo4j import

## Build MetagenomicKG
We constructed an automatic pipeline to rebuild MetagenomicKG via [Snakemake](https://snakemake.readthedocs.io/en/stable). Since MetagenomicKG uses [RTX-KG2](https://github.com/RTXteam/RTX-KG2), which includes UMLS data, you need to contact authors to demonstrate that you have accepted [the license terms](https://www.nlm.nih.gov/databases/umls.html) in order to get access to download KG2. We will provide you with a password so you can download the file [here](https://www.dropbox.com/scl/fi/o6458g9ai4ietb4kqp7sx/kg2c-v2.8.4-tsv.tar.gz?rlkey=0bbpesjmz5zct1axt0146xpdg&st=mjaqt9ml&dl=0). **Note**: you do _not_ need this file to view/interact with the [MetagenomicsKG](https://zenodo.org/records/10819216), it's just for if you want to _rebuild_ it. Once you have access, please download and put the `kg2c-tsv.tar.gz` file to `./data/RTX_KG2` folder.

After downloading the RTX-KG2 TSV files is done, you can run the pipeline via:
```bash
snakemake --cores 16 -s run_buildKG_pipeline.smk targets
``` 

Once it is completed, you can find merged node and edge TSV files (`KG_nodes_v6.tsv` and `KG_edges_v6.tsv`) from `./data/merged_KG` folder.

## Replicate Use Case1 Hypothesis Generation and Exploration
Once you login the neo4j instance of MetagenomicKG (see login info in `About Graph` Section), you can use the following Neo4j Cypher Query to replicate what we show in the paper. In this query, we find at most 10 paths of protein - pathogen with name 'Staphylococcus aureus' - KO - pathway - disease - drug.
```
MATCH p=(n0:`biolink:Protein`)-[]-(n1:`biolink:OrganismTaxon`)-[]-(n2:`biolink:BiologicalEntity`)-[]-(n3:`biolink:Pathway`)-[]-(n4:`biolink:Disease`)-[]-(n5:`biolink:Drug`)
WHERE n1.is_pathogen = "True" and ANY (n1_names IN n1.all_names WHERE n1_names contains 'Staphylococcus aureus')
RETURN p LIMIT 10
```

## Replicate Use Case2 Sample-specific Graph Embeddings
To replicate the results of use case 2, you can simply run the Snakemake pipelie via:
```bash
snakemake --cores 16 -s run_usecase2_pipeline.smk targets
```

## Replicate Use Case3 Pathogen Predictions
To replicate the results of use case 3, you can simply run the Snakemake pipelie via:
```bash
snakemake --cores 16 -s run_usecase3_pipeline.smk targets
``` 

## Contact
If you have any questions or need help, please contact @chunyuma or @ShaopengLiu1 or  @dkoslicki.
