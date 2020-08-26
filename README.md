# Pathogen Detection Project

## Table of contents

- [Goal](#goal-of-this-project)
- [Current Idea](#current-idea)


## Goal of This Project
The goal of this project is to develop a algorithemic method to identify the causative miroorganism for diseases via the taxonomic profile generated from some metagenoimcs tools (e.g. [metalign](https://github.com/nlapier2/Metalign) and etc.) and biomedical knowledge graph (e.g. kg2) developed by [ARAX](https://github.com/RTXteam/RTX) project.  

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
