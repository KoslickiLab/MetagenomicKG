## Benchmarking analysis

This folder contains descriptions and scripts for obtaining and processing data from the Human Microbiome Project (HMP) for Usecase 2; and the baseline model for pathogen detection for Usecase 3.

</br>

### Contents

1. [Download HMP datasets](#hmp)

2. [Build sourmash sketches](#sketch)

3. Usecase 2: graph-embeddings for metagenomic samples

   1. [Get taxonomic profiles of HMP samples by Sourmash gather](#profile)

4. Usecase 3: data preparation

   1. [Identify non-pathogen genomes from healthy samples by YACHT](#negative)
   2. [Generate ANI matrix for data split](#ani)

5. Usecase 3: baseline models

   1. PaPrBaG
   2. DeePac
   3. Random Forest with 4-to-6-mer frequencies
   4. SVM with 4-to-6-mer embeddings

   



</br>

### Download HMP datasets <a name="hmp"></a>

---

All human metagenomic samples analyzed in this study were obtained from the Human Microbiome Project (HMP) [data portal](https://portal.hmpdacc.org/search/s?facetTab=cases), utilizing the search parameters “fastq” for format and “WGS raw seq set” for type. For reproducibility, we provides our original metadata from HMP at [Zenodo](https://zenodo.org/records/10806617) in file `other_data.zip`. Outcome:

1. a folder named `checked_fq_files` storing all downloaded raw sequence data
2. `filepath_all_bz2.txt` a txt file storing all absolute paths of those raw data

```
# go to an empty folder where you want to download the data, and then download the 2 metadata from Zenodo

# clean metadata
bash <path-2-download_HMP/clean_metadata.sh>

# download all samples
mkdir downloaded_data
cd downloaded_data
# it takes a long time to download them
nohup bash <path-2-download_HMP/download_hmp_datasets.sh> & 

# delete files with size <200MB
find . -maxdepth 1 -type f -size -200M -exec rm {} \;

# get filepaths of all data
realpath ./checked_fq_files/*.bz2 > filepath_all_bz2.txt
```

</br>



### Build sourmash sketches <a name="sketch"></a> 

---

We use [Sourmash](https://sourmash.readthedocs.io/en/latest/) to build k-mer sketches with parameter `k=31, abund=1000` for all datasets for downstream analysis. Output:

1. dir `single_sketch`: a folder storing individual k-mer sketches
2. dir `runlog`: running log files for individual data

```
# go to an emtpy folder where you want to save these k-mer sketches, then copy the previous "filepath_all_bz2.txt" here as input.

# the only required parameter is parallel job number: 16
nohup bash <path-2-build_sourmash_sketch/build_sketch_for_hmp_bz2_data.sh> 16 &
```



</br>

### Get taxonomic profiles by Sourmash gather <a name="profile"></a>

---

There are 1675 bz2 files being downloaded, we will then use `sourmash gather` to taxonomically profile them based on GTDB r214. Output:

1. dir `single_gather_out`: a folder storing all individual profiling results
2. dir `runlog_single_gather`: a folder storing all individual run logs
3. CSV file `merged_funiqweighted_fillna0.csv`: merged relative abundance of taxonomic profiling results from all samples

```
# continue from the previous step

# prepare ref data: download pre-built GTDB r214 k31 reference from Sourmash
wget https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-k31.sbt.zip

# prepare query data: 
realpath ./single_sketch/*.zip > filepath_hmp_sketches.txt

# generate taxonomic profile, require parallel number: 16
nohup bash <path-2-taxonomic_profile/run_sourmash_gather_for_profile.sh> 16 &

# merge all outputs
cd single_gather_out
conda activate metagenomickg_env
python <path-2-taxonomic_profile/merge_output.py> 
```





</br>

### Identify non-pathogen genomes from healthy samples by YACHT  <a name="negative"></a> 

---

We use [YACHT](https://academic.oup.com/bioinformatics/article/40/2/btae047/7588873), a state-of-the-art metagenomic detection tool, to identify microbial genomes from healthy samples. These detected genomes are labelled as non-pathogen and used as negative samples in Usecase3: pathogen detection. Output:

1. dir `yacht_output`: a folder storing all individual YACHT output
2. dir `yacht_runlog`: a folder storing all individual YACHT run log
3. merged_genomes.txt: merged lists of genomes identified by YACHT

```
# go to an empty dir where you want to store your results, and copy the previous "filepath_hmp_sketches.txt" file here as input


# YACHT setup: download pre-built reference database based on GTDB. This may take a while depends on your internet
wget https://zenodo.org/records/10113574/files/gtdb-rs214-reps.k31_0.9995_pretrained.zip


# run YACHT to identify genomes present in HMP samples with 95% confidence, require parallel number: 16
nohup bash <path-2-yacht_detection/yacht_run.sh> 16 &

# generate a merged list
cd yacht_output
conda activate metagenomickg_env
python <path-2-yacht_detection/merge_yacht_output.py>
mv merged_genomes_*.txt ../
cd ..
```

</br>



### Get ANI matrics for training data to setup random and species-mask data split <a name="ani"></a> 

---

Besides the non-pathgen genomes identified from healthm HMP samples, we also added known pathogens from KEGG and BV-BRC as positive data in our training datasets for pathogen detection. We then use Sourmash to estimate pariwise ANI matrix based on which to perform clustering analysis. 

```
# go to an empty dir where you want to store your results

# we need to prepare 2 files here
# 1. "filepath_positive_data.txt": filepaths of all positive samples from KEGG and BV-BRC 
# 2. "filepath_negative_data.txt": filepaths of all negative smaples identified above

conda activate metagenomickg_env

# build sketch
cat filepath_positive_data.txt filepath_negative_data.txt > filepath_merged.txt
sourmash sketch dna -o mkg_posi_and_nega_genomes.sig.zip -p k=21,abund --from-file filepath_merged.txt

# generate ANI matrics
sourmash compare -k 21 --csv ani_matrix_for_all_genomes.csv --estimate-ani mkg_posi_and_nega_genomes.sig.zip

# prepare label data for plot
awk '{print $1"\t""nonp"}' filepath_negative_data.txt > file_label.txt
awk '{print $1"\t""patho"}' filepath_positive_data.txt >> file_label.txt

# make heatmap
python <path-2-data_split/make_heatmap.py>
```







