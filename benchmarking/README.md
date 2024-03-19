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

5. Usecase 3: baseline models for pathogen detection

   1. [PaPrBaG](#paprbag)
   2. [DeePac](#deepac)
   3. [RF and SVM with 4-to-6-mer frequencies](#rf)




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



</br>

## Usecase 3: pathogen detection

Our training data contains bacterial genomes from HMP healthy samples (negative) and known pathogens from KEGG and BV-BRC (positive). They are separated into 3 subgroups for model training and testing purpose:

1. training set: genomes and labels are used for model training
2. validation set: genomes and labels are used for internal model validation during training cycles
3. testing set: serves as external data to evaludate model performance by comparing predictions results with known labels

We intentionally applied 2 data split methods:

1. regular 10-folder cross validation: genomes are randomly assigned to each set in a rolling basis
2. species maskout: the testing set contains genomes that are NOT close (measured by low ANI) to training genomes

The evaluation matric used is:

```
 Accuracy (ACC), 
 Positive Accuracy (PACC), 
 Negative Accuracy (NACC), 
 Area Under the Receiver Operating Characteristic Curve (AUROC), 
 Average Precision (AP), 
 F1 score (F1) 
```

</br>

### Pathogen prediction by PaPrBag <a name="paprbag"></a>

---

[PaPrBaG](https://github.com/crarlus/paprbag) is a R package that would first extract sequence-based features out of input reads or contigs and then apply a random forest classifer for pathogen prediction on read or contig level. The genome-wide prediciton is achieved by a majority vote from all contigs within genome. Here we give an example of processing cluster maskout data, but we can apply to arbitrary genome lists by switching the input nodes file. 

#### Installation in R

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
install.packages("devtools")
devtools::install_github("crarlus/paprbag")

# it seems R can't properly handle the environment for "devtools", 
# use "conda install conda-forge::r-devtools" is an alternative way to install it inside an environment
```

</br>

#### Extract features from sequence data

1. Input data: our MKG nodes file used in the graph model contains all necessary information. Here we only need its column 1 (file id) and column 7 (pathogen status)
   1. genome list: 2-column tsv file with `file identifier` and `label`. For example: `GCA_002404795.1	False`; they are split into 3 subfiles: training, validation (not used here), and test
   2. filepath: a txt file containing fasta file paths of all genomes (downloaded from NCBI by file id). Genomes can be downloaded by `datasets download` and then use `realpath` to get its path.
2. Output:
   1. dir `PaPrBag_output_cluster_mask_${time}`: a dir containing intermediate outputs for PaPrBaG
   2. Under the directory above:
      1. `training.tsv`: a "|" deliminated table of training data ready to train a RF classifer
      2. `validation.tsv`: a "|" deliminated table of test data (the name is confusing, though) ready test the RF classifer

```
# use cluster maskout sample as example, we can dereive 3 tsv file storing nodes information in the graph
cut -f 1,7 train_set.tsv | sed '1d' > label_cluster_mask_train.tsv
cut -f 1,7 test_set.tsv | sed '1d' > label_cluster_mask_test.tsv
cut -f 1,7 valid_set.tsv | sed '1d' > label_cluster_mask_valid.tsv

# we only need genome id and label here
training_list=label_cluster_mask_train.tsv
test_list=label_cluster_mask_test.tsv

# prepare genomes into proper format for PaPrBaG (get contig-level labels and merge all input contigs)
bash <path-2-run_paprbag/prepare_genome_files_as_input.sh> ${training_list} ${test_list} <filepath_w_absolute_paths_of_all_genomes>  output_cluster_mask

# run PaPrBaG to extract features
cd PaPrBag_output_cluster_mask_${time_tag}
Rscript <path-2-run_paprbag/paprbag_extract_feature.R>
```

</br>

#### Train a RF classifier and examine the results

In our benchmarking analysis, we found the default `evaluation` function in R didn't give a fixed output label order. For a consistent prediction and for uniform performance calculation, we moved all the downstream analysis in Python. Output:

1. `py_contig_level_RF_prediction.tsv`: this is the direct output from PaPrBaG that predicts pathogen status on contig level
2. `py_RF_performance.tsv`: this is the genome-level results summary for this input sample

```
# continue of previous part
conda activate metagenomickg_env

# train a RF classifier and make contig-level prediction
python <path-2-run_paprbag/run_py_RF.py>

# aggregate contig-level predictions into genome-level, and get performance
python <path-2-run_paprbag/get_py_RF_performance.py>
```



</br>

### Pathogen prediction by DeePac <a name="deepac"></a>

---

[DeePac](https://gitlab.com/rki_bioinformatics/DeePaC) is a python package and a CLI tool for predicting labels (e.g. pathogenic potentials) from short DNA sequences (e.g. Illumina reads) with interpretable reverse-complement neural networks. Similar to PaPrBag, it gives fragment-level predictions and then do a majority vote for genome prediction by average probability. Here we give an example of processing cluster maskout data, but we can apply to arbitrary genome lists by switching the input nodes file. 

#### Environment setup

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda create -n deepac python=3.9
conda activate deepac
conda install deepac
conda install anaconda::ipython

# need an older version of numpy (with depreciated method 'typeDict')
pip install numpy==1.22.4
```

</br>

#### Prepare training npy data for DeePac

The input data is same as PaPrBaG, we need:

1. filepath: a txt file containing fasta file paths of all genomes (downloaded from NCBI by file id). Genomes can be downloaded by `datasets download` and then use `realpath` to get its path.
2. genome list: 2-column tsv file with `file identifier` and `label`. For example: `GCA_002404795.1	False`; they are split into 3 subfiles: training, validation, and test. They can be extracted from the nodes file used in graph models (col1 and col7)

**Warning:**

1. it's highly recommended to train the model on GPU, o.w. additional computational cost and model transformation are required 
2. the training data is **EXTREMELY large** because DeePac:
   1. chop genomes into overlapping fragments of length 250 (for efficiency, we set this overlap to 0)
   2. generate 1-hot encoding for ALL fragments
   3. for example, a 60MB fasta file will be transformed in to 3GB training data in npy format in our testing
   4. therefor, the deepac process is MEM-extensive (expect max 500GB for the analysis below)
3. **Please specify resouce usage for model training** in the model training part, parameters:
   1. `-r` for rapid CNN mode, use `-s` for sensitive LSTM mode
   2. `-n` thread number
   3. `-g` specify GPU to use, default ALL

```
# the input data is the same as PaPrBaG, we just need a files containing genome id and label
training_list=label_cluster_mask_train.tsv
valid_list=label_cluster_mask_valid.tsv 
test_list=label_cluster_mask_test.tsv

# prepare input data (time consuming and MEM extensive)
nohup bash <path-2-run_deepac/prepare_train_valid_metric_from_input_lists.sh> ${traiing_list} ${valid_list} cluster_mask <filepath file> <path-2-run_deepac/build_training_metric_from_fragments.py> &

# model training (time consuming, computationally expensive)
conda activate deepac
cd DeePac_cluster_mask_${time_tag}/input_npy
nohup deepac train -r -T train.npy -t label_train.npy -V valid.npy -v label_valid.npy -R model -n 96 -g 2 &


# make prediction
# models and training outputs are stored in input_npy/logs/model-logs, you may want to pick different one than the last round of training
model_file=$(realpath ./input_npy/logs/model-logs/model-e015.h5)

conda activate deepac
bash <path-2-run_deepac/make_prediction_for_test_data_with_selected_model.sh> ${model_file} ${test_list} <filepath file>


# summarize results
# go to the DeePac_cluster_mask_${time_tag} folder where you can find "input_npy" dir
python <path-2-run_deepac/get_deepac_performance.py> <path-2-run_deepac/metrics.py>
```

</br>

### Pathogen prediction by 4-to-6-mer frequencies <a name="rf"></a>

---

In the PaPrBaG manuscript, the authors benchmarked on lots of feature combinations and found that the 4-mer frequencies give best discriminative power. Similarly, 4-mer frequencies are widely used for contig binning in metagenomic analysis. In order to get whole-genome-level prediction comparisons, we collected 4-to-6-mer frequency information and use them to train RF and SVM models.

#### Environment setup

```
conda create -n kmer_model
conda activate kmer_model
conda install conda-forge::biopython
conda install anaconda::scikit-learn
conda install conda-forge::gensim
```



</br>

#### Analysis

The input data is same as PaPrBaG, we need:

1. filepath: a txt file containing fasta file paths of all genomes (downloaded from NCBI by file id). Genomes can be downloaded by `datasets download` and then use `realpath` to get its path.
2. genome list: 2-column tsv file with `file identifier` and `label`. For example: `GCA_002404795.1	False`; they are split into 3 subfiles: training, validation, and test. They can be extracted from the nodes file used in graph models (col1 and col7)

Output data:

1. `py_genome_level_RF_prediction.tsv`: genome level prediction probability with true label
2. `py_genome_level_SVM_prediction.tsv`: genome level prediction probability with true label

```
# go to the dir where you want to store results
conda activate kmer_model

# the input data is the same as PaPrBaG, we just need a files containing genome id and label
training_list=label_cluster_mask_train.tsv
test_list=label_cluster_mask_test.tsv
filepath=<path-2-all-fasta-genomes>


############ RF with 4-mer abundance
# build kmer abundance vector (last column "label")
cut -f 1 ${training_list} | grep -f - ${filepath} > temp_filepath_training.txt
cut -f 1 ${test_list} | grep -f - ${filepath} > temp_filepath_test.txt
### vector for training data
# 2nd parameter is k_length
python <path-2-run_kmer_based_models/make_kmer_vector.py> temp_filepath_training.txt 4 ${training_list} k4_vector_training.csv
python <path-2-run_kmer_based_models/make_kmer_vector.py> temp_filepath_test.txt 4 ${test_list} k4_vector_test.csv

# make prediction
python <path-2-run_kmer_based_models/train_RF_and_predict.py> k4_vector_training.csv k4_vector_test.csv ${test_list}

# then we can directly summarize results as it's genome-level


############# SVM with 6-mer embedding
# 1st, we need to train Word2vec models based on input sequences
# the 2nd parameter is for k length, used 6 here
cat temp_filepath_training.txt temp_filepath_test.txt > merged_files_for_word2vec.txt
python <path-2-run_kmer_based_models/word2vec_train_model.py> merged_files_for_word2vec.txt 6 out_word2vec_model

# generate sample specific embeddings:
python <path-2-run_kmer_based_models/word2vec_generate_embedding.py> merged_files_for_word2vec.txt out_word2vec_model 6 both_pos_neg_embeddings.csv

# after adding labels, we can do similar model training/prediction as the previous RF model
python <path-2-run_kmer_based_models/train_SVM_and_predict_based_on_word2vec_embeddings.py> both_pos_neg_embeddings.csv ${training_list} ${test_list} 
```

