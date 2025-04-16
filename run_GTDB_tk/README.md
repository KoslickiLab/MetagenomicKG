# How to run GTDB_tk to generate taxonomic assignment
This folder contains description and all scripts about how to run [GTDB_tk](https://github.com/Ecogenomics/GTDBTk) to assign taxonomic labels to KEGG and BV-BRC genomes. For convenience, the results for KEGG and BV-BRC pathogen genomes can be downloaded at [Zenodo](https://zenodo.org/records/10806617) in file `taxonomy_assignment_by_GTDB_tk.tar.gz`.

</br>

## Virtual Environment Installation

We recommend using a virtual environment to ensure a clean and isolated workspace for reproducibility. This can be accomplished using either [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Mamba](https://github.com/mamba-org/mamba) (a faster alternative to Conda).

### Using Conda
To create your Conda environment, follow these steps:

```bash
# Clone the YACHT repository
git clone https://github.com/KoslickiLab/MetagenomicKG.git
cd MetagenomicKG/run_GTDB_tk/

# Create a new virtual environment named 'gtdb_tk_env', current version is 2.4.0
conda env create -f ../envs/gtdb_tk_env.yml

# Activate the newly created environment
conda activate gtdb_tk_env
```

</br>

### Download reference database

---

1. this process can take several hours. So it's recommended to run in background.

2. `download-db.sh` is a built-in script in GTDB-tk that downloads a specific version of GTDB database. You can check or edit this version by modifying the links inside this script. By default, our analysis utilizes GTDB-tk v2.4.0 and GTDB r220.

```
# this will download ~100GB ref data into local disk and takes long time to run
nohup download-db.sh & 
```

</br>

### Running GTDB-tk with batch file

---

GTDB-tk can process all genomes withing a directory or recorded in a batch file, [check here](https://ecogenomics.github.io/GTDBTk/commands/classify_wf.html) for more details. We suggest the latter way for convience. A batch file is a tab deliminated file with two column file indicating the location of each genome and the desired genome identifier.

```
# process all genomes
gtdbtk classify_wf --batchfile <batchfile> --out_dir <out_dir> --skip_ani_screen
```

</br>

### Process all KEGG and BV-BRC genomes

1. download BV-BRC pathogen genomes from its [FTP server](https://www.bv-brc.org/docs/quick_references/ftp.html)
2. KEGG genomes are obtained from its FTP server as well, but this step requires KEGG's license
3. prepare a batch file (`gtdb-tk_input_batchfile.tsv`) with absolute paths of genome files as 1st column, and arbitrary file id as 2nd colum
4. to accelerate this process, we use the `split` command with `GNU parallel` as following:

```
# go to an empty dir where you want to store all the results

# run the bash script:
bash <path-2-parallel_gtdbtk.sh> <batchfile> 4 8 
```

</br>

### Merge all results together

---

```
head -1 $(find . -maxdepth 2 -name "gtdbtk.*.summary.tsv" | head -1) | cut -f 1-12,14,17 > merged_all.tsv

for output_file in $(find . -maxdepth 2 -name "gtdbtk.*.summary.tsv"); do
 wc -l ${output_file}
 cut -f 1-12,14,17 ${output_file} | sed '1d' >> merged_all.tsv
done
```



