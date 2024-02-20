# How to run GTDB_tk to generate assignment data
This folder contains description and all scripts about how to run GTDB_tk to generate assignment data

## Virtual Environment Installation
We recommend using a virtual environment to ensure a clean and isolated workspace for reproducibility. This can be accomplished using either [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [Mamba](https://github.com/mamba-org/mamba) (a faster alternative to Conda).

### Using Conda
To create your Conda environment, follow these steps:

```bash
# Clone the YACHT repository
git clone https://github.com/KoslickiLab/MetagenomicKG.git
cd MetagenomicKG

# Create a new virtual environment named 'gtdb_tk_env'
conda env create -f ../envs/gtdb_tk_env.yml

# Activate the newly created environment
conda activate gtdb_tk_env
```