import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from gensim.models import Word2Vec
import numpy as np
import pandas as pd


if len(sys.argv) < 5:
    print("Usage: python script.py query_file_path model_file k_length out_csv_name")
    sys.exit(1)

input_file_paths = sys.argv[1]
model_file = sys.argv[2]
k_length = int(sys.argv[3])
out_csv = sys.argv[4]

# load model
model = Word2Vec.load(model_file)

# load input files
fasta_files = []
with open(input_file_paths, 'r') as file:
    for line in file:
        # You can use line.strip() to remove newline characters
        fasta_files.append(line.strip())


def kmerize(sequence, k):
    out_kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k].upper()
        canonical_kmer = min(kmer, str(Seq(kmer).reverse_complement()))
        out_kmers.append(canonical_kmer)
    return out_kmers


def get_kmer_list_of_a_genome(file, k):
    out_list = []
    for record in SeqIO.parse(file, 'fasta'):
        out_list.extend(kmerize(str(record.seq), k))
    return out_list

def get_genome_vector(model, kmer_list_from_genome):
    """Get the vector representation of a genome."""
    vectors = []

    for k_mer in kmer_list_from_genome:
        if k_mer in model.wv:
            vectors.append(model.wv[k_mer])
        else:
            # Skip the k-mer if it's not in the model
            continue

    # Return the mean of the vectors
    if vectors:
        return np.mean(vectors, axis=0)
    else:
        return np.zeros(model.vector_size)


# build df
row_name=[]
row_content=[]
for file in fasta_files:
    row_name.append(os.path.basename(file))
    temp_kmer_list = get_kmer_list_of_a_genome(file, k=k_lengt
h)
    temp_vector = get_genome_vector(model, temp_kmer_list)
    row_content.append(temp_vector)


# write to csv
df = pd.DataFrame(row_content, index=row_name)
df.index.name = "filename"
df.to_csv(out_csv, index=True)
