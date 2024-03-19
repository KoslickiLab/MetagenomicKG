import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from gensim.models import Word2Vec

if len(sys.argv) < 4:
    print("Usage: python script.py train_file_paths k_length model")
    sys.exit(1)

train_file_paths = sys.argv[1]
k_length = int(sys.argv[2])
out_name = sys.argv[3]

fasta_files = []
with open(train_file_paths, 'r') as file:
    for line in file:
        # You can use line.strip() to remove newline characters
        fasta_files.append(line.strip())

# word2vec relies on context, so the relative order matters, we need stream all, but not list all possible
def kmerize(sequence, k):
    out_kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k].upper()
        canonical_kmer = min(kmer, str(Seq(kmer).reverse_complement()))
        out_kmers.append(canonical_kmer)
    return out_kmers


# Prepare data for Word2Vec
sentences = []  # This will hold all k-mers across all genomes
for file_path in fasta_files:
    print(file_path)
    for record in SeqIO.parse(file_path, 'fasta'):
        sentences.append(kmerize(str(record.seq), k_length))  # each contig becomes a sentence

# Train a Word2Vec model
model = Word2Vec(sentences, vector_size=300, window=10, min_count=3, workers=16)

# Save the model
model.save(out_name)
