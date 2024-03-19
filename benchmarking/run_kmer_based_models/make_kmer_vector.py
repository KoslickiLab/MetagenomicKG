from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import itertools, os, sys
import pandas as pd

if len(sys.argv) < 5:
        print("Usage: python script.py filepath_list k_length input_data_w_label out_name")
        sys.exit(1)

target_file = sys.argv[1]
k_length = int(sys.argv[2])
label_file = sys.argv[3]
out_name = sys.argv[4]

def generate_kmers(k):
        bases = ['A', 'C', 'G', 'T']
        kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
        return list(set(min(kmer, str(Seq(kmer).reverse_complement())) for kmer in kmers))

def compute_kmer_frequency(fasta_file, kmer_size, kmer_list):
        kmer_counts = Counter()
        for record in SeqIO.parse(fasta_file, "fasta"):
                sequence = record.seq
                for i in range(len(sequence) - kmer_size + 1):
                        kmer = str(sequence[i:i+kmer_size]).upper()
                        canonical_kmer = min(kmer, str(Seq(kmer).reverse_complement()))
                        kmer_counts[canonical_kmer] += 1
        return [kmer_counts[kmer] for kmer in kmer_list]


with open(target_file, 'r') as file:
    fna_files = [line.strip() for line in file]

out_kmer_list = generate_kmers(k_length)
temp_out = []  # record each vector and then merge them

# loop through files
for file in fna_files:
        print(file)
        temp_vector = compute_kmer_frequency(file, k_length, out_kmer_list)
        temp_out.append(temp_vector)

# write to csv
df = pd.DataFrame(temp_out, index=fna_files, columns=out_kmer_list)
df.index = [os.path.basename(path) for path in df.index]
df.index = df.index.map(lambda x: '_'.join(x.split('_')[:2]))
df.index.name="filename"


# add label information
label = pd.read_csv(label_file, header=None, sep="\t", index_col=0, names=['label'])
out_df = pd.merge(df, label, left_index=True, right_index=True, how='inner')
out_df.index.name="filename"
out_df.to_csv(out_name, index=True)
