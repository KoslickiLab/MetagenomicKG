#!/bin/bash

# make predictions for deepac

# need 3 positional parameters:
# 1. model file
# 2. test_file: genome_id and label
# 3. filepath: genome fastas
# run this pipe inside the folder, it will store results at $PWD

model_file=$1
test_file=$2
filepath=$3

model_file=$(realpath ${model_file})
test_file=$(realpath ${test_file})
filepath=$(realpath ${filepath})

date

# active conda inside script
temp=$(which conda)
conda_path=$(echo ${temp%/*bin/conda})
if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
        . ${conda_path}/etc/profile.d/conda.sh
else
        echo "ERROR: conda path can't be corrected identified!!!"
        exit 1
fi
unset conda_path
conda activate deepac

mkdir -p prediction
cd prediction

# soft link all fasta file here
for f_id in $(cut -f 1 ${test_file}); do
 fasta_file=$(grep ${f_id} ${filepath})
 ln -s ${fasta_file} ./${f_id}.fasta
done

# build npy matric
deepac gwpa fragment -g . -o . -s 250
rm *.fasta

# make prediction
echo -n > merged_deepac_prediction.tsv
for f_id in $(cut -f 1 ${test_file}); do
 deepac predict -c ${model_file} -a ${f_id}_fragmented_genomes.npy -g 0 > output_${f_id}.txt
 prediction=$(grep "Mean prediction:" output_${f_id}.txt | sed 's/Mean prediction: \([0-9.]*\).*/\1/')
 label=$(grep ${f_id} ${test_file} | awk '{print $2}')
 echo -e "${f_id}\t${label}\t${prediction}" >> merged_deepac_prediction.tsv
 unset f_id label prediction
done

rm *npy

echo done
