#! /bin/bash

# run YACHT to detect genomes from HMP samples with 95% confidence. May adjust this cutoff if needed.

parallel_num=$1 # number of jobs to run in parallel

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

conda activate metagenomickg_env


# ref data
json_file=$(realpath ./gtdb-rs214-reps.k31_0.9995_pretrained/gtdb-rs214-reps.k31_0.9995_config.json)

run_yacht_on_sketch ()
{
 query_sig=$1
 name=$(echo ${query_sig%_scale_1000.sig.zip})
 name=$(echo ${name##*/})
 /usr/bin/time -av -o runLog_yacht_run_${name} yacht run --json $2 --sample_file ${query_sig} --num_threads 16 --keep_raw --significance 0.95 --min_coverage_list 0.5 0.1 0.05 0.01 --out ./result_${name}.xlsx
}

export -f run_yacht_on_sketch

cat filepath_hmp_sketches.txt | parallel -j ${parallel_num} run_yacht_on_sketch {} ${json_file}

mkdir -p yacht_runlog
mv runLog* ./yacht_runlog/
mkdir -p yacht_output
mv result*.xlsx ./yacht_output

date

