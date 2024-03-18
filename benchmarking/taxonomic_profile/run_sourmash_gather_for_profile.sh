#! /bin/bash

# run sourmash gather to get taxonomic profiles on input sketches

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


# ref-db
ref_db=$(realpath ./gtdb-rs214-k31.sbt.zip)

# query data
input_filepaths=$(realpath ./filepath_hmp_sketches.txt)


run_sourmash_gather() {
 query_file=$1
 ref_db=$2
 echo ${query_file}

 # pick name
 filename=$(echo ${query_file##*/})
 out_name=$(echo ${filename} | cut -d"_" -f 1-3)

 # run gather with time command
 /usr/bin/time -av -o runlog_gather_k31_s1000_DNA_${out_name}.txt \
  sourmash gather -k 31 -o sourmash_gather_out_s1000_k31_${out_name}.csv \
  --threshold-bp 2000  ${query_file} ${ref_db}
 unset query_file ref_db filename out_name
}

export -f run_sourmash_gather

cat ${input_filepaths} | parallel -j ${parallel_num} run_sourmash_gather {} ${ref_db} 

mkdir -p single_gather_out
mkdir -p runlog_single_gather
mv runlog_*.txt ./runlog_single_gather
mv sourmash_gather_out*.csv ./single_gather_out

date
echo "done"




