#! /bin/bash

# parameters:
batchfile=$1
thread=$2  # number of CPUs in gtdb-tk
parallel_num=$3 # number of jobs to run in parallel


# activate conda environment
temp=$(which conda)
conda_path=$(echo ${temp%/*bin/conda})
if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
        . ${conda_path}/etc/profile.d/conda.sh
else
        echo "ERROR: conda path can't be corrected identified!!!"
        exit 1
fi
unset conda_path

conda activate gtdb_tk_env

# split input to subfiles of size 100
split -l 100 ${batchfile}
mkdir -p finished

# do parallel
run_parallel_gtdb_tk() {
 batch_file=$1
 out_dir=results_${batch_file}
 /usr/bin/time -av -o runlog_gtdbtk_${batch_file} gtdbtk classify_wf --batchfile ${batch_file} --out_dir ${out_dir} --skip_ani_screen --cpus ${thread}
 mv ${out_dir} ./finished/
 mv runlog_gtdbtk_${batch_file} ./finished/
 mv ${batch_file} ./finished/
}

export -f run_parallel_gtdb_tk

ls -1 x* | parallel -j ${parallel_num} run_parallel_gtdb_tk {}

date
echo "Pipe done"
