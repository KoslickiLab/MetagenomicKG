#! /bin/bash
# Download Genomes/Nucleotide Seqs and Run AMRFinderPlus

# set up output directory
if [ $# -eq 0 ]; then
    echo "Usage: run.sh <cpu_num>"
    exit 1
fi
current_dir=$(pwd)


## define a bash function
function run_amrfinder() {
    # go to the seq directory
    cd $1/$2;

    # find the genome id
    genome_id=`echo $2 | awk -F '/' '{print $2}'`;

    if [ -f amrfinder_results.txt ]; then
        rm -rf amrfinder_results.txt;
    fi
    if [ -f run.log ]; then
        rm -rf run.log;
    fi
    if [ -f ${genome_id}.PATRIC_temp.gff ]; then
        rm -rf ${genome_id}.PATRIC_temp.gff;
    fi
    
    # check if .fna, .faa, .gff files exist
    if [ -f ${genome_id}.fna ] && [ -f ${genome_id}.PATRIC.faa ] && [ -f ${genome_id}.PATRIC.gff ]; then
        awk 'BEGIN {FS=OFS="\t"} ($2 == "PATRIC") {sub(/^[^|]*\|/, "", $1); match($9, /ID=([^;]+)/, arr); $9 = $9 ";Name=" arr[1]} 1' ${genome_id}.PATRIC.gff > ${genome_id}.PATRIC_temp.gff
        amrfinder -n ${genome_id}.fna -p ${genome_id}.PATRIC.faa -g ${genome_id}.PATRIC.gff --plus --threads 20 > amrfinder_results.txt 2> run.log;
        if [ `wc -l amrfinder_results.txt | awk '{print $1}'` -le 1 ]; then
            amrfinder -n ${genome_id}.fna --plus --threads 20 > amrfinder_results.txt 2> run.log;
            if [ `wc -l amrfinder_results.txt | awk '{print $1}'` -le 1 ]; then
                echo ${genome_id} >> $1/empty_result_list.txt;
                rm -rf amrfinder_results.txt;
            else
                echo ${genome_id} >> $1/finished_taxids_list.txt;
            fi
        fi
        rm -rf ${genome_id}.PATRIC_temp.gff;
        return;
    else
        if [ -f ${genome_id}.fna ]; then
            amrfinder -n ${genome_id}.fna --plus --threads 20 > amrfinder_results.txt 2> run.log;
            if [ `wc -l amrfinder_results.txt | awk '{print $1}'` -le 1 ]; then
                echo ${genome_id} >> $1/empty_result_list.txt;
                rm -rf amrfinder_results.txt;
            else
                echo ${genome_id} >> $1/finished_taxids_list.txt;
            fi
        fi
        return;
    fi
}
export -f run_amrfinder

## download data in parallel
cd $current_dir;
find genomes/ -maxdepth 1 -type d | parallel -j $1 --link run_amrfinder $current_dir {};
