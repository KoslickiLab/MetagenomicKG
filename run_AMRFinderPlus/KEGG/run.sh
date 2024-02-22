#! /bin/bash
# Download and Run AMRFinderPlus

# set up output directory
if [ $# -eq 0 ]; then
    echo "Usage: run.sh <cpu_num>"
    exit 1
fi
current_dir=$(pwd)


function download_seq_and_run_amrfinder() {

    # download protein fasta and gff files
    if [ ! -d $1/seqs/$2 ]; then
        mkdir -p $1/seqs/$2;
    fi
    cd $1/seqs/$2;
    accession=$2;
    datasets download genome accession ${accession} --include genome,protein,gff3;
    unzip ncbi_dataset.zip;
    if [ -d ncbi_dataset/data/${accession} ]; then
        mv ncbi_dataset/data/${accession}/* .;
        # check if all required files exist
        if [ ! -f "${accession}"*_genomic.fna ] && [ ! -f protein.faa ] && [ ! -f genomic.gff ]; then
            # if genome file exits but protein or gff file does not exist, then run amrfinder without protein and gff files
            if [ -f ${accession}_*_genomic.fna ]; then
                amrfinder -n ${accession}_*_genomic.fna --plus --threads 20 > amrfinder_results.txt 2> run.log;
                # check if the output file has title only
                if [ `wc -l amrfinder_results.txt | awk '{print $1}'` -le 1 ]; then
                    echo ${accession} >> $1/empty_result_list.txt;
                    rm -rf run.log amrfinder_results.txt;
                else
                    echo ${accession} >> $1/finished_accession_list.txt;
                fi
            else
                echo ${accession} >> $1/empty_result_list.txt;
            fi
        else
            amrfinder -n ${accession}_*_genomic.fna -p protein.faa -g genomic.gff --plus --threads 20 > amrfinder_results.txt 2> run.log;
            # check if the output file has title only
            if [ `wc -l amrfinder_results.txt | awk '{print $1}'` -le 1 ]; then
                echo ${accession} >> $1/empty_result_list.txt;
                rm -rf run.log amrfinder_results.txt;
            else
                echo ${accession} >> $1/finished_accession_list.txt;
            fi
        fi
    fi
    rm -rf README.md ncbi_dataset.zip ncbi_dataset;
}
export -f download_seq_and_run_amrfinder

# create a sequence directory if not exist
if [ ! -d $current_dir/seqs ]; then
    mkdir $current_dir/seqs;
fi

## download data in parallel
less $current_dir/KEGG_microbe_assembly_mapping.tsv | sed '1d' | cut -f 2 | parallel -j $1 --link download_seq_and_run_amrfinder $current_dir {};
