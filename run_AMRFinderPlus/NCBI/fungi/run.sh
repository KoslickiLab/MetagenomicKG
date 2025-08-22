#! /bin/bash
# Download Genomes/Nucleotide Seqs and Run AMRFinderPlus

# set up output directory
if [ $# -eq 0 ]; then
    echo "Usage: run.sh <cpu_num>"
    exit 1
fi
current_dir=$(pwd)


## define a bash function
function download_seq_and_run_amrfinder() {
    # go to the seq directory
    cd $1/fungi_seqs;
    # create an individual folder for each taxid
    if [ ! -d $1/fungi_seqs/$2 ]; then
        mkdir -p $1/fungi_seqs/$2;
    fi

    # go to the taxid folder
    cd $1/fungi_seqs/$2;

    # download data
    ncbi-genome-download -s refseq -F fasta,gff,protein-fasta -l all -t $2 fungi;
    if [ ! -d $1/fungi_seqs/$2/refseq ]; then
        ncbi-genome-download -s genbank -F fasta,gff,protein-fasta -l all -t $2 fungi;
        if [ ! -d $1/fungi_seqs/$2/genbank ]; then
            echo $2 >> $1/empty_result_list.txt;
            return;
        fi
        filename=`find ./genbank -type f -name "*.fna.gz" | head -1`
    else
        filename=`find ./refseq -type f -name "*.fna.gz" | head -1`
    fi
    if [ -n "$filename" ]; then
        dir_path=`dirname $filename`;
        amrfinder -n ${dir_path}/*.fna.gz -p ${dir_path}/*.faa.gz -g ${dir_path}/*.gff.gz --plus --threads 20 > amrfinder_results.txt 2> run.log;
        # check if the output file has title only (check if the line number is smaller or equal to 1)
        if [ `wc -l amrfinder_results.txt | awk '{print $1}'` -le 1 ]; then
            echo $2 >> $1/empty_result_list.txt;
            rm -rf run.log amrfinder_results.txt;
        else
            echo $2 >> $1/finished_taxids_list.txt;
        fi
    else
        echo $2 >> $1/empty_result_list.txt;
    fi
    rm -rf refseq genbank;
}
export -f download_seq_and_run_amrfinder

## Create a folder to store the seq data
if [ ! -d $current_dir/fungi_seqs ]; then
    mkdir -p $current_dir/fungi_seqs;
fi

## Find the unique taxi ids that we need to download its genomes/nucleotide seqs
cd $current_dir;
less fungi_hierarchy.tsv | cut -f 3 | sed '1d' | sort -u > parent_taxids.txt;
less fungi_hierarchy.tsv | cut -f 6 | sed '1d' | sort -u > child_taxids.txt;
comm -13 <(sort  parent_taxids.txt) <(sort child_taxids.txt) > unique_download_taxids.txt


## download data in parallel
# Create a temporary file with filtered IDs that haven't been processed yet
if [ -f $current_dir/empty_result_list.txt ] && [ -f $current_dir/finished_taxids_list.txt ]; then
    echo "Filtering out already processed genomes..."
    cat $current_dir/unique_download_taxids.txt | grep -v -f $current_dir/empty_result_list.txt | grep -v -f $current_dir/finished_taxids_list.txt > $current_dir/unique_download_taxids.txt.filtered
    echo "Original genome count: $(wc -l < $current_dir/unique_download_taxids.txt)"
    echo "Filtered genome count: $(wc -l < $current_dir/unique_download_taxids.txt.filtered)"
    
    # Use the filtered list for processing
    cat $current_dir/unique_download_taxids.txt.filtered | parallel -j $1 --link download_seq_and_run_amrfinder $current_dir {};
    
    # Clean up temporary file
    rm -f $current_dir/unique_download_taxids.txt.filtered
else
    echo "'empty_result_list.txt' and 'finished_taxids_list.txt' do not exist, processing all genomes."
    less $current_dir/unique_download_taxids.txt | parallel -j $1 --link download_seq_and_run_amrfinder $current_dir {};
fi


