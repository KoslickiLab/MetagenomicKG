#! /bin/bash

# Prepare input genomes as to fit PaPrBaG format

input_train=$1 # training genome list w. 2 cols: genome_id and label
input_test=$2 # test genome list w. 2 cols: genome_id and label
filepath=$3 # txt file containg paths to all genomes
outname=$4 # out_dir name
[ -z "$outname"} ] && outname="test"


# transfer path to abs path
input_train=$(readlink -f ${input_train})
input_test=$(readlink -f ${input_test})
filepath=$(readlink -f ${filepath})


# create output dir
time_tag=$(date +%Y%m%d_%H%M)
mkdir PaPrBag_${outname}_${time_tag}
cd PaPrBag_${outname}_${time_tag}


# we need separate training and test dir because PaPrBag will read all FASTA files under a given dir
mkdir merged_training_data
mkdir merged_test_data

###### 1. merge all train data into 1 file and prepare a contig-level label
echo -n > temp_merged_training_contigs.fasta
echo -n > temp_contig_label_training.tsv

# merge fasta contigs and create labels
while IFS="\t" read -r line
do
  echo "$line"
  file_id=$(echo $line | awk '{print $1}')
  file_label=$(echo $line | awk '{print $2}')
  fasta_path=$(grep ${file_id} ${filepath})
  echo $fasta_path
  # merge fasta, may use softlink to replace merge, but will create lots of i
ntermediate rds file per sample
  cat ${fasta_path} >> temp_merged_training_contigs.fasta
  # merge label file
        grep "^>" $fasta_path | sed 's/^>//g' | awk -v label=${file_label} -v
 fid=${file_id} '{print $0"\t"label"\t"fid}' >> temp_contig_label_training.ts
v
        unset file_id file_label fasta_path
done < ${input_train}

mv temp_merged_training_contigs.fasta ./merged_training_data
mv temp_contig_label_training.tsv ./merged_training_data


###### 2. merge all validation data into 1 and make a label file, because the
 prediciton is also at contig level
echo -n > temp_merged_validation_contigs.fasta
echo -n > temp_contig_label_validation.tsv

# merge fasta contigs and create labels
while IFS="\t" read -r line
do
  echo "$line"
  file_id=$(echo $line | awk '{print $1}')
  file_label=$(echo $line | awk '{print $2}')
  fasta_path=$(grep ${file_id} ${filepath})
  echo $fasta_path
  # merge fasta, may use softlink to replace merge, but will create lots of i
ntermediate rds file per sample
  cat ${fasta_path} >> temp_merged_validation_contigs.fasta
  # merge label file
        grep "^>" $fasta_path | sed 's/^>//g' | awk -v label=${file_label} -v
 fid=${file_id} '{print $0"\t"label"\t"fid}' >> temp_contig_label_validation.
tsv
        unset file_id file_label fasta_path
done < ${input_test}

mv temp_merged_validation_contigs.fasta  ./merged_test_data
mv temp_contig_label_validation.tsv ./merged_test_data

# we need the genome label file for validation
cp ${input_test} ./merged_test_data/test_genome_label.tsv


date





