#! /bin/bash

# Prepare npy data to fit DeePac input


# need 5 positional parameters:
# 1. training file with label
# 2. validation file with label
# 3. keyword for folder naming
# 4. filepath for all genomes
# 5. py wrappter to merge npy

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

file_train=$1
file_valid=$2
folder_keyword=$3
filepath=$4
prepare_metric_from_frag=$5


###### abs path
file_train=$(realpath ${file_train})
file_valid=$(realpath ${file_valid})
filepath=$(realpath ${filepath})
prepare_metric_from_frag=$(realpath ${prepare_metric_from_frag})

###### create folder
time_tag=$(date +%Y%m%d_%H%M)
mkdir DeePac_${folder_keyword}_${time_tag}
cd DeePac_${folder_keyword}_${time_tag}


###### merge pos / neg fastas together for easier process
merge_pos_neg_fasta () {
 input_file=$1
 out_name=$2

 # merge pos data
 echo -n > merged_pos_${out_name}.fasta
 for f_id in $(awk '$2=="True" {print $1}' $input_file); do
  echo -e "${f_id}\tTrue"
  fasta_file=$(grep ${f_id} ${filepath})
  cat ${fasta_file} >> merged_pos_${out_name}.fasta
 done

 # merge neg data
 echo -n > merged_neg_${out_name}.fasta
 for f_id in $(awk '$2=="False" {print $1}' $input_file); do
  echo -e "${f_id}\tFalse"
  fasta_file=$(grep ${f_id} ${filepath})
  cat ${fasta_file} >> merged_neg_${out_name}.fasta
 done
}
export -f merge_pos_neg_fasta


###### generate merged fasta
mkdir merged_fasta
cd merged_fasta
merge_pos_neg_fasta  ${file_train}  training
merge_pos_neg_fasta  ${file_valid} valid
cd ..


###### generate npy matric by deepac
deepac gwpa fragment -g merged_fasta -o input_npy -s 250

###### generate label vector
cd input_npy
rm *.fasta
python ${prepare_metric_from_frag}
rm merged*.npy

###### ready to train the model!
date
echo "Done"



