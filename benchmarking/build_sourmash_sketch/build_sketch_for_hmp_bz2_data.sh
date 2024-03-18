#! /bin/bash

# build Sourmash sketches for HMP bz2 data. You can find more data types in the Functional_profile manuscript and GitHub

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


# process bz2 file 
process_tarbz_sourmash ()
{
 tar_bz_file=$1
 scaleFactor=1000
 filename=$(echo ${tar_bz_file##*/})
 out_name=$(echo ${filename} | cut -d"_" -f 1-3)
 tar -xvf ${tar_bz_file} > temp_file_${out_name}.txt
 # 1st line is folder
 format=$(sed -n '2p' temp_file_${out_name}.txt | rev | cut -d"." -f 1 | rev)
 # merge reads into 1 file
 for file in $(sed '1d' temp_file_${out_name}.txt); do
   cat ${file} >> merged_${out_name}.${format}
   rm ${file}
 done
 # check if folder become empty after merging
 dirname=$(sed -n '1p' temp_file_${out_name}.txt)
 rmdir ${dirname} || echo -e "${out_name}\t${dirname}" >> error_records.txt
 # run sourmash
 /usr/bin/time -av -o runlog_scale_${scaleFactor}_${filename}.txt sourmash sketch dna -p k=31,abund,scaled=${scaleFactor} -o ${out_name}_scale_${scaleFactor}.sig.zip merged_${out_name}.${format}
 # clean temp files
 echo ${out_name} >> file_note_all_bz2.txt
 cat temp_file_${out_name}.txt >> file_note_all_bz2.txt
 rm temp_file_${out_name}.txt merged_${out_name}.${format}
 unset tar_bz_file scaleFactor filename out_name dirname format
}

export -f process_tarbz_sourmash



# run in parallel
cat filepath_all_bz2.txt | parallel -j ${parallel_num} process_tarbz_sourmash {}

mkdir -p single_sketch
mv *.sig.zip ./single_sketch
mkdir -p log_file
mv runlog* ./log_file

date
echo "Pipe done"




