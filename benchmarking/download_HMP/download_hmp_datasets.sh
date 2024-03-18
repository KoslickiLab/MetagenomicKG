# may encounter failed urls, so add this extra step to check
sed 1d ../all_data_cleaned_metadata.tsv | cut -f 1,3,5  > file_2_download.tsv
input_file=file_2_download.tsv

date
time_tag=$(date +%h%d_%H%M)
echo -e "f_id\tmd5\tlink" > error_record_${time_tag}.tsv

cat ${input_file} | while read line; do
  # echo "$line"
  file_id=$(echo ${line} | awk '{print $1}')
  md5=$(echo ${line} | awk '{print $2}')
  link=$(echo ${line} | awk '{print $3}' | cut -d"," -f 1)
  default_name=$(echo ${link##*/})

  # download
  wget -nv -O ${file_id}_${default_name} ${link}

  if [[ "$?" != 0 ]]; then
   echo "Error downloading file"
   rm ${file_id}_${default_name}
   echo -e "${line}\tFailed_download" >>  error_record_${time_tag}.tsv
  fi
done

rm file_2_download.tsv
mv error_record_${time_tag}.tsv ../failed_download_record.tsv

# check md5 for all downloaded files
echo -e "f_id\tmeta_md5\tfile_md5" > ../md5_not_match_record.tsv
for file in $(ls -1 ); do
 echo $file
 f_id=$(echo $file | cut -d"_" -f 1-3)
 meta_md5=$(grep $f_id ../all_data_cleaned_metadata.tsv | cut -f 3)
 file_md5=$(md5sum $file | awk '{print $1}')
 [ "$meta_md5" != "$file_md5" ] && echo -e "${file}\t${meta_md5}\t${file_md5}" >> ../md5_not_match_record.tsv
done

# make 2 folders
mkdir md5_not_match
sed '1d' ../md5_not_match_record.tsv | cut -f 1 | xargs mv -t ./md5_not_match

mkdir checked_fq_files
mv f* ./checked_fq_files


echo "Finish downloading"
echo "start time is ${time_tag} (hd_HM)"
echo "end time is: $(date +%h%d_%H%M)"
date
