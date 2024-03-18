#! /bin/bash

# clean metadata and prepare for download

metadata="hmp_manifest_metadata_5c3ac73263.tsv"
manifest=hmp_manifest_56af75f76a.tsv
paste ${manifest} ${metadata} > merged_metadata.tsv

# remove private data, scaffold data, keep Healthy only
grep -v "Data not accessible" merged_metadata.tsv | grep -v '.scaffolds.fa.bz2' | grep Healthy > hmp_download_list.tsv

# clean file
cat <(head -1 merged_metadata.tsv) hmp_download_list.tsv > temp
sed 's/[)(]//g' temp | cut -f -5,7,10- | awk -F"\t" 'NR==1{print "f_uid""\t"$0}; NR>1{
 gsub("Healthy Human Subjects","HHS",$10);
 gsub("Inflammatory Bowel Disease Multi-omics Database IBDMDB","IBD",$10);
 gsub("prediabetes","T2D",$10);
 gsub("Human Microbiome Project HMP","HMP",$11);
 gsub("Integrative Human Microbiome Project","ihmp",$11);
print "f"NR-1"_"$11"_"$10"\t"$0}' OFS="\t"  > all_data_cleaned_metadata.tsv 
rm temp hmp_download_list.tsv merged_metadata.tsv
