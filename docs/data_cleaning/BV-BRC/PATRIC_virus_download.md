## Virus genome from BVBRC

Followup of PATRIC database, location:

```
/data/shared_data/metagenomics_graph/pathogen_database/BV-BRC/raw_data/Virus_genomes
```





Download metadata

```
wget -O genome_metadata.txt   ftp://ftp.bvbrc.org/viruses/genome_metadata
wget -O genome_summary.txt ftp://ftp.bvbrc.org/viruses/genome_summary

#
head -1 genome_summary.txt | awk -F '\t' 'BEGIN {OFS="\t"} {print $1, $3, $4, $5}' > cleaned_genome_summary.txt
awk -F '\t' 'BEGIN {OFS="\t"} {print $1, $3, $4, $5}'  genome_summary.txt >> cleaned_genome_summary.txt

wc -l cleaned_genome_summary.txt
cut -f 2 cleaned_genome_summary.txt | sed '1d' | sort | uniq -c | wc -l
cut -f 1 cleaned_genome_summary.txt | sed '1d' | sort | uniq -c | wc -l

```

There are 6M strain-level records and 51.7k species-level records.



Merge metadata

```
head -1 genome_metadata.txt |  awk -F '\t' 'BEGIN {OFS="\t"} {print $1, $4, $6, $16, $18, $19, $46, $64}'  > cleaned_metadata.txt
awk -F '\t' 'BEGIN {OFS="\t"} {print $1, $4, $6, $16, $18, $19, $46, $64}' genome_metadata.txt >> cleaned_metadata.txt

# merge metadate
# by md5sum check: the 2 files have exact same id colume, so can be paste togethe directly
md5sum <(cut -f 1 cleaned_metadata.txt)
md5sum <(cut -f 1 cleaned_genome_summary.txt)
paste cleaned_genome_summary.txt cleaned_metadata.txt | awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $3, $4, $7, $8, $9, $10, $11, $12}'  > output_Viruses_cleaned_metadata.tsv

# only keep human-related disease records for now
grep 'Human\|HUMAN\|human' output_Viruses_cleaned_metadata.tsv > output_subset_human_related_cleaned_metadata.tsv
cat <(head -1 output_Viruses_cleaned_metadata.tsv) output_subset_human_related_cleaned_metadata.tsv > _temp_aa && mv _temp_aa output_subset_human_related_cleaned_metadata.tsv

### didn't use disease filter (keep all) because most of them have NO disease annotation
```



Download data

```
mkdir -p genomes
cut -f 1 output_subset_human_related_cleaned_metadata.tsv | sed '1d' | sort -u > download_ids.txt
wc -l download_ids.txt
cut -d"." -f 1 download_ids.txt | sort -u | wc -l 
awk -F"\t" '$6' output_subset_human_related_cleaned_metadata.tsv | sed '1d' | cut -f 1 > download_assembly_avail_list.txt

awk -F"\t" '$6' output_subset_human_related_cleaned_metadata.tsv | sort -u | wc -l
awk -F"\t" '$6' output_subset_human_related_cleaned_metadata.tsv | cut -d"." -f 1 | sort -u | wc -l

################
# For human-related virus:
# There are 3M strain-level records and 11k species.
################


# seems doesn't support id-based download for virus, virus seqs are merged into higher-taxa level

# download all 
mkdir -p genomes
cd genomes
wget ftp://ftp.bvbrc.org/viruses/*.fna
cat *.fna > merged_virus_seq.fna
grep "^>" merged_virus_seq.fna | rev | cut -d"]" -f 2 | cut -d"|" -f 1 | rev | sed 's/ //g' | sort -u > merged_virus_list.txt

grep -f ../download_ids.txt merged_virus_list.txt | wc -l

# Only 65890 / 165580 ids can be found from these fna files, need find other sources.


wc -l merged_virus_list.txt
cut -d"." -f 1 merged_virus_list.txt | sort -u | wc -l 
################
# For what we've downloaded
# There are 2037k strain-level records and 155k species.
################

```

Note:

1. lots of strain-level records have NO seqeuence available OR only partial sequence available

   
