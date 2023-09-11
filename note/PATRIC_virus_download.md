## Virus genome from BVBRC

Followup of PATRIC database



Download metadata

```
wget -O genome_metadata.txt   ftp://ftp.bvbrc.org/viruses/genome_metadata
wget -O genome_summary.txt ftp://ftp.bvbrc.org/viruses/genome_summary

# keep only complete or WGS data
head -1 genome_summary.txt | awk -F '\t' 'BEGIN {OFS="\t"} {print $1, $3, $4, $5}' > cleaned_genome_summary.txt
awk -F '\t' 'BEGIN {OFS="\t"} {print $1, $3, $4, $5}'  genome_summary.txt | awk '$3=="Complete"' >> cleaned_genome_summary.txt

wc -l cleaned_genome_summary.txt
cut -f 2 cleaned_genome_summary.txt | sed '1d' | sort | uniq -c | wc -l
cut -f 1 cleaned_genome_summary.txt | sed '1d' | sort | uniq -c | wc -l

```

There are 6M strain-level records and 51.7k species-level records.



Merge metadata

```
head -1 genome_metadata.txt |  awk -F '\t' 'BEGIN {OFS="\t"} {print $1, $4, $6, $16, $18, $19, $46, $64}'  > cleaned_metadata.txt
awk -F "\t" '$5=="Complete"' genome_metadata.txt |  awk -F '\t' 'BEGIN {OFS="\t"} {print $1, $4, $6, $16, $18, $19, $46, $64}' >> cleaned_metadata.txt

# merge metadate
# by md5sum check: the 2 files have exact same id colume, so can be paste togethe directly
md5sum <(cut -f 1 cleaned_metadata.txt)
md5sum <(cut -f 1 cleaned_genome_summary.txt)
paste cleaned_genome_summary.txt cleaned_metadata.txt | awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $3, $4, $7, $8, $9, $10, $11, $12}'  > output_Viruses_cleaned_metadata.tsv

# only keep human-related disease records for now
grep 'Human\|HUMAN\|human' output_Viruses_cleaned_metadata.tsv > output_subset_human_related_cleaned_metadata.tsv
cat <(head -1 output_Viruses_cleaned_metadata.tsv) output_subset_human_related_cleaned_metadata.tsv > _temp_aa && mv _temp_aa output_subset_human_related_cleaned_metadata.tsv


### didn't use disease filter (keep all) because most of them have NO disease annotation
sed '1d' output_subset_human_related_cleaned_metadata.tsv | cut -f 1 | wc -l
sed '1d' output_subset_human_related_cleaned_metadata.tsv | cut -f 1 | cut -d"." -f 1 | sort | uniq -u | wc -l
```

There are 210k strain-level records and 85 species.



Download data

```
mkdir -p genomes
cut -f 1 output_subset_human_related_cleaned_metadata.tsv | sed '1d' > download_ids.txt

# seems doesn't support id-based download for virus, virus seqs are merged into higher-taxa level

# download all 
wget ftp://ftp.bvbrc.org/viruses/*.fna
cat *.fna > merged_virus_seq.fna
grep -f ../download_ids.txt merged_virus_list.txt | wc -l

# Only 65890 / 165580 ids can be found from these fna files, need find other sources.


wget ftp://ftp.bvbrc.org/RELEASE_NOTES/viruses/genome_metadata
```




