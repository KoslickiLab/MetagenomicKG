## PATRIC information:

Shaopeng Liu (sml6467@psu.edu)





Conclusion: PATRIC or BV-BRC contains strain level records that link microbes/viruses to disease and AMR. The corresponding genomes can be directly retrived from NCBI by accession id OR taxon id + strain name. So it can be directly augmented to KEGG at strain-level. However, we still need phenotype data.



1. Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3965095/

2. 1. information: bacteria only, but integrated into a larger on called BV-BRC with bac + virus
   2. multiple data source: nt, transcriptome, protein-protein interaction, 3d structure
   3. taxonomy: NCBI taxonomy, updated every month

3. Document: [https://www.bv-brc.org/](https://docs.patricbrc.org/user_guides/ftp.html)

4. 1. Data API: https://www.bv-brc.org/api/doc/
   2. FTP DB: https://www.bv-brc.org/docs/quick_references/ftp.html
   3. CLI (sudo needed): https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html#installation-on-debian-ubuntu-mint-linux

5. Metadata exploration: GPU server /data/shared_data/PATRIC

6. 1. (Only downloaded metadata for now due to shortage of space)
   2. notes:
      1. WGS (collection of contigs from WGS data, but MAY not be complete) vs complete: https://www.biostars.org/p/380919/
      2. genome id: aa.bb where aa is taxon id and bb is strain identifier



==PATRIC DB contains genomes that are not in NCBI, and also additional data==

```
https://www.bv-brc.org/docs/quick_references/ftp.html

for i in `cat genome_list`; do wget -qN "ftp://ftp.bvbrc.org/genomes/$i/$i.fna";
done

# or use wget to download the folder
wget -r --no-parent ftp://ftp.bvbrc.org/genomes/$i/
```

Other data (the tutorial missed some "PATRIC" word in suffix)

| file                    | description                                                  | note                                                    |
| ----------------------- | ------------------------------------------------------------ | ------------------------------------------------------- |
| ==.fna==                | contig-level DNA sequence                                    |                                                         |
| ==.PATRIC.faa==         | gene-level protein sequence                                  |                                                         |
| .PATRIC.features.tab    | gene annotation with figfam_id, plfam_id, pgfam_id           | No gene name                                            |
| .PATRIC.ffn             | gene-level DNA sequence                                      | Match to PATRIC.faa file                                |
| .PATRIC.frn             | FASTA nucleotide sequences for RNAs                          |                                                         |
| .PATRIC.gff             | Genome annotations in GFF file format                        | feature file is better                                  |
| ==.PATRIC.pathway.tab== | Metabolic pathway assignments (KEGG)                         | use it!                                                 |
| ==.PATRIC.spgene.tab==  | Specialty gene assignements (i.e. AMR genes, virulance factors, essential genes, etc) | for AMR, VRF                                            |
| .PATRIC.subsystem.tab   | Subsystem assignments                                        | description of gene functions, not sure how to use them |





---

### FTP Data exploration:

Location GPU server: /data/shared_data/PATRIC



Download metadata

```
wget -O genome_summary.txt  ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_summary
wget -O genome_metadata.txt ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_metadata
wget -O genome_lineage.txt ftp://ftp.bvbrc.org/RELEASE_NOTES/genome_lineage
wget ftp://ftp.bvbrc.org/RELEASE_NOTES/PATRIC_genomes_AMR.txt
```



#### Genome status for WGS or Complete only

```
head -1 genome_summary.txt | cut -f 1,3,4,5,17,18 > cleaned_genome_summary.txt
cut -f 1,3,4,5,17,18  genome_summary.txt | awk '$3=="Complete" || $3=="WGS"' >> cleaned_genome_summary.txt

wc -l cleaned_genome_summary.txt
cut -f 2 cleaned_genome_summary.txt | sed '1d' | sort | uniq -c | wc -l
cut -f 1 cleaned_genome_summary.txt | sed '1d' | sort | uniq -c | wc -l

```



Though many records have no accession ID, **you can find them in PATRIC's own database**, but I'm not sure about the source of them.

| BV-BRC records | Number of genomes | BioProj id / accession id | With disease information | BioProj id / accession id |
| -------------- | ----------------- | ------------------------- | ------------------------ | ------------------------- |
| Species-level  | 90k               |                           | 3k                       |                           |
| Strain-level   | 633k              | 573k / 350k               | 18k                      | 6k, 3.7k                  |

Summary:

1. BV(PATRIC) database records strain-level genomes and metadata, which is great
2. the DB uses NCBI taxonomy and contains bacteria and viral genomes
3. data identifier: 
   1. the combination of taxon id and strain name is an unique identifier to retrive genomes from NCBI
   2. col18 is assembly_accession, but contains both RefSeq and GenBank (use 1 if we only want 1 resource)
   3. the data id itself is unique for PATRIC DB and can be used to download genomes from its FTP server



#### Metadata

```
head -1 genome_metadata.txt |  cut -f 1,4,6,16,18,46,64   > cleaned_metadata.txt
awk -F "\t" '$5=="Complete" || $5=="WGS"' genome_metadata.txt |  cut -f 1,4,6,16,18,46,64 >> cleaned_metadata.txt

# merge metadate
# by md5sum check: the 2 files have exact same id colume, so can be paste togethe directly
md5sum <(cut -f 1 cleaned_metadata.txt)
md5sum <(cut -f 1 cleaned_genome_summary.txt)
paste cleaned_genome_summary.txt cleaned_metadata.txt | cut -f 1,3-6,9-  > output_PATRIC_cleaned_metadata.tsv


# only keep human-related disease records for now
grep Human output_PATRIC_cleaned_metadata.tsv > output_subset_human_related_cleaned_metadata.tsv
cat <(head -1 output_PATRIC_cleaned_metadata.tsv) output_subset_human_related_cleaned_metadata.tsv > _temp_aa && mv _temp_aa output_subset_human_related_cleaned_metadata.tsv

awk -F"\t" '$10' output_subset_human_related_cleaned_metadata.tsv > output_subset_human_disease_cleaned_metadata.tsv





### final metadata status
# disease number / strain number
wc -l output_subset_human_disease_cleaned_metadata.tsv  # 15519

# species number:
cut -f 1 output_subset_human_disease_cleaned_metadata.tsv | sed '1d' | cut -d"." -f 1 | sort | uniq -c | wc -l  #1549

# disease distribution
sed '1d' output_subset_human_disease_cleaned_metadata.tsv | cut -f 10  | sort | uniq -u > disease_count.txt

```

Summary:

1. used to retrive genomes based on id

2. provide additional metadata, but couldn't find a universal standard to use

3. **didn't find phenotype data**

4. some cols that might be useful (majority are empty)

   1. col52: patient infor
   2. col53: AMR status (Resistant/intermediate/suspective/etc.)
   3. col54: AMR evidence (putative or AMR panel)

5. ==How about those records with NO accession ID? Examples:==

   1. Some have BioProject id but No accession ID (missing records): https://www.ncbi.nlm.nih.gov/biosample/SAMN29593090/

   2. Some have Accession ID but No BioProject id (dup records):

      ```
      grep -w 40214.31 genome_metadata.txt
      grep -w 40214.23 genome_metadata.txt
      
      
      40214.31	Acinetobacter johnsonii strain Aj2199		40214	WGS	Aj2199				2016-04-26T00:00:00Z				GCA_001632325.1	LVIB00000000					154		3799071	41.39	3917				peritoneal fluid		2014	Argentina	Argentina: Buenos Aires						Human, Homo sapiens
      
      40214.23	Acinetobacter johnsonii strain Aj2199		40214	WGS	Aj2199				2016-04-26T00:00:00Z		PRJNA315996	SAMN04572947	GCF_001632325.1	LVIB00000000		CSUFIllumina HiSeq	65.0x	spades v. 3.1			154		3799071	41.39	3919			peritoneal fluid		2014	Argentina	Argentina: Buenos Aires						Human, Homo sapiens			infection								Acinetobacter johnsoni clinical strain (Aj2199) that was co-producing PER-2 and OXA-58	collected_by:Marisa Almuzara
      ```

   3. Some have neither accession id nor bioproject id: hard to locate data

      Example:

      ```
      562.96458 -> Escherichia coli VREC0659
      
      1773.29523 -> Mycobacterium tuberculosis 6323-01
      
      485.15075 ->	Neisseria gonorrhoeae AT13_495
      ```




6. Final data status for human-related records

   | Item           | Record |
   | -------------- | ------ |
   | Strain number  | 2018   |
   | Species number | 1539   |
   | Disease count  | 156    |

   





#### Genome lineage

Gives a full lineage information for each strain-level record.



#### PATRIC_genomes_AMR.txt

Records antibiotic information for strain-level records.





### Download all genome data and get FMH similarity

---

Build sketches

```bash
mkdir -p PATRIC_genomes 
cd PATRIC_genomes

# download all related information
cut -f 1 ../output_subset_human_disease_cleaned_metadata.tsv | sed '1d' > download_ids.txt

# use nohup
date > record.txt

for file in $(cat download_ids.txt); do
 wget -r --no-parent ftp://ftp.bvbrc.org/genomes/${file}/ || echo -e "${file}\tFailed......" >> record.txt
done




# build sourmash sketch
mkdir -p sourmash_sketch
cd sourmash_sketch

readlink -f ../PATRIC_genomes/*.fna > patric_genome_list.txt

for scaleFactor in 100 1000; do
 # run protein sketch
 /usr/bin/time -av -o runlog_scale_${scaleFactor}_patric_protein.txt \
  sourmash sketch translate -f \
  -p protein,k=5,k=7,k=11,k=15,k=25,abund,scaled=${scaleFactor} \
  -o PATRIC_protein_scale_${scaleFactor}.sig.zip \
  --from-file patric_genome_list.txt
  
 # run dna sketch
 /usr/bin/time -av -o runlog_scale_${scaleFactor}_patric_dna.txt \
  sourmash sketch dna -f \
  -p k=21,k=31,k=51,abund,scaled=${scaleFactor} \
  -o PATRIC_dna_scale_${scaleFactor}.sig.zip \
  --from-file patric_genome_list.txt
done






# sourmash compare
# PATRIC * KEGG gene matrix is too large (requrie 1.4PB MEM to process by command line tools)
kegg_ko_faa_db=/data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/KOs_sketched_scaled_100.sig.zip

for kSize in 5,7,11,15; do
 sourmash compare --protein -k ${kSize} -p 50 --ignore-abundance --containment --csv compare_protein_k_${kSize}_with_KEGG_KOs.csv sourmash_protein_scale_100.sig.zip ${kegg_ko_faa_db}
done

```



Compare sketches manually if db too large

```
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
conda activate sourmash


# file locations
py_code="/data/shared_data/PATRIC/script/py_wrapper_sourmash_for_ciji.py"
DATADIR=/data/shared_data/PATRIC/sourmash_sketch
cd ${DATADIR}
kegg_dna_db=/data/shared_data/KEGG_data/sourmash_sketches/output/kegg_genes_KO.fna_scale_10.db.zip
patric_dna_db=/data/shared_data/PATRIC/sourmash_sketch/sourmash_dna_scale_10.db.zip



/usr/bin/time -av -o runlog_compare_dna_PATRIC_KEGG.txt python ${py_code} -r ${kegg_dna_db} -q ${patric_dna_db} -k 21,31,51 -m CI -o PATRIC_KEGG_CI_scale10. 

echo "pipe done"
date
```







Note: Sourmash containment format

The containment matrix is organized such that the value in row A for column B is the containment of the B’th sketch in the A’th sketch, i.e.

```
C(A, B) = B.contained_by(A)
```

This is NOT containment index, it's reversed: proportion of B in A. 

==This is Containment index of column B in row A==

```
a1.fna                  a2.fna                  6%_a1.fna
1.0                     0.00026920206061940946  1.0
0.00037823619130891825  1.0                     0.0012168002349683211
0.058535734723801615    5.069389453222646e-05   1.0
```




