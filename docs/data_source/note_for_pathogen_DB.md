## Integrate pathogen genomes with KEGG

#### Pathogen databases

1. BV-BRC
1. PHI
1. MicroPhenoDB



File information

```
Where I got this:
1. MicroPhenoDB: http://www.liwzlab.cn/microphenodb/#/download
2. PHI: http://www.phi-base.org/index.jsp
3. BV-BRC: https://www.bv-brc.org/docs/quick_references/ftp.html



MicroPhenoDB:
1. core_table.txt: association table between microbe species name (col2) and disease name (col3)
2. EFO.txt: disease annotation file for name <-> Experimental Factor Ontology (EFO)
3. GeneID.txt (not necessary): core genes found from each genome with accession id and gene id
4. NCIT.txt: gives NCI annotation of species, can be used to retrieve taxid

Key genes are "computed" from marker gene list. So it's equivalent for us to run a sourmash gather.



PHI:
1. pathogen_species_list.csv: taxid and pathogen name
2. diseases.csv: list of disease names
3. phi-base_current.csv: core table, check here for annotation: http://www.phi-base.org/helpLink.htm. This is a gene-level mapping file.
    1. PathogenID: taxid
    2. Pathogenspecies: name
    3. Pathogenstrain: strain information (if any)
    4. disease: disease name
    5. Hostdescription
    6. GeneFunction
    7. GOannotation
    8. Pathway: if any


BV-BRC:
in raw_data folder
1. output_subset_human_disease_cleaned_metadata.tsv: genome_id, host, and disease
2. PATRIC_genomes_AMR.txt: genome_id with AMR
3. PATRIC_genomes: ALL genomic features of selected genome_id, including gene-level seq, pathway, etc.

```



