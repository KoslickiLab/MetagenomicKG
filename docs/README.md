## Documents and development notes

This folder contains internal documents and development notes regarding data source, data cleaning, and data analysis. It's for internal development purpose and for future lab members who wants to rebuild, edit, or change the data sources. 



### Contents:

1. [Data source](#source)
2. [Data cleaning](#clean)
3. [Data processing](#process)





### Data source  <a name="source"></a>

---

We haven't assign version for MetagenomicKG (recommended), now let's call it v1.00.

1. Taxonomy
   1. GTDB: downloads from its server [here](https://data.gtdb.ecogenomic.org/releases/latest/).
   2. NCBI: downloads from its server [here](https://support.nlm.nih.gov/knowledgebase/article/KA-03474/en-us) (we only uses its virus taxonomy).
2. Functional annotation
   1. KEGG: downloads from its FTP server (licensed required)
   2. AMR and VF: downloads from NCBI reference database [here](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/). 
3. Pathogen 
   1. BV-BRC: downloads from its server [here](https://www.bv-brc.org/docs/quick_references/ftp.html). 
   2. MicroPhenoDB: downloads from its server [here](http://www.liwzlab.cn/microphenodb/#/download).
4. Existing biomedical knowledge graph
   1. RTX-KG2: check usage [here](https://github.com/RTXteam/RTX-KG2).



</br>



### Data cleaning  <a name="clean"></a>

---

| Database        | Dir       | Contents                                                     |
| --------------- | --------- | ------------------------------------------------------------ |
| BV-BRC (PATRIC) | BV-BRC    | Metadata exploration and data cleaning process               |
| VFDB and CARD   | VFDB_CARD | Exploration of database for virulence factor (VFDB) and AMR genes (CARD). We don't use them in MKG, but uses NCBI curated references instead. |



</br>



### Data process  <a name="process"></a>

---

Most of works here are finalized in the `benchmarking` folder. Keep this lot for now for future notes.

