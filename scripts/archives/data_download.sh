
#!/bin/bash

## set up working path
work_folder=$(pwd)

## set up directory
if [ ! -d "${work_folder}/data/Pathogen_data" ]
then
    mkdir ${work_folder}/data/Pathogen_data
fi

# ## Download PATRIC data (https://www.patricbrc.org/)
# if [ ! -d "${work_folder}/data/Pathogen_data/PATRIC" ]
# then
#     mkdir ${work_folder}/data/Pathogen_data/PATRIC
# fi
# loc=${work_folder}/data/Pathogen_data/PATRIC
# cd ${loc}
### the files downloaded from the links below are not completed, so I have to download data directly from https://www.patricbrc.org/ by diseases
# wget ftp://ftp.patricbrc.org/RELEASE_NOTES/genome_metadata
# wget ftp://ftp.patricbrc.org/RELEASE_NOTES/genome_lineage
# p3-all-genomes --attr genome_id,genome_name,taxon_id,disease > all_genomes_with_meta.txt

# ## Download MicroPhenoDB data (http://lilab2.sysu.edu.cn/microphenodb/#/home)
# if [ ! -d "${work_folder}/data/Pathogen_data/MicroPhenoDB" ]
# then
#     mkdir ${work_folder}/data/Pathogen_data/MicroPhenoDB
# fi
# loc=${work_folder}/data/Pathogen_data/MicroPhenoDB
# cd ${loc}
# ### this database only allows us to download data from its website


# ## Download globalrph data (https://globalrph.com/bacteria/)
# if [ ! -d "${work_folder}/data/Pathogen_data/globalrph" ]
# then
#     mkdir ${work_folder}/data/Pathogen_data/globalrph
# fi
# loc=${work_folder}/data/Pathogen_data/globalrph
# cd ${loc}
# ### this database only allows us to download data from its website

# ## Download NCBI data (https://www.ncbi.nlm.nih.gov/pathogens/organisms/)
# if [ ! -d "${work_folder}/data/Pathogen_data/NCBI" ]
# then
#     mkdir ${work_folder}/data/Pathogen_data/NCBI
# fi
# loc=${work_folder}/data/Pathogen_data/NCBI
# cd ${loc}
# ### this database only allows us to download data from its website

# ## Download VFDB data (https://www.ncbi.nlm.nih.gov/pathogens/organisms/)
# if [ ! -d "${work_folder}/data/Pathogen_data/VFDB" ]
# then
#     mkdir ${work_folder}/data/Pathogen_data/VFDB
# fi
# loc=${work_folder}/data/Pathogen_data/VFDB
# cd ${loc}
# wegt http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz

## Download CARD data from BV-BRC (ftp://ftp.bvbrc.org/)
if [ ! -d "${work_folder}/data/Pathogen_data/CARD" ]
then
    mkdir ${work_folder}/data/Pathogen_data/CARD
fi
loc=${work_folder}/data/Pathogen_data/CARD
cd ${loc}
wget ftp://ftp.bvbrc.org/specialty_genes/referenceDBs/CARD.faa