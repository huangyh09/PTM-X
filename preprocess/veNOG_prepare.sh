#!/bin/bash

data_dir=$HOME/research/PTM-X/data

### Part1: download veNOG data
path=`pwd`
# mkdir $data_dir/eggNOG4/
# cd $data_dir/eggNOG4/
# wget http://eggnogdb.embl.de/download/eggnog_4.5/data/veNOG/veNOG.raw_algs.tar.gz
# tar -xvzf veNOG.raw_algs.tar.gz
# cd $path
# echo downloaded veNOG data!

### Part 2: Download Uniprot and Emsemble id maps from BioMart
## - go to bioMart, choose Human genes (GRCh38.p10), 
## - choose feature: UniProtKB/Swiss-Prot ID and Protein stable ID (Unique results only)
## - link: https://www.ensembl.org/biomart/martview/adcdf2493e514f1915b479493df03831


### Part3: mapping the veNOG id to Ensembl protein id by searching from 
### the veNOG file, and then map to Uniprot id

OUT_FILE=$data_dir/eggNOG4/veNOG_Ensembl_Uniprot.txt
mart_FILE=$data_dir/eggNOG4/mart_export.txt
veNOG_dir=$data_dir/eggNOG4/veNOG_raw_algs/
python map_veNOG_Uniprot.py -o $OUT_FILE -d $veNOG_dir -m $mart_FILE




