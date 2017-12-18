#!/bin/bash

BIN_dir=`pwd`/PTMsites_src
data_dir=$HOME/research/PTM-X/data

### step 1: download PTM data (stored 08/11/2017, downloaded 06/12/2017)
# path=`pwd`
# mkdir $data_dir/PTMsites/phosphosite
# cd $data_dir/PTMsites/

# wget http://www.phosphosite.org/downloads/Acetylation_site_dataset.gz
# wget http://www.phosphosite.org/downloads/Methylation_site_dataset.gz
# wget http://www.phosphosite.org/downloads/O-GalNAc_site_dataset.gz
# wget http://www.phosphosite.org/downloads/O-GlcNAc_site_dataset.gz
# wget http://www.phosphosite.org/downloads/Phosphorylation_site_dataset.gz
# wget http://www.phosphosite.org/downloads/Sumoylation_site_dataset.gz
# wget http://www.phosphosite.org/downloads/Ubiquitination_site_dataset.gz

# gzip -d *_dataset.gz
# cd $path
# echo downloaded PTM data from PhosphoSitePlus!


### step 2: combine multiple PTMs
# LIST=$BIN_dir/PTMsites_file.lst
# DataDir=$data_dir/PTMsites
# for SPECIES in human mouse rat
# do
#   OUTFILE=PTMsites_$SPECIES.txt
#   python $BIN_dir/PTMsites_process.py --data_dir=$DataDir/phosphosite --file_list=$LIST --species=$SPECIES --out_file=$DataDir/$OUTFILE
#   echo processed the PTM sites for $SPECIES from PhosphoSitePlus!
# done


### step 3: orthology across human, mouse and rat
# ## protein sequence and orthologous data
# mkdir $data_dir/ortholog
# mkdir $data_dir/ortholog/InParanoid.raw
# cd $data_dir/ortholog/InParanoid.raw
# wget http://inparanoid.sbc.su.se/download/current/sequences/processed/9606.fasta
# wget http://inparanoid.sbc.su.se/download/current/sequences/processed/10090.fasta
# wget http://inparanoid.sbc.su.se/download/current/sequences/processed/10116.fasta

# wget http://inparanoid.sbc.su.se/download/current/Orthologs_other_formats/H.sapiens/InParanoid.H.sapiens-M.musculus.tgz
# wget http://inparanoid.sbc.su.se/download/current/Orthologs_other_formats/H.sapiens/InParanoid.H.sapiens-R.norvegicus.tgz

# tar -xvzf InParanoid.H.sapiens*tgz


# ## step 4: combine orthology
# python PTMsites_src/inParanoid_merge.py
# cd $data_dir/ortholog
# mkdir InParanoid.fasta InParanoid.align
# python $BIN_dir/ortholog2fasta.py

# ## MSA tool: MUSCLE
# mkdir $data_dir/ortholog/tools
# cd $data_dir/ortholog/tools
# wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz

# tar -xvzf muscle3.8.31_i86linux64.tar.gz

## step 5: align with MUSCLE
# FASTAdir=$data_dir/ortholog/InParanoid.fasta
# ALIGNdir=$data_dir/ortholog/InParanoid.align
# MUSCLE=$data_dir/ortholog/tools/muscle3.8.31_i86linux64

# for f in `ls $FASTAdir`
# do
#   fname=`echo "$f" |awk '{split($1,a,"."); print a[1]}'`
#   $MUSCLE -in $FASTAdir/$f -out $ALIGNdir/$fname.fa
# done


## step 6: merge PTM data for human, mouse and rat
python $BIN_dir/PTM_species_merge.py -d $data_dir/




