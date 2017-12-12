#!/bin/bash

PTMX_dir=$HOME/research/PTM-X
BIN_dir=`pwd`/PTMsites_src

### step 1: download PTM data (stored 08/11/2017, downloaded 06/12/2017)
# path=`pwd`
# mkdir $PTMX_dir/data/PTMsites/phosphosite
# cd $PTMX_dir/data/PTMsites/

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
# DataDir=$PTMX_dir/data/PTMsites
# for SPECIES in human mouse rat
# do
# 	OUTFILE=PTMsites_$SPECIES.txt
# 	python $BIN_dir/PTMsites_process.py --data_dir=$DataDir/phosphosite --file_list=$LIST --species=$SPECIES --out_file=$DataDir/$OUTFILE
# 	echo processed the PTM sites for $SPECIES from PhosphoSitePlus!
# done


### step 3: orthology across human, mouse and rat
# ## protein sequence and orthologous data
# mkdir $PTMX_dir/data/ortholog
# mkdir $PTMX_dir/data/ortholog/InParanoid.raw
# cd $PTMX_dir/data/ortholog/InParanoid.raw
# wget http://inparanoid.sbc.su.se/download/current/sequences/processed/9606.fasta
# wget http://inparanoid.sbc.su.se/download/current/sequences/processed/10090.fasta
# wget http://inparanoid.sbc.su.se/download/current/sequences/processed/10116.fasta

# wget http://inparanoid.sbc.su.se/download/current/Orthologs_other_formats/H.sapiens/InParanoid.H.sapiens-M.musculus.tgz
# wget http://inparanoid.sbc.su.se/download/current/Orthologs_other_formats/H.sapiens/InParanoid.H.sapiens-R.norvegicus.tgz

# tar -xvzf InParanoid.H.sapiens*tgz


# ## step 4: combine orthology
# python PTMsites_src/inParanoid_merge.py
# cd $PTMX_dir/data/ortholog
# mkdir InParanoid.fasta InParanoid.align
# python $BIN_dir/ortholog2fasta.py

# ## MSA tool: MUSCLE
# mkdir $PTMX_dir/data/ortholog/tools
# cd $PTMX_dir/data/ortholog/tools
# wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz

# tar -xvzf muscle3.8.31_i86linux64.tar.gz

# ## step 5: align with MUSCLE
# FASTAdir=$PTMX_dir/data/ortholog/InParanoid.fasta
# ALIGNdir=$PTMX_dir/data/ortholog/InParanoid.align

# MUSCLE=$PTMX_dir/data/ortholog/tools/muscle3.8.31_i86linux64

# files=`ls $FASTAdir`

# for f in $files
# do
# 	fname=`echo "$f" |awk '{split($1,a,"."); print a[1]}'`
# 	$MUSCLE -in $FASTAdir/$f -out $ALIGNdir/$fname.afa
# done


## step 6: merge PTM data for human, mouse and rat
python $BIN_dir/PTM_species_merge.py -d $PTMX_dir/data/




