#!/bin/bash

# dat_dir=$HOME/research/PTM/PTM-X_V2.1/
dat_dir=$HOME/research/PTM-X/

sample_file=$dat_dir/interface/crosstalkSamples/positive_samples.txt
feature_file=$dat_dir/interface/positive_features_all.txt

python feature_extractor.py -i $sample_file -o $feature_file -d $dat_dir/data/


sample_file=$dat_dir/interface/crosstalkSamples/negative_samples.txt
feature_file=$dat_dir/interface/negative_features_all.txt

python feature_extractor.py -i $sample_file -o $feature_file -d $dat_dir/data/