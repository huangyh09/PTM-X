#!/bin/bash
# a demo file to use PTM-X to extract features and predict on test data

### fetch features for positive and negative samples
dat_dir=$HOME/research/PTM-X/

sample_positive=$dat_dir/interface/crosstalkSamples/positive_samples.txt
sample_negative=$dat_dir/interface/crosstalkSamples/negative_samples.txt

feature_positive=$dat_dir/interface/positive_features.txt
feature_negative=$dat_dir/interface/negative_features.txt

PTMX-feature -p 20 -i $sample_positive -o $feature_positive -d $dat_dir/data/ #--verbose
PTMX-feature -p 20 -i $sample_negative -o $feature_negative -d $dat_dir/data/

## alternative option
# python PTMXtalk/feature_extractor.py -p 20 -i $sample_positive -o $feature_positive -d $dat_dir/data/ #--verbose
# python PTMXtalk/feature_extractor.py -p 20 -i $sample_negative -o $feature_negative -d $dat_dir/data/


### predict based on features
result_file=$dat_dir/interface/test_results.txt
PTMX-predict -i $feature_positive -o $result_file --positive $feature_positive --negative $feature_negative

## alternative option
# python PTMXtalk/predict.py -i $feature_positive -o $result_file --positive $feature_positive --negative $feature_negative