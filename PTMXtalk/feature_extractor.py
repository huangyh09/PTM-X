#!/usr/bin/python2.7

# This file is to extract features for a given file of samples.
# The features include: 
# 1) co-evolution on sequence
# 2) co-conservation on modifications
# 3) co-occurrence (from Rongting)

import os
import numpy as np
from optparse import OptionParser, OptionGroup
from utils.cooccur import fetch_coOccur
from utils.coevolution import fetch_seqCoEvol, fetch_PTMcoEvol

def main():
    #0. parse command line options
    parser = OptionParser()
    
    parser.add_option("-i", "--input-file", dest="sample_file", 
        help="The list file with samples for studying")
    parser.add_option("-o", "--output-file", dest="feature_file", 
        help="The output file with extracted features", default=None)
    parser.add_option("-d", "--data-dir", dest="data_dir", 
        help="The diroctory of the data for feature extraction", default=None)

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default="1",
        help="Number of subprocesses [default: %default]")
    group.add_option("--no-PTM-check", action="store_true", dest="no_PTM_check", 
        default=False, help="Don't check residue match for disorder file")
    group.add_option("--verbose", action="store_true", dest="verbose", 
        default=False, help="Print out log info")
    parser.add_option_group(group)
    
    # main arguments
    (options, args) = parser.parse_args()
    data_dir = options.data_dir
    sample_file = options.sample_file
    feature_file = options.feature_file
    if data_dir is None:
        data_dir = os.path.dirname(sample_file) + "/../../data"
    if os.access(data_dir, os.R_OK) == False:
        print("[PTM-X] Error: data-dir is not readable. Please check!")
        exit()
    if feature_file is None:
        feature_file = os.path.dirname(sample_file) + "/sample_feature.tsv"
    if os.access(os.path.dirname(feature_file), os.W_OK) == False:
        print("[PTM-X] Error: output-file is not writable. Please check!")
        exit()

    # optional arguments
    nproc = options.nproc
    verbose = options.verbose
    residue_check = (options.no_PTM_check == False)
    
    # load sample for feature obtaining
    samples = np.genfromtxt(sample_file, delimiter='\t', skip_header=0, 
        dtype="str")

    # fetch features
    feature = np.zeros((samples.shape[0], 4))

    ## sequence co-evolution features
    feature[:,0:2] = fetch_seqCoEvol(samples, data_dir, verbose, nproc)

    ## modification co-evolution features
    PTM_species_file = data_dir + "/PTMsites/PTM_data_3species.txt"
    feature[:,2] = fetch_PTMcoEvol(samples, PTM_species_file, verbose=verbose)

    ## co-occurence features
    coocur_file = data_dir + "/co-occur/PTM_coocurrence_data.txt"
    feature[:,3] = fetch_coOccur(samples, coocur_file, verbose, nproc)[:,1]

    # save features
    fid = open(feature_file, "w")
    head_line = "Prot_id1\tResidue1\tPTM1\tProt_id2\tResidue2\tPTM2\t"
    head_line = head_line + "seqCoEvol\tmotifCoCons\tPTMcoCons\tcoOcurr\n"
    fid.writelines(head_line)
    for i in range(samples.shape[0]):
        val_out = [x for x in samples[i,:]] + ["%.3f" %x for x in feature[i,:]]
        fid.writelines("\t".join(val_out) + "\n")
    fid.close()

    print("[PTM-X] saved features into file!")


if __name__ == '__main__':
    main()
