#!/usr/bin/python2.7

# This file is to extract features for a given file of samples.
# The features include: 
# 1) disordered region 
# 2) co-evolution on sequence
# 3) co-evolution on modifications
# 4) domain interaction (from Rongting)

import os
import numpy as np
from optparse import OptionParser, OptionGroup
from utils.disorder import fetch_disorder
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
    samples = np.genfromtxt(sample_file, delimiter='\t', 
        skip_header=0, dtype="str")

    # fetch features
    feature = np.zeros((samples.shape[0], 4))
    feature[:,:] = None

    ## disorder features
    disorder_file = data_dir + "/disorder/result.txt"
    disorder_val = fetch_disorder(samples, disorder_file, residue_check, verbose)
    feature[:,0] = disorder_val["both_site_in"] #0 or 1

    ## sequence co-evolution features
    feature[:,1] = fetch_seqCoEvol(samples, verbose, data_dir)

    ## modification co-evolution features
    PTM_species_file = data_dir + "/PTMsites/PTMsites_3species.txt"
    feature[:,2] = fetch_PTMcoEvol(samples, PTM_species_file, verbose)

    ## co-occurence features
    #feature[:,3] = fetch_coOccur(samples, verbose)

    # save features
    fid = open(feature_file, "w")
    head_line = "Prot_id1\tResidue1\tPTM1\tProt_id2\tResidue2\tPTM2\t"
    head_line = head_line + "disordered\tseqCoEvol\tPTMcoEvol\tcoOcurr\n"
    for i in range(samples.shape[0]):
        val_out = [x for x in samples[i,:]] + ["%.3f" %x for x in feature[i,:]]
        fid.writelines("\t".join(val_out) + "\n")
    fid.close()

    print("[PTM-X] saved features into file!")


    # # please check the input!!!
    # unip_ids = np.unique(samples[:,0])

    # feature = np.zeros((samples.shape[0], 5)) 
    # feature[:] = None
    # for i in range(unip_ids.shape[0]):
    #     print(i+1, unip_ids.shape[0], unip_ids[i])

    #     # protein id
    #     unip_id = unip_ids[i]
    #     pro_idx = np.where(unip_ids[i] == samples[:,0])[0]

    #     # residues and locations
    #     res1 = np.zeros(pro_idx.shape[0], "str")
    #     res2 = np.zeros(pro_idx.shape[0], "str")
    #     loc1 = np.zeros(pro_idx.shape[0], "int")
    #     loc2 = np.zeros(pro_idx.shape[0], "int")

    #     for j in range(pro_idx.shape[0]):
    #         res1[j] = samples[pro_idx[j],1][0]
    #         res2[j] = samples[pro_idx[j],3][0]

    #         loc1[j] = samples[pro_idx[j],1][1:]
    #         loc2[j] = samples[pro_idx[j],3][1:]

    #     # PTM types
    #     PTM1 = samples[pro_idx,2]
    #     PTM2 = samples[pro_idx,4]

    #     # sequence distance
    #     dist_seq2 = dist_1d(unip_id, loc1, loc2, res1, res2)
    #     feature[pro_idx, 0] = dist_seq2

    #     # structual distance
    #     dist_str2 = dist_3d(unip_id, loc1, loc2, res1, res2)
    #     feature[pro_idx, 1] = dist_str2

    #     # disorder location
    #     disoreder_bool = is_disorder(unip_id, loc1, loc2, res1, res2, threshold=2)
    #     feature[pro_idx, 2] = disoreder_bool

    #     # site co-evolution
    #     site_ce_val = site_coevolve(unip_id, loc1, loc2, res1, res2)
    #     feature[pro_idx, 3] = site_ce_val

    #     # modificatio co-evolution
    #     modi_ce_val = modify_coevolve(unip_id, loc1, loc2, res1, res2)
    #     feature[pro_idx, 4] = modi_ce_val

    # print "Feature shape is:" 
    # print (feature.shape)

    # # save the features
    # f = h5py.File(feature_file, "w")
    # f["feature"] = feature
    # f.close()


if __name__ == '__main__':
    main()
