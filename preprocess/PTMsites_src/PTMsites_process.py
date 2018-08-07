#!/usr/bin/python2.7
# run python PTMsites_process.py

import os
import numpy as np
from optparse import OptionParser

def load_file(PTM_file, species_name, keys):
    data = np.loadtxt(PTM_file, delimiter='\t', skiprows=3, dtype="str")
    key_idx = np.array([],"int")
    for i in range(len(keys)):
        _idx = np.where(data[0,:] == keys[i])[0]
        if _idx.shape[0] == 0:
            print("There is no keywords of %s " %keys[i] + "in the file %s!" %PTM_file)
        if _idx.shape[0] > 1:
            print("There is multiple keywords of %s " %keys[i] + "in the file %s!" %PTM_file)
        key_idx = np.append(key_idx, _idx[0])
    spc_idx = np.where(data[:,key_idx[0]]==species_name)[0]
    #print(np.unique(data[:,key_idx[0]], return_counts=True))
    RV = data[spc_idx, :][:, key_idx]
    return RV

if __name__ == '__main__':
    #0. parse command line options
    parser = OptionParser()
    parser.add_option("--data_dir",dest="data_dir", 
        help="The diroctory of the PTM sites data")
    parser.add_option("--file_list",dest="file_list", 
        help="The list file that contains the files waiting for processing")
    parser.add_option("--species",dest="species", 
        help="The species wanted to obtained from the full data, e.g., human or mouse")
    parser.add_option("--out_file",dest="out_file", 
        help="The file for saving processed data",default='untitled_PTMsite_file.txt')

    (options, args) = parser.parse_args()
    data_dir = options.data_dir
    file_list = options.file_list
    species = options.species
    out_file = options.out_file

    # define the keys that we will use
    keys = ["ORGANISM", "PROTEIN", "ACC_ID", "MOD_RSD", "MOD_RSD", "SITE_+/-7_AA", 
            "LT_LIT","MS_LIT","MS_CST"]

    #keys = ["ORG", "PROTEIN", "ACC_ID", "MOD_TYPE", "MOD_RSD", "MODSITE_SEQ",
    #        "PUBMED_LTP", "PUBMED_MS2", "CST_MS2"]


    # load the list file that contains the processing files
    fid = open(file_list,"r")
    all_files = fid.readlines()
    fid.close()

    # load all files that contains the PTM sites
    PTM_file = all_files[0].split()[0]
    PTMsites = np.array([], dtype="S50").reshape(-1, len(keys))
    for i in range(0,len(all_files)):
        PTM_file = all_files[i].split()[0]
        PTM_type = os.path.basename(all_files[i]).split("_")[0]
        PTMsites_tmp = load_file(os.path.join(data_dir, PTM_file), species, keys)
        PTMsites_tmp[:,3] = PTM_type
        PTMsites_tmp[:,4] = np.array([x.split("-")[0] for x in PTMsites_tmp[:,4]])
        PTMsites = np.append(PTMsites, PTMsites_tmp, axis=0)
        if(PTMsites_tmp.shape[0] > 0):
            print("%d %s for %s included!" %(len(PTMsites_tmp), PTM_type, species))
        else:
            print("%d %s for %s included!" %(0, PTM_type, species))

    # obtain the location of the PTMs, and sort the PTMs by the protein name, 
    # then by the PTM location
    rsd_loc = np.zeros(PTMsites.shape[0],"int")
    for i in range(PTMsites.shape[0]):
        rsd_loc[i] = PTMsites[i,4][1:]
    idx = np.lexsort((rsd_loc, PTMsites[:,1]))
    PTMsites = PTMsites[idx,:]

    # save the data into txt file
    keys[3] = "MOD_TYPE"
    fid = open(out_file,"w")
    key_line = "\t".join(keys) + "\n"
    fid.writelines(key_line)

    for i in range(PTMsites.shape[0]):
        data_line = "\t".join(list(PTMsites[i,:])) + "\n"
        fid.writelines(data_line)
    fid.close()
