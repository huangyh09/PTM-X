#!/usr/bin/python2.7

# This file is to generate positive and negative samples of PTM cross-talk paris 
# The positive samples are directly from the collected data set, the negative 
# samples will be generated from the same set of protein pairs as positive 
# samples.


import numpy as np
from optparse import OptionParser

def main():
    #0. parse command line options
    parser = OptionParser()
    parser.add_option("--data_dir", "-d", dest="data_dir", 
        help="The diroctory of collected crosstalk samples",
        default="../interface/crosstalkSamples/")
    (options, args) = parser.parse_args()
    data_dir = options.data_dir

    PTM_tmp = np.genfromtxt(data_dir + "../../data/PTMsites/PTMsites_human.txt", 
        delimiter='\t', skip_header=1, dtype="str")
    idx_tmp = np.where(PTM_tmp[:, 6] != "")[0] #PTM with literature support
    PTMsite = PTM_tmp[idx_tmp, :]

    # load released crosstalk data
    # idx_use: protein1_id, res1, PTM1, protein2_id, res2, PTM2, pro_pair
    idx_use = [1,2,3,5,6,7,8]     
    ct_data =  np.genfromtxt(data_dir + "PTM_crosstalk_Release.txt", 
        delimiter='\t', skip_header=1, dtype="str")[:,idx_use]

    pro_ids, _idx = np.unique(ct_data[:,6], return_index=True)

    # generating negative samples
    samp_data = []
    for i in range(pro_ids.shape[0]):
        pro1, pro2 = pro_ids[i].split('+')

        # cross-talk PTM sites
        ct_idx = np.where(ct_data[:,6] == pro_ids[i])[0]
        ct_sites1 = np.unique(ct_data[ct_idx,1])
        ct_sites2 = np.unique(ct_data[ct_idx,4])
        
        # map all PTM with the same proteins
        pro_idx1 = np.where(PTMsite[:, 2] == pro1)[0]
        uni_site1, _idx1 = np.unique(PTMsite[pro_idx1, 4], return_index=True)
        PTM_use1 = PTMsite[pro_idx1[np.sort(_idx1)]]

        pro_idx2 = np.where(PTMsite[:, 2] == pro2)[0]
        uni_site2, _idx2 = np.unique(PTMsite[pro_idx2, 4], return_index=True)
        PTM_use2 = PTMsite[pro_idx2[np.sort(_idx2)]]

        #combine PTM pair
        for h in range(PTM_use1.shape[0]):
            for k in range(PTM_use2.shape[0]):
                if PTM_use1[h,4] in ct_sites1 and PTM_use2[k,4] in ct_sites2:
                    continue

                samp_temp = list(PTM_use1[h, [2,4,3]]) + list(PTM_use2[k, [2,4,3]])
                samp_temp[2] = samp_temp[2].lower()
                samp_temp[5] = samp_temp[5].lower()
                samp_data.append(samp_temp)

    # save postive and negative the samples
    fid = open(data_dir + "positive_samples.txt", "w")
    for i in range(ct_data.shape[0]):
        fid.writelines("\t".join(list(ct_data[i, :6])) + "\n")
    fid.close()

    fid = open(data_dir + "negative_samples.txt", "w")
    for i in range(len(samp_data)):
        fid.writelines("\t".join(samp_data[i]) + "\n")
    fid.close()

    print("[PTM-X] samples generated: %d postive and %d negative." 
        %(ct_data.shape[0], len(samp_data)))


if __name__ == '__main__':
    main()

