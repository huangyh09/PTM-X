#!/usr/bin/python2.7
# run python sample_generator.py

import numpy as np
from optparse import OptionParser

if __name__ == '__main__':
    #0. parse command line options
    parser = OptionParser()
    
    parser.add_option("--sample_type",dest="sample_type", help="The type of sample",default="positive")
    parser.add_option("--data_dir",dest="data_dir", help="The diroctory of the PTM sites data",
                      default="../interface/crosstalkSamples/")
    parser.add_option("--out_file",dest="out_file", help="The file for saving processed data",
                      default='untitled_samples.txt')

    (options, args) = parser.parse_args()
    data_dir = options.data_dir
    out_file = options.out_file
    sample_type = options.sample_type

    # load human PTM data
    PTM_dir = "../data/PTMsites/"
    PTMsite = np.loadtxt(PTM_dir + "PTMsites_human.txt", delimiter='\t', skiprows=1, dtype="str")

    # load released crosstalk data
    idx_use = [1,2,3,5,6,7,8]     # protein1_id, res1, PTM1, protein2_id,res2, PTM2,pro_pair
    ct_data =  np.loadtxt(data_dir + "PTM_crosstalk_Release.txt", delimiter='\t', skiprows=1, dtype="str")[:,idx_use]
    
    pro_ids, _idx = np.unique(ct_data[:,6],return_index=True)

    # collect samples
    samp_data = np.zeros((0,6),"str")
    if sample_type == "positive":
    	idp_use = [0,1,2,3,4,5]
        samp_data = ct_data[:,idp_use]
    else :
        for i in range(pro_ids.shape[0]):
            # cross-talk PTM sites
            ct_idx = np.where(ct_data[:,6] == pro_ids[i])[0]
            ct_sites1 = ct_data[ct_idx,1]
            ct_sites2 = ct_data[ct_idx,4]
            ct_sites1 = np.unique(ct_sites1)
            ct_sites2 = np.unique(ct_sites2)

            
            pro1=pro_ids[i].split('+')[0]
            pro2=pro_ids[i].split('+')[1]

            # all PTM sites
            all_idx1 = np.where(PTMsite[:,2] == pro1)[0]
            all_idx2 = np.where(PTMsite[:,2] == pro2)[0]
            # wanted PTM sites,only negative_pub is considered
            if sample_type == "negative_pub": 
                pub_idx1 = np.where(PTMsite[all_idx1,6] !="")[0]
                use_idx1 = pub_idx1
                pub_idx2 = np.where(PTMsite[all_idx2,6] !="")[0]
                use_idx2 = pub_idx2

            # if sample_type == "negative_pub2": 
            # 	pub_idx = np.where(PTMsite[all_idx,6] !="")[0]
            #     pub_idx2 = np.where(PTMsite[all_idx[pub_idx],6].astype("int") >= 2)[0]
            #     use_idx = pub_idx[pub_idx2]
            # if sample_type == "negative_pub3": 
            # 	pub_idx = np.where(PTMsite[all_idx,6] !="")[0]
            #     pub_idx3 = np.where(PTMsite[all_idx[pub_idx],6].astype("int") >= 3)[0]
            #     use_idx = pub_idx[pub_idx3]
            # elif sample_type == "negative_MS":
            #     MS1_idx = np.where(PTMsite[all_idx,7] !="")[0]
            #     MS2_idx = np.where(PTMsite[all_idx,8] !="")[0]
            #     use_idx = np.unique(np.append(MS1_idx, MS2_idx))
            # elif sample_type == "negative_all":
            #     pub_idx = np.where(PTMsite[all_idx,6] !="")[0]
            #     MS1_idx = np.where(PTMsite[all_idx,7] !="")[0]
            #     MS2_idx = np.where(PTMsite[all_idx,8] !="")[0]
            #     MS_idx = np.append(MS1_idx, MS2_idx)
            #     use_idx = np.unique(np.append(pub_idx, MS_idx))

            use_site1 = PTMsite[all_idx1[use_idx1],4]
            use_PTM1 = PTMsite[all_idx1[use_idx1],3]
            uni_val_1, _idx1 = np.unique(use_site1,return_index=True)
            use_site1 = use_site1[np.sort(_idx1)]
            use_PTM1 = use_PTM1[np.sort(_idx1)]

            use_site2 = PTMsite[all_idx2[use_idx2],4]
            use_PTM2 = PTMsite[all_idx2[use_idx2],3]
            uni_val_2, _idx2 = np.unique(use_site2,return_index=True)
            use_site2 = use_site2[np.sort(_idx2)]
            use_PTM2 = use_PTM2[np.sort(_idx2)]

            for h in range(use_site1.shape[0]):
            	for k in range(use_site2.shape[0]):
            		if (ct_sites==use_site1[h]).sum() + (ct_sites==use_site2[k]).sum() < 2:
            			samp_temp =[pro1,use_site1[h],use_PTM1[h],pro2,use_site2[k],use_PTM2[k]]
            			samp_temp = np.array(samp_temp).reshape(1,-1)
            			samp_data = np.append(samp_data, samp_temp, axis=0)
    print(samp_data.shape)

    # output the samples
    fid = open(data_dir + out_file,"w")
    for i in range(samp_data.shape[0]):
        fid.writelines("\t".join(list(samp_data[i,:])) + "\n")
    fid.close()