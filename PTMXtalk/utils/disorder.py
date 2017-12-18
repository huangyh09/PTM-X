# This function is to check the whether the sites are located in disordered 
# regions and both in a same disordered region
# for import: from fun3_disorder_v1 import is_disorder

import numpy as np

def load_disorder_file(disorder_file):
    """
    This function is to load and process the disorder prediction results file.
    It will return each disorder segments with the protein ids, no matter from
    which prediction method the segments comes.
    """
    fid = open(disorder_file,"r")
    all_lines = fid.readlines()
    fid.close()

    dis_prot_ID, dis_regions, dis_seq = [], [], []
    pro_tmp = None
    for i in range(len(all_lines)):
        # ignore comments and sequences
        if all_lines[i].count("> sp|") != 1:
            continue
        
        a_line = all_lines[i].rstrip().split(" ")
        curr_pro = a_line[1].split("|")[1]

        # first protein
        if pro_tmp is None:
            pro_tmp = curr_pro
            seq_tmp = all_lines[i+1].rstrip()
            reg_tmp = [x.split(",")[0].split("-") for x in a_line[2:]]

        # the same protein
        elif pro_tmp == curr_pro:
            reg_tmp += [x.split(",")[0].split("-") for x in a_line[2:]]

        # a new protein
        elif pro_tmp != curr_pro:
            dis_prot_ID.append(pro_tmp)
            dis_regions.append(reg_tmp)
            dis_seq.append(seq_tmp)
            pro_tmp = curr_pro
            seq_tmp = all_lines[i+1].rstrip()
            reg_tmp = [x.split(",")[0].split("-") for x in a_line[2:]]

        # last protein
        if i == len(all_lines)-2:
            dis_prot_ID.append(pro_tmp)
            dis_regions.append(reg_tmp)
            dis_seq.append(seq_temp)

    return dis_prot_ID, dis_regions, dis_seq


def count_disorder(PTM, dis_regions, prot_seq=None, verbose=True):
    """
    This function is to count how many times a protein residue is predicted to 
    be in a disordered region. There are three predictors, thus the output can 
    be 0, 1, 2, 3. 

    If the prot_seq is supplied, check the residue. Otherwise, not.
    """
    PTM_res = PTM[0]
    PTM_site = int(PTM[1:])
    dis_regions = np.array(dis_regions, dtype=int)

    # check the residue of the PTM, if given prot_seq
    if prot_seq is not None:
        if PTM_site >= len(prot_seq):
            if verbose: print("PTM sites out side of protein sequence.")
            return None
        if prot_seq[PTM_site-1].upper() != PTM_res:
            if verbose: print("PTM residue doesn't match the sequence.")
            return None

    # count the disordered regions
    idx = (PTM_site >= dis_regions[:,0]) * (PTM_site <= dis_regions[:,1])
    count = sum(idx)

    return count

def fetch_disorder(samples, dis_file="../../data/disorder/result.txt", 
    residue_check=False, verbose=True):
    """
    Fetch the counts for disorders for a set of PTM crosstalk samples.
    """
    # load disorder file
    dis_prot_ID, dis_regions, dis_seq = load_disorder_file(dis_file)

    # fetch counts
    dis_count = np.zeros((len(samples), 2), dtype=float)
    for i in range(len(dis_count)):
        for j in range(2):
            _pro = samples[i,0+3*j]
            _PTM = samples[i,1+3*j]

            _cnt = dis_prot_ID.count(_pro)
            if _cnt == 0:
                dis_count[i,j] = None
            else:
                _idx = dis_prot_ID.index(_pro)
                if residue_check:
                    prot_seq = dis_seq[_idx]
                else:
                    prot_seq = None
                dis_count[i,j] = count_disorder(_PTM, dis_regions[_idx], 
                    prot_seq, verbose)
            if _cnt > 1 and verbose:
                print("mutiple %s items in the file, use the first here." %_pro)

    # check if both sites in disordered region
    threshold = 2
    both_in = np.zeros(samples.shape[0])
    idx_mat = dis_count == dis_count
    idx_use = idx_mat[:,0] * idx_mat[:,1]
    both_in[idx_use] = ((dis_count[idx_use,:] >= threshold).sum(axis=1) == 2)

    mean_val = np.mean(both_in[idx_use])
    print("[PTM-X] fetched disordered regions for %d samples. mean: %.3f, "
        "NA: %.1f%%." %(len(samples), mean_val, 100-np.mean(idx_mat)*100))
    
    RV = {}
    RV["disorder_count"] = dis_count
    RV["both_site_in"] = both_in
    return RV

