# This function is to get the residue co-evolution score 
# for import: from fun2_coevolution import site_coevolve, modify_coevolve

import os
import numpy as np
from sklearn.metrics import hamming_loss
from sklearn.metrics.cluster import normalized_mutual_info_score

#TODO: check X for unknown amino acid.

def get_msa_site(pro_id, seq, site, res, verbose=True):
    cnt = -1
    msa_site = -1
    for i in range(len(seq)):
        if seq[i] != "-" and seq[i] != "?":
            cnt += 1
        if cnt == site:
            msa_site = i
            break
    if msa_site == -1:
        if verbose:
            print(str(site + 1) + " is out the length of %s." %pro_id)
    elif seq[msa_site] != res:
        #print("The %d site" %(site+1) + " of %s" %pro_id + " is %s." %seq[msa_site])
        msa_site = -1
    return msa_site


def get_msa(NOG_id, Ensm_id, PTM, rm_paralog=True, verbose=True, 
        align_dir="../../data/eggNOG4/veNOG/"):

    msa_seq, msa_species, msa_prot = [], [], []  
    # check the veNOG file
    # print(align_dir, NOG_id)
    if os.path.isfile(align_dir + NOG_id + ".fa") == False:
        if verbose:
            print("No file of %s.fa in the path %s." %(NOG_id, align_dir))
        return msa_seq, msa_species, msa_prot

    # process the fasta file
    fasta_file = align_dir + NOG_id + ".fa"
    names, seqs, species = readFastaEntry(fasta_file, rm_paralog, Ensm_id)

    # get the sequence of pro_id in use
    pro_idx = np.where(names == ("9606."+Ensm_id))[0][0]
    human_seq = seqs[pro_idx]
    seq_len = len(human_seq) - human_seq.count("-")

    site, res = int(PTM[1:])-1, PTM[0]
    msa_site = get_msa_site(Ensm_id, human_seq, site, res, verbose)
    if msa_site == -1:
        if verbose:
            print("The %s of %s doesn't match in MSA file" %(PTM, Ensm_id))
    else:
        msa_seq = [_seq[msa_site] for _seq in seqs]
        msa_prot = [_name.split(".")[1] for _name in names]
        msa_species = species
    
    return msa_seq, msa_species, msa_prot


def get_seq_nMI(uni_prot_pair, PTM_pair, rm_paralog=True, verbose=True, 
    data_dir="../../data/"):
    """
    standardize pro_id: Ensembl protein id; site1/2: int; res1/2: char
    output: save multiple sequence alignment for the two sites into a txt file
            return co-evolution score, float
            (optional)show the Phylogenetic Trees
    """
    RV = None

    # load the id mapping file
    veNOG_Ensm_Unip_file = data_dir + "/eggNOG4/veNOG_Ensm_Unip.txt"
    veNOG_Ensm_Unip = np.loadtxt(veNOG_Ensm_Unip_file, delimiter='\t', 
        skiprows=1, dtype="str")

    # map protein ids
    idx_both = [np.where(veNOG_Ensm_Unip[:,2] == uni_prot_pair[0])[0],
                np.where(veNOG_Ensm_Unip[:,2] == uni_prot_pair[1])[0]]
    for i in [0,1]:
        if len(idx_both[i]) == 0:
            if verbose:
                print("No MSA file in veNOG of protein %s!" %uni_prot_pair[i])
            return RV
        if len(idx_both[i]) > 1:
            if verbose:
                print("Multiple veNOG files for %s!" %uni_prot_pair[i])

    NOG_ids  = [str(veNOG_Ensm_Unip[idx_both[0], 0][0]),
                str(veNOG_Ensm_Unip[idx_both[1], 0][0])]
    Ensm_ids = [str(veNOG_Ensm_Unip[idx_both[0], 1][0]),
                str(veNOG_Ensm_Unip[idx_both[1], 1][0])]

    # get seq alignment
    align_dir = data_dir + "/eggNOG4/veNOG/"
    msa_seq1, msa_species1, msa_prot1 = get_msa(NOG_ids[0], Ensm_ids[0], 
        PTM_pair[0], rm_paralog, verbose, align_dir)
    msa_seq2, msa_species2, msa_prot2 = get_msa(NOG_ids[1], Ensm_ids[1], 
        PTM_pair[1], rm_paralog, verbose, align_dir)

    # shared species
    _idx = id_mapping(msa_species1, msa_species2, uniq_ref_only=False)
    idx1 = np.arange(len(_idx))[_idx>=0]
    idx2 = _idx[_idx>=0]

    msa_seq_use1 = np.array(msa_seq1)[idx1.astype(int)]
    msa_seq_use2 = np.array(msa_seq2)[idx2.astype(int)]

    # get coEvolution
    if len(idx1) > 0:
        msa_seq_use1 = msa_seq_use1 == PTM_pair[0][0]
        msa_seq_use2 = msa_seq_use2 == PTM_pair[1][0]
        nMI = normalized_mutual_info_score(msa_seq_use1, msa_seq_use2)

        # print(np.array(msa_species1)[idx1.astype(int)])
        # print(np.array(msa_species2)[idx2.astype(int)])
        # print(msa_seq_use1)
        # print(msa_seq_use2)
        
        # print(msa_species1[idx1.astype(int)])
        # print(msa_species2[idx2.astype(int)])
    else:
        nMI = None
    return nMI


def fetch_seqCoEvol(samples, verbose, data_dir):
    coevol_val = np.zeros(samples.shape[0])

    for i in range(len(coevol_val)):
        coevol_val[i] = get_seq_nMI(samples[i,[0,3]], samples[i,[1,4]], 
            verbose=verbose, data_dir=data_dir)

    idx = coevol_val == coevol_val
    mean_val = np.mean(coevol_val[idx])
    print("[PTM-X] fetched sequence co-evolution for %d samples. mean: %.3f, "
          "NA: %.1f%%." %(len(samples), mean_val, 100-np.mean(idx)*100))
    return coevol_val


def readFastaEntry(fasta_file, rm_paralog=False, pro_id=""):
    # open file and read lines
    fid = open(fasta_file,"r")
    all_lines = fid.readlines()
    fid.close()

    # process all lines
    names, species, sequences = [], [], []
    seq = ""
    for i in range(len(all_lines)):
        line = all_lines[i].split()[0]
        if line.startswith( ">" ):
            names.append(line.split(">")[1])
            species.append(line.split(">")[1].split(".")[0])
            if i == 0:
                continue
            sequences.append(seq)
            seq = ""
        else:
            seq = seq + line
    sequences.append(seq)
    names = np.array(names)
    species = np.array(species)
    sequences = np.array(sequences)
    if rm_paralog == False:
        return names, sequences, species

    # remove paralog data
    human_idx = np.where(names == ("9606."+pro_id))[0][0]
    human_seq = []
    human_seq[:0] = str(sequences[human_idx])
    species_uni = np.unique(species)
    kp_idx = np.zeros(species_uni.shape[0],"int")
    for i in range(species_uni.shape[0]):
        if species_uni[i] == "9606":
            kp_idx[i] = human_idx
        else :
            # choose the most similar paralogous protein via hamming distance
            _idx = np.where(np.array(species)==species_uni[i])[0]
            if _idx.shape[0] > 0:
                hm_loss = np.zeros(_idx.shape[0],"float")
                for j in range(_idx.shape[0]):
                    seq_tmp = []
                    seq_tmp[:0] = str(sequences[_idx[j]])
                    hm_loss[j] = hamming_loss(human_seq, seq_tmp)
                _idx = _idx[np.argsort(hm_loss)]
            kp_idx[i] = _idx[0]
            # kp_idx[i] = species.index(species_uni[i])
    return names[kp_idx], sequences[kp_idx], species[kp_idx]
    

def id_mapping(IDs1, IDs2, uniq_ref_only=True):
    """
    Mapping IDs2 to IDs1. IDs1 (ref id) can have repeat values, but IDs2 need 
    to only contain unique ids.
    Therefore, IDs2[rv_idx] will be the same as IDs1.
    
    Parameters
    ----------
    IDs1 : array_like or list
        ids for reference.
    IDs2 : array_like or list
        ids waiting to map.
        
    Returns
    -------
    RV_idx : array_like, the same length of IDs1
        The index for IDs2 mapped to IDs1. If an id in IDs1 does not exist 
        in IDs2, then return a None for that id.
    """
    idx1 = np.argsort(IDs1)
    idx2 = np.argsort(IDs2)
    RV_idx1, RV_idx2 = [], []
    
    i, j = 0, 0
    while i < len(idx1):
        if j == len(idx2) or IDs1[idx1[i]] < IDs2[idx2[j]]:
            RV_idx1.append(idx1[i])
            RV_idx2.append(None)
            i += 1
        elif IDs1[idx1[i]] == IDs2[idx2[j]]:
            RV_idx1.append(idx1[i])
            RV_idx2.append(idx2[j])
            i += 1
            if uniq_ref_only: 
                j += 1
        elif IDs1[idx1[i]] > IDs2[idx2[j]]:
            j += 1
            
    origin_idx = np.argsort(RV_idx1)
    RV_idx = np.array(RV_idx2)[origin_idx]
    return RV_idx



def fetch_PTMcoEvol(samples, PTM_species_file, verbose):
    PTMcoEvol_val = np.zeros(samples.shape[0])

    data = np.genfromtxt(PTM_species_file, delimiter='\t', 
        skip_header=1, dtype="str")
    PTM_val = data[:, 3:6].astype(float)
    PTM_ids = [x[0]+":"+x[1] for x in data[:,0:2]]

    samples_ids1 = [x[0]+":"+x[1] for x in samples[:,0:2]]
    samples_ids2 = [x[0]+":"+x[1] for x in samples[:,3:5]]

    idx1 = id_mapping(samples_ids1, PTM_ids, uniq_ref_only=False)
    idx2 = id_mapping(samples_ids2, PTM_ids, uniq_ref_only=False)

    for i in range(len(PTMcoEvol_val)):
        if idx1[i] >= 0 and idx2[i] >= 0:
            PTMcoEvol_val[i] = np.mean(PTM_val[int(idx1[i])] * 
                                       PTM_val[int(idx2[i])])
        else:
            PTMcoEvol_val[i] = None

    idx = PTMcoEvol_val == PTMcoEvol_val
    mean_val = np.mean(PTMcoEvol_val[idx])
    print("[PTM-X] fetched PTM co-evolution for %d samples. mean: %.3f, " 
          "NA: %.1f%%." %(len(samples), mean_val, 100-np.mean(idx)*100))
    return PTMcoEvol_val