# This function is to get the residue co-evolution score 
# for import: from fun2_coevolution import site_coevolve, modify_coevolve

import os
import numpy as np
from sklearn.metrics.cluster import *
from sklearn.metrics import hamming_loss

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


def get_seq_coEvol(uni_prot_pair, PTM_pair, rm_paralog=True, verbose=True, 
    data_dir="../../"):
    """
    standardize pro_id: Ensembl protein id; site1/2: int; res1/2: char
    output: save multiple sequence alignment for the two sites into a txt file
            return co-evolution score, float
            (optional)show the Phylogenetic Trees
    """
    RV = None

    # load the id mapping file
    veNOG_Ensm_Unip_file = data_dir + "/data/PTMsites/veNOG_Ensm_Unip.txt"
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
    align_dir = data_dir + "/data/eggNOG4/veNOG/"
    msa_seq1, msa_species1, msa_prot1 = get_msa(NOG_ids[0], Ensm_ids[0], 
        PTM_pair[0], rm_paralog, verbose, align_dir)
    msa_seq2, msa_species2, msa_prot2 = get_msa(NOG_ids[1], Ensm_ids[1], 
        PTM_pair[1], rm_paralog, verbose, align_dir)

    _idx = id_mapping(msa_species1, msa_species2) #shared species
    idx1 = np.arange(len(_idx))[_idx>=0]
    idx2 = _idx[_idx>=0]

    msa_seq_use1 = np.array(msa_seq1)[idx1.astype(int)]
    msa_seq_use2 = np.array(msa_seq2)[idx2.astype(int)]

    # print(np.array(msa_species1)[idx1.astype(int)])
    # print(np.array(msa_species2)[idx2.astype(int)])
    # print(msa_seq_use1)
    # print(msa_seq_use2)

    # get coEvolution
    if len(idx1) > 0:
        nMI = normalized_mutual_info_score(msa_seq_use1, msa_seq_use2)
    else:
        nMI = None
    return nMI


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
    
def get_align_site(pro_id, alignseq, site, res, data_dir="../../"):
    '''the pro_id is the uniport id here
    the site is the absolute site - 1'''
    RV = -1
    # load fasta file with uniprot id
    fasta_file = data_dir + "/data/humanProteinSeq/human_proseq_unip.fasta"
    fid = open(fasta_file,"r")
    all_lines = fid.readlines()
    fid.close()

    # find out the protein with uniprot id
    line_idx = -1
    for i in range(len(all_lines)):
        if all_lines[i][0] == ">" and all_lines[i][1:10].count(pro_id) == 1:
            line_idx = i
            break

    # get the uniprot sequence
    seq = ""
    for i in range(line_idx+1, len(all_lines)):
        if all_lines[i][0] == ">":
            break
        else :
            seq = seq + all_lines[i].split()[0]
    if seq == "":
        print("There is no protein with id as %s." %pro_id)
        return RV

    # check uniprot residue 
    if seq[site] != res:
        # print("The %d residue " %(site+1) + "on %s " %pro_id + "is not %s!" %res)
        return RV

    # map the residue to the aligned sequence
    seq_align = []
    seq_align[:] = alignseq
    seq_align = np.array(seq_align)
    idx_no_gap = np.where(seq_align != '-')[0]
    seq_align = seq_align[idx_no_gap]

    len_sur = 10
    left = min(len_sur, site)
    right = min(len(seq)-1-site, len_sur)

    for i in range(len(seq_align)):
        left_tmp = min(len_sur, i)
        right_tmp = min(len(seq_align)-1-i, len_sur)

        left_use = min(left, left_tmp)
        right_use = min(right, right_tmp)

        seq_sur_unip = []
        seq_sur_unip[:] = seq[site-left_use : site+right_use+1]

        seq_sur_alig = []
        seq_sur_alig[:] = seq_align[i-left_use : i+right_use+1]

        mapped_len = sum(np.array(seq_sur_unip) == np.array(seq_sur_alig))

        if (seq[site] == seq_align[i] and 
            mapped_len >= min(len(seq_sur_unip), len_sur+1)):
            if RV != -1:
                print "multiple matched surrouding sequence"
            RV = idx_no_gap[i]

    return RV



def site_coevolve(unip_id, site1, site2, res1, res2, rm_paralog=True, 
	data_dir="../../"):
    """
    standardize pro_id: Ensembl protein id; site1/2: int; res1/2: char
    output: save multiple sequence alignment for the two sites into a txt file
            return co-evolution score, float
            (optional)show the Phylogenetic Trees
    """
    RV = np.zeros(site1.shape[0])
    RV[:] = None

    # load the id mapping file
    veNOG_Ensm_Unip_file = data_dir + "/data/PTMsites/veNOG_Ensm_Unip.txt"
    veNOG_Ensm_Unip = np.loadtxt(veNOG_Ensm_Unip_file, delimiter='\t', 
    	skiprows=1, dtype="str")

    # id mapping and co-evolution
    idx_pro = np.where(veNOG_Ensm_Unip[:,2] == unip_id)[0]
    if idx_pro.shape[0] == 0:
        print("No MSA file in veNOG of protein %s!" %unip_id)
        return RV
    elif idx_pro.shape[0] > 1:
        if np.unique(veNOG_Ensm_Unip[idx_pro,0]).shape[0] == 1:
            pass
            # print("Multiple Ensembl ids in one veNOG file for %s!" %unip_id)
        else:
            print("Multiple veNOG files for %s!" %unip_id)
    # We will use the first one for simplisity
    NOG_id = veNOG_Ensm_Unip[idx_pro,0]
    Ensm_id = veNOG_Ensm_Unip[idx_pro,1]

    pro_id, NOG_id = Ensm_id[0], NOG_id[0]

    # load veNOG and Ensembl id map file
    #align_dir = "../data/eggNOG4/veNOG.align/"
    align_dir = data_dir + "/data/eggNOG4/veNOG/"

    # check the veNOG file
    if os.path.isfile(align_dir + NOG_id + ".fa") == False:
        print("No file of %s" %NOG_id + ".fa in the path %s" %align_dir)
        return RV

    # process the fasta file
    fasta_file = align_dir + NOG_id + ".fa"
    names, seqs, species = readFastaEntry(fasta_file, rm_paralog, pro_id)

    # get the sequence of pro_id in use
    pro_idx = np.where(names == ("9606."+pro_id))[0][0]
    human_seq = seqs[pro_idx]
    seq_len = len(human_seq) - human_seq.count("-")

    for i in range(site1.shape[0]):
        # check the site and residue in aligned sequences
        msa_site1 = get_msa_site(pro_id, human_seq, site1[i]-1, res1[i])
        msa_site2 = get_msa_site(pro_id, human_seq, site2[i]-1, res2[i])

        if msa_site1 == -1 or msa_site2 == -1:
            msa_site1 = get_align_site(unip_id, human_seq, site1[i]-1, res1[i], data_dir)
            msa_site2 = get_align_site(unip_id, human_seq, site2[i]-1, res2[i], data_dir)

        # if msa_site1 == -1 or msa_site2 == -1:
        #     print("The %dth or %dth site" %(site1[i]+1, site2[i]+1) + 
        #           " of %s" %pro_id + " does not match with residue.")
        #     continue

        unmatch = False
        if msa_site1 == -1:
            print("The %d residue %s of %s doesn't match in MSA file" 
                %(site1[i], res1[i], pro_id))
            unmatch = True
        if msa_site2 == -1:
            print("The %d residue %s of %s doesn't match in MSA file" 
                %(site2[i], res2[i], pro_id))
            unmatch = True
        if unmatch:
            continue

        # get the consevation sequences
        con_res1, con_res2 = [], []
        for j in range(len(names)):
            con_res1.append(seqs[j][msa_site1])
            con_res2.append(seqs[j][msa_site2])

        # get the site co-evolution score via normalized mutual information
        RV[i] = normalized_mutual_info_score(con_res1, con_res2)
        #print con_res1, con_res2

    return RV




def id_mapping(IDs1, IDs2):
    """
    Mapping IDs2 to IDs1, both of which should only contain unique ids.
    
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
            j += 1
        elif IDs1[idx1[i]] > IDs2[idx2[j]]:
            j += 1
            
    origin_idx = np.argsort(RV_idx1)
    RV_idx = np.array(RV_idx2)[origin_idx]
    return RV_idx
