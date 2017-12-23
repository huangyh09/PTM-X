# This function is to get the co-evolution scores on sequence, motif and 
# modification existence. 

import os
import numpy as np
import scipy.stats as st
import multiprocessing
from .base import FastaMSA, id_mapping
from sklearn.metrics import hamming_loss
from sklearn.metrics.cluster import normalized_mutual_info_score


def get_msa_site(pro_id, seq, site, res, verbose=True):
    """TODO: Please check the surrounding residues
    """
    cnt = -1
    msa_site = -1
    for i in range(len(seq)):
        if seq[i] != "-": # and seq[i] != "?" # and seq[i] != "X"
            cnt += 1
        if cnt == site:
            msa_site = i
            break
    if msa_site == -1:
        if verbose:
            print("%s is too short for %s%d." %(pro_id, res, site+1))
    elif seq[msa_site] != res:
        #check near sites
        for ii in [-1, 1, -2, 2]:
            near_site = msa_site+ii
            if near_site>=0 and near_site<len(seq) and seq[near_site]==res:
                msa_site = near_site
                break
        if verbose and seq[msa_site] != res:
            idx_min = max(msa_site-4, 0)
            idx_max = min(msa_site+5, len(seq))
            ref_motif = seq[idx_min : idx_max]
            print("%s %s%d doesn't match %s" %(pro_id, res, site+1, ref_motif))
        msa_site = -1
    return msa_site


def get_msa(NOG_id, Ensm_id, PTM, rm_paralog=True, verbose=True, 
        align_dir="../../data/eggNOG4/veNOG/", motif_half=3):
    """
    """
    msa_seq, msa_species, msa_prot, motif_cons = [], [], [], []
    fasta_file = "%s/veNOG.%s.meta_raw.fa" %(align_dir, NOG_id)

    if os.path.isfile(fasta_file) == False:
        if verbose:
            print("No file of %s.fa in the path %s." %(NOG_id, align_dir))
        return msa_seq, msa_species, msa_prot, motif_cons

    # process the fasta file. 
    # Please remove paralog at this level. 
    names, seqs, species = readFastaEntry(fasta_file, Ensm_id, rm_paralog)

    # get the sequence of pro_id in use
    pro_idx = np.where(names == ("9606."+Ensm_id))[0][0]
    human_seq = seqs[pro_idx]
    seq_len = len(human_seq) - human_seq.count("-")

    site, res = int(PTM[1:])-1, PTM[0]
    msa_site = get_msa_site(Ensm_id, human_seq, site, res, verbose)
    # if msa_site == -1:
    #     if verbose:
    #         names, seqs, species = readFastaEntry(fasta_file, Ensm_id, False)
    #         idx = np.where(species == "9606")[0]
    #         if len(idx) >= 2:
    #             print(np.array(names)[idx])
    #             for ii in idx:
    #                 msa_site = get_msa_site(names[ii], seqs[ii], site, res, False)
    #                 if msa_site != -1:
    #                     print("%d Ensemble proteins: found PTM in other protein" 
    #                         %(len(idx)))
    #         #print(species)

    #         # print("The %s of %s doesn't match in MSA file" %(PTM, Ensm_id))
    #         pass
    if msa_site != -1:
        msa_seq = [_seq[msa_site] for _seq in seqs]
        msa_prot = [_name.split(".")[1] for _name in names]
        msa_species = species

        motif_start = max(msa_site-motif_half, 0)
        motif_stop  = min(msa_site+motif_half, len(human_seq)-1)
        
        ref_motif = []
        ref_motif[:] = human_seq[motif_start:motif_stop+1]

        for _seq in seqs:
            _motif = []
            _motif[:] = _seq[motif_start:motif_stop+1]
            motif_cons.append(np.mean(np.array(ref_motif) == np.array(_motif)))
    
    return msa_seq, msa_species, msa_prot, motif_cons


def get_seqCoEvol(uniprot_pair, PTM_pair, veNOG_Ensm_Unip, align_dir=None, 
    rm_paralog=True, verbose=True):
    """
    example: T15 on acc1 and K28 on acc2:
    uniprot_pair=[acc1, acc2], PTM_pair=[T15, K28]
    """
    # map protein ids
    idx_both = [np.where(veNOG_Ensm_Unip[:,2] == uniprot_pair[0])[0],
                np.where(veNOG_Ensm_Unip[:,2] == uniprot_pair[1])[0]]
    for i in [0,1]:
        if len(idx_both[i]) == 0:
            if verbose:
                print("No MSA file in veNOG of protein %s!" %uniprot_pair[i])
            return None, None
        if len(idx_both[i]) > 1:
            if verbose:
                print("Multiple veNOG files for %s!" %uniprot_pair[i])

    # get seq alignment
    for ii in idx_both[0]:
        veNOG_id, Ensem_id = veNOG_Ensm_Unip[ii, 0:2]
        msa_seq1, msa_species1, msa_prot1, motif_cons1 = get_msa(veNOG_id, 
            Ensem_id, PTM_pair[0], rm_paralog, verbose, align_dir)
        if len(msa_seq1) > 0:
            break
    for ii in idx_both[1]:
        veNOG_id, Ensem_id = veNOG_Ensm_Unip[ii, 0:2]
        msa_seq2, msa_species2, msa_prot2, motif_cons2 = get_msa(veNOG_id, 
            Ensem_id, PTM_pair[1], rm_paralog, verbose, align_dir)
        if len(msa_seq2) > 0:
            break

    # shared species
    _idx = id_mapping(msa_species1, msa_species2, uniq_ref_only=False)
    idx1 = np.arange(len(_idx))[_idx>=0]
    idx2 = _idx[_idx>=0]

    msa_seq_use1 = np.array(msa_seq1)[idx1.astype(int)]
    msa_seq_use2 = np.array(msa_seq2)[idx2.astype(int)]
    motif_cons1 = np.array(motif_cons1)[idx1.astype(int)]
    motif_cons2 = np.array(motif_cons2)[idx2.astype(int)]

    # get coEvolution
    if len(idx1) > 0:
        msa_seq_use1 = msa_seq_use1 == PTM_pair[0][0]
        msa_seq_use2 = msa_seq_use2 == PTM_pair[1][0]
        hamming_kept = np.mean(msa_seq_use1 == msa_seq_use2)
        motif_cocnsv = np.mean(motif_cons1 * motif_cons2)

        ## alternative methods (used)
        # nMI = hamming_loss(msa_seq_use1, msa_seq_use2)
        # nMI = np.mean(msa_seq_use1 * msa_seq_use2)
        # nMI = st.pearsonr(motif_cons1, motif_cons2)[0]
        # nMI = normalized_mutual_info_score(msa_seq_use1, msa_seq_use2)
    else:
        hamming_kept, motif_cocnsv = None, None
    return hamming_kept, motif_cocnsv


def fetch_seqCoEvol(samples, data_dir, verbose=False, nproc=1):
    align_dir = data_dir + "/eggNOG4/veNOG_raw_algs/"
    veNOG_Ensm_Unip_file = data_dir + "/eggNOG4/veNOG_Ensembl_Uniprot.txt"
    veNOG_Ensm_Unip = np.genfromtxt(veNOG_Ensm_Unip_file, delimiter='\t', 
        skip_header=1, dtype="str")
    
    coevol_val = np.zeros((samples.shape[0], 2))
    if nproc <= 1:
        for i in range(len(coevol_val)):
            coevol_val[i,:] = get_seqCoEvol(samples[i,[0,3]], samples[i,[1,4]], 
                veNOG_Ensm_Unip, align_dir, rm_paralog=True, verbose=verbose)
    else:
        pool = multiprocessing.Pool(processes=nproc)
        result = []
        for i in range(len(coevol_val)):
            result.append(pool.apply_async(get_seqCoEvol, (samples[i,[0,3]], 
                samples[i,[1,4]], veNOG_Ensm_Unip, align_dir, True, verbose)))
        pool.close()
        pool.join()
        coevol_val = np.array([res.get() for res in result], dtype=float)

    idx = coevol_val[:,0] == coevol_val[:,0]
    mean_val = np.mean(coevol_val[idx,0])
    print("[PTM-X] fetched seq co-evolution for %d samples. mean: %.3f, "
          "nan: %.1f%%." %(len(samples), mean_val, 100-np.mean(idx)*100))
    return coevol_val


def readFastaEntry(fasta_file, pro_id="", rm_paralog=False):
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


def fetch_PTMcoEvol(samples, PTM_species_file, shift_loc=[0, 1, -1], 
    verbose=False):
    PTMcoEvol_val = np.zeros(samples.shape[0])
    samples_ids1 = [x[0]+":"+x[1] for x in samples[:,0:2]]
    samples_ids2 = [x[0]+":"+x[1] for x in samples[:,3:5]]

    PTM_dat = np.genfromtxt(PTM_species_file, delimiter='\t', 
        skip_header=1, dtype="str")
    PTM_loc = np.array([x[1:] for x in PTM_dat[:, 1]]).astype(int)
    PTM_ids_full = []
    for _shift in shift_loc:
        PTM_ids = []
        for ii in range(len(PTM_dat)):
            _PTM_id = PTM_dat[ii,0]+":"+PTM_dat[ii,1][0]+str(PTM_loc[ii]+_shift)
            PTM_ids.append(_PTM_id)
        PTM_ids_full.append(PTM_ids)
    
    idx1_full, idx2_full = [], []
    for i in range(len(shift_loc)):
        idx1 = id_mapping(samples_ids1, PTM_ids_full[i], uniq_ref_only=False)
        idx2 = id_mapping(samples_ids2, PTM_ids_full[i], uniq_ref_only=False)
        idx1_full.append(idx1)
        idx2_full.append(idx2)

    for i in range(len(PTMcoEvol_val)):
        idx_1 = np.array([x[i] for x in idx1_full])
        _idx1 = np.where(idx_1 >= 0)[0]

        idx_2 = np.array([x[i] for x in idx2_full])
        _idx2 = np.where(idx_2 >= 0)[0]
        if len(_idx1) == 0 or len(_idx2) == 0:
            PTMcoEvol_val[i] = None
        else:
            idx_1 = int(idx_1[np.min(_idx1)])
            idx_2 = int(idx_2[np.min(_idx2)])

            states1 = (PTM_dat[idx_1, [2,7,12]] != "nan").astype(float)
            states2 = (PTM_dat[idx_2, [2,7,12]] != "nan").astype(float)
            PTMcoEvol_val[i] = np.mean(states1 * states2)

    idx = PTMcoEvol_val == PTMcoEvol_val
    mean_val = np.mean(PTMcoEvol_val[idx])
    print("[PTM-X] fetched PTM co-evolution for %d samples. mean: %.3f, " 
          "nan: %.1f%%." %(len(samples), mean_val, 100-np.mean(idx)*100))
    return PTMcoEvol_val


    