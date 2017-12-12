#!/usr/bin/python2.7

# This file is to merge the PTM existence of human, mouse and rat.
# The imputation of the PTM on rat is removed from here, as this can be better 
# don by using a proper classifier, e.g., random forest or neural net, by using 
# the sequence features.

import numpy as np
from optparse import OptionParser
from sklearn.metrics import hamming_loss

def readFastaEntry(fasta_file):
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
    names, sequences = np.array(names), np.array(sequences)
    return names, sequences


def get_msa_site(seq, sites):
    if len(sites) == 0:
        return np.array([],"int")

    loc, res = np.array([],"int"), np.array([],"str")
    for i in range(len(sites)):
        loc = np.append(loc, int(sites[i][1:]))
        res = np.append(res, sites[i][0])
    sort_idx = np.argsort(loc)
    loc = loc[sort_idx]
    res = res[sort_idx]

    cnt = 0
    msa_site = []
    site_idx = 0
    for i in range(len(seq)):
        if seq[i] != "-" :
            cnt += 1
        else :
            continue
        if cnt == loc[site_idx]:
            if seq[i] == res[site_idx]:
                msa_site.append(i+1)
            if site_idx == (len(sites)-1):
                break
            site_idx += 1
    return np.array(msa_site,"int")


def main():
    #0. parse command line options
    parser = OptionParser()
    parser.add_option("-l", "--seq_len", type="int", dest="seq_len", 
        help="The sequence around the PTM site", default="4")
    parser.add_option("-d", "--data_dir", dest="data_dir", 
        help="The directory of the data.")
    (options, args) = parser.parse_args()
    seq_len = options.seq_len
    data_dir = options.data_dir

    out_file = data_dir + "/PTMsites/PTMsites_3species.txt"

    #1. load the PTM data
    human_PTM_file = data_dir + "/PTMsites/PTMsites_human.txt"
    mouse_PTM_file = data_dir + "/PTMsites/PTMsites_mouse.txt"
    rat_PTM_file = data_dir + "/PTMsites/PTMsites_rat.txt"
    
    human_PTM = np.loadtxt(human_PTM_file, "str", delimiter="\t", skiprows=1)
    mouse_PTM = np.loadtxt(mouse_PTM_file, "str", delimiter="\t", skiprows=1)
    rat_PTM = np.loadtxt(rat_PTM_file, "str", delimiter="\t", skiprows=1)

    #2. orthologous data and ids
    msa_dir = data_dir + "/ortholog/InParanoid.align/"
    map_file = data_dir + "/ortholog/human_mouse_rat.txt"
    orth_map = np.loadtxt(map_file, "str", delimiter="\t", skiprows=1)
    pro_uni = orth_map[:,0]

    fid = open(out_file, "w")
    head_line = "pro_id\tresidue\tseq_len\thm_state\tms_state\trt_state\t"
    head_line = head_line + "ms_rt_sim\thm_rt_sim\thm_ms_sim\n"
    fid.writelines(head_line)

    for i in range(pro_uni.shape[0]):
        print(i,pro_uni[i])
        _idx = np.where(orth_map[:,0] == pro_uni[i])[0]
        if _idx.shape[0] == 0:
            continue

        names, seqs = readFastaEntry(msa_dir+pro_uni[i]+".afa")
        hm_pro = orth_map[_idx[0],0]
        ms_pro = orth_map[_idx[0],1]
        rt_pro = orth_map[_idx[0],2]
        hm_seq = seqs[names == hm_pro][0]
        ms_seq = seqs[names == ms_pro][0]
        rt_seq = seqs[names == rt_pro][0]

        hm_idx = np.where(human_PTM[:,2] == pro_uni[i])[0]
        ms_idx = np.where(mouse_PTM[:,2] == ms_pro)[0]
        rt_idx = np.where(rat_PTM[:,2] == rt_pro)[0]
        hm_PTM_sites = np.unique(human_PTM[hm_idx,4])
        ms_PTM_sites = np.unique(mouse_PTM[ms_idx,4])
        rt_PTM_sites = np.unique(rat_PTM[rt_idx,4])

        # this is the real site not that -1
        hm_msa_locs = get_msa_site(hm_seq, hm_PTM_sites) 
        ms_msa_locs = get_msa_site(ms_seq, ms_PTM_sites)
        rt_msa_locs = get_msa_site(rt_seq, rt_PTM_sites)

        # sites_uniq = np.unique(np.append(np.append(hm_msa_locs, ms_msa_locs), 
        # rt_msa_locs))
        sites_uniq = np.unique(hm_msa_locs)

        for j in range(sites_uniq.shape[0]):
            state_all = np.zeros(3,"int")
            seq_sim = np.zeros(3,"float")

            if (hm_msa_locs == sites_uniq[j]).sum() > 0:
                state_all[0] = 1
            if (ms_msa_locs == sites_uniq[j]).sum() > 0:
                state_all[1] = 1
            if (rt_msa_locs == sites_uniq[j]).sum() > 0:
                state_all[2] = 1

            left_seq, right_seq = [], []
            if sites_uniq[j]-1 >= seq_len:
                left_site = sites_uniq[j]-1 - seq_len
            else :
                left_site = 0
                left_seq = ["_"] * (seq_len - (sites_uniq[j]-1))

            if sites_uniq[j]-1 + seq_len <= len(seqs[0]) -1:
                right_site = sites_uniq[j]-1 + seq_len
            else :
                right_site = len(seqs[0]) -1
                right_seq = ["_"] * (sites_uniq[j] + seq_len - len(seqs[0]))

            hm_seq_tmp, ms_seq_tmp, rt_seq_tmp = [], [], []
            hm_seq_tmp[:0] = hm_seq[left_site:right_site+1]
            ms_seq_tmp[:0] = ms_seq[left_site:right_site+1]
            rt_seq_tmp[:0] = rt_seq[left_site:right_site+1]
            seqs_tmp = [left_seq + hm_seq_tmp + right_seq,
                        left_seq + ms_seq_tmp + right_seq,
                        left_seq + rt_seq_tmp + right_seq]

            seq_sim[0] = ((seqs_tmp[1][seq_len] == seqs_tmp[2][seq_len]) * 
                          (1-hamming_loss(seqs_tmp[1], seqs_tmp[2])))
            seq_sim[1] = ((seqs_tmp[2][seq_len] == seqs_tmp[0][seq_len]) * 
                          (1-hamming_loss(seqs_tmp[2], seqs_tmp[0])))
            seq_sim[2] = ((seqs_tmp[0][seq_len] == seqs_tmp[1][seq_len]) * 
                          (1-hamming_loss(seqs_tmp[0], seqs_tmp[1])))

            PTM_site_tmp = [hm_seq[sites_uniq[j]-1] + 
                str(sites_uniq[j]-hm_seq[:sites_uniq[j]].count("-"))]
            out_list  = [pro_uni[i]] + PTM_site_tmp
            out_list += [str(2*seq_len+1)] + ["%d" %x for x in state_all]
            out_list += ["%.2f" %x for x in seq_sim]

            fid.writelines("\t".join(out_list)+"\n")
    fid.close()

if __name__ == '__main__':
    main()
    
    