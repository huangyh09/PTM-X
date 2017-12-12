#!/usr/bin/python2.7
# run python ortholog2fasta.py

import numpy as np

def get_name_idx(fasta_lines):
    names = [fasta_lines[0][1:-1]]
    start = 0
    idx = []

    for i in range(1, len(fasta_lines)):
        if fasta_lines[i][0] == ">":
            idx.append([start, i])
            names.append(fasta_lines[i][1:-1])
            start = i
    idx.append([start, len(fasta_lines)])
    return np.array(names), np.array(idx, "int")

def main():
    data_dir = "/afs/inf.ed.ac.uk/user/s13/s1333321/research/PTM-X/data/ortholog/"
    fasta_dir = "/homes/huangh/crosstalkPTM/data/ortholog/InParanoid.fasta/"
    hm_ms_rt_file = data_dir + "human_mouse_rat.txt"
    ortholog = np.loadtxt(hm_ms_rt_file, "str", delimiter="\t", skiprows=1)

    human_pro_file = data_dir + "9606.fasta"
    mouse_pro_file = data_dir + "10090.fasta"
    rat_pro_file = data_dir + "10116.fasta"

    fid = open(human_pro_file,"r")
    human_pro = fid.readlines()
    fid.close()
    
    fid = open(mouse_pro_file,"r")
    mouse_pro = fid.readlines()
    fid.close()

    fid = open(rat_pro_file,"r")
    rat_pro = fid.readlines()
    fid.close()

    hm_name, hm_idx = get_name_idx(human_pro)
    ms_name, ms_idx = get_name_idx(mouse_pro)
    rt_name, rt_idx = get_name_idx(rat_pro)

    for i in range(ortholog.shape[0]):
        _idx1 = np.where(ortholog[i,0]==hm_name)[0]
        _idx2 = np.where(ortholog[i,1]==ms_name)[0]
        _idx3 = np.where(ortholog[i,2]==rt_name)[0]
        if (_idx1.shape[0] * _idx2.shape[0] * _idx3.shape[0]) == 0:
            print("No protein for all species!")
            continue
        hm_seq = human_pro[hm_idx[_idx1[0],0] : hm_idx[_idx1[0],1]]
        ms_seq = mouse_pro[ms_idx[_idx2[0],0] : ms_idx[_idx2[0],1]]
        rt_seq = rat_pro[rt_idx[_idx3[0],0] : rt_idx[_idx3[0],1]]

        fid = open(data_dir + "InParanoid.fasta/" + ortholog[i,0]+".fasta","w")
        fid.writelines(hm_seq)
        fid.writelines(ms_seq)
        fid.writelines(rt_seq)


if __name__ == '__main__':
    main()
    