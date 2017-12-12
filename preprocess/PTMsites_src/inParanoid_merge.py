#!/usr/bin/python2.7
# run python inParanoid_merge.py

import numpy as np

def data_load(data_file):
    fid = open(data_file, "r")
    all_lines = fid.readlines()
    fid.close()
    
    # delete the orthogous protein with score lower than 90%
    lines_new = []
    for i in range(len(all_lines)):
        line_tmp = all_lines[i].split()
        if len(line_tmp) == 6 and int(line_tmp[5][:-1]) > 90:
            lines_new.append(line_tmp)
    return np.array(lines_new,"str")
    
def rm_no_ortholog(data_set, spc_num=2):
    # delete the protein with no orthologous protein
    del_idx = np.array([],"int")
    start_idx = 0
    tmp_num = data_set[0,0]
    for i in range(len(data_set)):
        if data_set[i,0] != tmp_num :
            # remove the protein with paralog
            if (i-start_idx != 2 or 
                np.unique(data_set[start_idx:i,2]).shape[0] < spc_num):
                del_idx = np.append(del_idx, np.arange(start_idx, i))
            start_idx = i
            tmp_num = data_set[i,0]
    RV = np.delete(data_set, del_idx, axis=0)
    return RV

def main():
    data_dir = "/afs/inf.ed.ac.uk/user/s13/s1333321/research/PTM-X/data/ortholog/"
    hm_ms_file = data_dir + "sqltable.H.sapiens-M.musculus"
    hm_rt_file = data_dir + "sqltable.H.sapiens-R.norvegicus"
    hm_ms_rt_file = data_dir + "human_mouse_rat.txt"

    # hm_ms_file = "/homes/huangh/crosstalkPTM/data/ortholog/sqltable.H.sapiens-M.musculus"
    # hm_rt_file = "/homes/huangh/crosstalkPTM/data/ortholog/sqltable.H.sapiens-R.norvegicus"
    # hm_ms_rt_file = "/homes/huangh/crosstalkPTM/data/ortholog/human_mouse_rat.txt"

    hm_ms_map = data_load(hm_ms_file)
    hm_rt_map = data_load(hm_rt_file)

    hm_ms_map = rm_no_ortholog(hm_ms_map)
    hm_rt_map = rm_no_ortholog(hm_rt_map)

    hm_idx1 = np.where(hm_ms_map[:,2] == "H.sapiens")[0]
    hm_idx2 = np.where(hm_rt_map[:,2] == "H.sapiens")[0]
    hm_pro1 = hm_ms_map[hm_idx1,4]
    hm_pro2 = hm_rt_map[hm_idx2,4]

    fid = open(hm_ms_rt_file,"w")
    fid.writelines("human\tmouse\trat\n")

    for i in range(hm_pro1.shape[0]):
        _idx = np.where(hm_pro1[i] == hm_pro2)[0]
        if _idx.shape[0] == 1:
            if (hm_ms_map[hm_idx1[i],0] != hm_ms_map[hm_idx1[i]+1,0] or 
                hm_rt_map[hm_idx2[_idx[0]],0] != hm_rt_map[hm_idx2[_idx[0]]+1,0]):
                print("something wrong!")
                exit()

            else :
                line_tmp = "\t".join([hm_ms_map[hm_idx1[i],4], 
                    hm_ms_map[hm_idx1[i]+1,4], 
                    hm_rt_map[hm_idx2[_idx[0]]+1,4]])
                fid.writelines(line_tmp + "\n")

    print(hm_ms_map.shape)
    print(hm_rt_map.shape)

    print(np.unique(hm_ms_map[:,0]).shape)
    print(np.unique(hm_rt_map[:,0]).shape)


if __name__ == '__main__':
    main()

