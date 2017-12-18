#!/usr/bin/python2.7
# merge the human-mouse and human-rat orthoglo from InParaniod

import sys
sys.path.append('../PTMXtalk/')
from utils.base import id_mapping

import numpy as np
from os.path import expanduser

def load_inParanoid(data_file, threshold=90, ref_org="H.sapiens", 
    multiple_ref=True):
    """
    When multiple ref proteins map to multiple proteins in the other species,
    Keep all ref proteins but only use the highest or the first protein in the 
    other species.
    """
    fid = open(data_file, "r")
    all_lines = fid.readlines()
    fid.close()

    data = []
    for _line in all_lines:
        if len(_line.split()) == 6:
            data.append(_line.split())
    data = np.array(data)

    idx_sort = np.argsort(data[:, 0])
    data_use = data[idx_sort, :]
    data_ids, data_idx = np.unique(data_use[:, 0], return_index=True)

    pair_ids = []
    for i in range(data_ids.shape[0]):
        if i == (data_ids.shape[0] - 1):
            idx_item = np.arange(data_idx[i], data_use.shape[0])
        else:
            idx_item = np.arange(data_idx[i], data_idx[i+1])

        organism_unique = np.unique(data_use[idx_item, 2])
        if len(organism_unique) < 2 or (ref_org not in organism_unique):
            continue

        scores = np.array([float(x[:-1]) for x in data_use[idx_item, 5]])
        _idx1 = np.where(data_use[idx_item, 2] == ref_org)[0]
        _idx2 = np.where(data_use[idx_item, 2] != ref_org)[0]
        idx_1 = idx_item[_idx1[np.argmax(scores[_idx1])]]
        idx_2 = idx_item[_idx2[np.argmax(scores[_idx2])]]
        if min(np.max(scores[_idx1]), np.max(scores[_idx2])) <= threshold:
            continue
        if multiple_ref == False:
            _idx1 = np.array([_idx1[np.argmax(scores[_idx1])]])
        for idx_1 in idx_item[_idx1]:
            pair_ids.append([data_use[idx_1, 4], data_use[idx_2, 4]])
    return np.array(pair_ids)


def main():
    data_dir = expanduser("~") + "/research/PTM-X/data/ortholog/InParanoid.raw/"
    hm_ms_file = data_dir + "sqltable.H.sapiens-M.musculus"
    hm_rt_file = data_dir + "sqltable.H.sapiens-R.norvegicus"
    hm_ms_rt_file = data_dir + "/../human_mouse_rat.txt"

    threshold = 90

    hm_ms_map = load_inParanoid(hm_ms_file, threshold)
    hm_rt_map = load_inParanoid(hm_rt_file, threshold)

    idx0 = id_mapping(hm_ms_map[:,0], hm_rt_map[:,0])
    idx1 = np.arange(len(idx0))[idx0>=0]
    idx2 = idx0[idx1].astype(int)

    hm_ms_map_out = hm_ms_map[idx1, :]
    hm_rt_map_out = hm_rt_map[idx2, :]
    # print(hm_ms_map_out)
    # print(hm_rt_map_out)

    fid = open(hm_ms_rt_file,"w")
    fid.writelines("human\tmouse\trat\n")
    for i in range(hm_ms_map_out.shape[0]):
        _list = np.append(hm_ms_map_out[i,:], hm_rt_map_out[i,1])
        fid.writelines("\t".join(list(_list)) + "\n")
    fid.close()


if __name__ == '__main__':
    main()

