#!/usr/bin/python2.7

# This file is to merge the PTM existence of human, mouse and rat.
# Check: H3, P68431 

import sys
sys.path.append('../PTMXtalk/')
from utils.base import FastaMSA, id_mapping

import os
import numpy as np
from optparse import OptionParser


def collapse_PTMs(PTM_data):
    """collapse PTMs by pooling all modifications on a shared the residue.
    """
    collapsed_PTM_data = []

    # sort by uniprot id
    idx_sort = np.argsort(PTM_data[:, 2])
    PTM_full = PTM_data[idx_sort, :]
    prot_ids, prot_idx = np.unique(PTM_full[:, 2], return_index=True)
    for i in range(prot_ids.shape[0]):
        if i == (prot_ids.shape[0] - 1):
            PTM_prot = PTM_full[prot_idx[i]:]
        else:
            PTM_prot = PTM_full[prot_idx[i]:prot_idx[i+1]]

        # sort by residue loc
        site_loc = np.array([int(x[1:]) for x in PTM_prot[:, 4]])
        idx_sort = np.argsort(site_loc)
        site_loc = site_loc[idx_sort]
        PTM_prot = PTM_prot[idx_sort, :]
        site_ids, site_idx = np.unique(site_loc, return_index=True)
        for j in range(site_ids.shape[0]):
            if j == (site_ids.shape[0] - 1):
                PTM_site = PTM_prot[site_idx[j]:]
            else:
                PTM_site = PTM_prot[site_idx[j]:site_idx[j+1]]

            PTM_site[np.where(PTM_site == "")] = "0"
            PTM_pubm = PTM_site[:,6:9].astype(int)
            PTM_pubm = "|".join([str(x) for x in PTM_pubm.sum(axis=0)])
            PTM_type = "|".join(list(np.unique(PTM_site[:,3])))
            PTM_save = [PTM_site[0,2], PTM_site[0,4], PTM_type, 
                        PTM_site[0,5], PTM_pubm]
            collapsed_PTM_data.append(PTM_save)

    return np.array(collapsed_PTM_data)

def get_query(query_ids, ref_data):
    query_data = np.zeros((len(query_ids), ref_data.shape[1]), dtype="S100")
    query_data[:,:] = ["nan"] * ref_data.shape[1]

    idx_eff = np.where(query_ids != "nan")[0]
    query_data[idx_eff, :2] = [x.split(":") for x in query_ids[idx_eff]]

    ref_ids = [":".join(x[:2]) for x in ref_data]
    idx_map = id_mapping(query_ids[idx_eff], ref_ids)
    idx_use = idx_map >= 0
    idx_eff = idx_eff[idx_use]
    idx_map = idx_map[idx_use].astype(int)
    query_data[idx_eff, :] = ref_data[idx_map, :]
    return query_data

def main():
    #0. parse command line options
    parser = OptionParser()
    parser.add_option("-d", "--data_dir", dest="data_dir", 
        help="The directory of the data.")
    (options, args) = parser.parse_args()
    data_dir = options.data_dir
    out_file = data_dir + "/PTMsites/PTM_data_3species.txt"

    #1. orthologous data and ids
    msa_dir  = data_dir + "/ortholog/InParanoid.align/"
    map_file = data_dir + "/ortholog/human_mouse_rat.txt"
    orth_map = np.loadtxt(map_file, "str", delimiter="\t", skiprows=1)

    #2. load the PTM data
    human_PTM_file = data_dir + "/PTMsites/PTMsites_human.txt"
    mouse_PTM_file = data_dir + "/PTMsites/PTMsites_mouse.txt"
    rat_PTM_file   = data_dir + "/PTMsites/PTMsites_rat.txt"    
    human_PTM = np.loadtxt(human_PTM_file, "str", delimiter="\t", skiprows=1)
    mouse_PTM = np.loadtxt(mouse_PTM_file, "str", delimiter="\t", skiprows=1)
    rat_PTM   = np.loadtxt(rat_PTM_file, "str", delimiter="\t", skiprows=1)

    human_PTM = collapse_PTMs(human_PTM)
    mouse_PTM = collapse_PTMs(mouse_PTM)
    rat_PTM   = collapse_PTMs(rat_PTM)
    prot_ids, prot_idx = np.unique(human_PTM[:, 0], return_index=True)

    mouse_query, rat_query = [], []
    for i in range(len(prot_ids)):
        if i == (prot_ids.shape[0] - 1):
            idx_now = range(prot_idx[i], human_PTM.shape[0])
        else:
            idx_now = range(prot_idx[i], prot_idx[i+1])
        if os.path.isfile(msa_dir + "/%s.fa" %prot_ids[i]) == False:
            # print("No file of %s.fasta in the path %s." 
            #     %(prot_ids[i], msa_dir))
            mouse_query += ["nan"] * len(idx_now)
            rat_query   += ["nan"] * len(idx_now)
        else:
            MSA_File = FastaMSA(msa_dir + "/%s.fa" %prot_ids[i])
            for j in idx_now:
                _res = human_PTM[j, 1][0]
                _loc = int(human_PTM[j, 1][1:])
                human_loc = MSA_File.get_msa_loc(prot_ids[i], _loc, _res)

                #residue didn't match
                if human_loc is None or human_loc == -1:
                    mouse_query.append("nan")
                    rat_query.append("nan")
                else:
                    _seq = MSA_File.seq[1][:human_loc]
                    _res = MSA_File.seq[1][human_loc-1]
                    _loc = len(_seq) - _seq.count("-")
                    mouse_query.append(MSA_File.ref[1] + ":" + _res + str(_loc))

                    _seq = MSA_File.seq[2][:human_loc]
                    _res = MSA_File.seq[2][human_loc-1]
                    _loc = len(_seq) - _seq.count("-")
                    rat_query.append(MSA_File.ref[2] + ":" + _res + str(_loc))
    mouse_query, rat_query = np.array(mouse_query), np.array(rat_query)

    rat_mapped = get_query(rat_query, rat_PTM)
    mouse_mapped = get_query(mouse_query, mouse_PTM)

    # save mapped data
    fid = open(out_file, "w")
    head_list  = ["pro", "res", "PTM", "seq", "pub"]
    head_line  = [x+"_human" for x in head_list]
    head_line += [x+"_mouse" for x in head_list]
    head_line += [x+"_rat" for x in head_list]
    fid.writelines("\t".join(head_line) + "\n")

    for i in range(human_PTM.shape[0]):
        a_line = "\t".join(list(human_PTM[i,:]) + list(mouse_mapped[i,:]) + 
            list(rat_mapped[i,:])) + "\n"
        fid.writelines(a_line)
    fid.close()


if __name__ == '__main__':
    main()
    
    