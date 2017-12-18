#!/usr/bin/python2.7

# This file is to extract the Ensembl id(s) in each file in the veNOG folder.
# Note, the any of the three ids can occur multiple times in the output file,
# to make the map file is complete.

import os
import numpy as np
from optparse import OptionParser

def get_veNOG_id(data_dir, org_id="9606"):
    """To get the Ensembl id for a given organism in all files in veNOG folder.
    """
    path_old = os.getcwd()
    os.chdir(data_dir)

    veNOG_Ensem = []
    file_cnt, pair_cnt = 0, 0
    for file_name in os.listdir(data_dir):
        file_cnt += 1
        veNOG_ID = file_name.split(".")[1]

        f_tmp = open(data_dir + "/" + file_name,"r")
        all_lines = f_tmp.readlines()
        f_tmp.close()

        # check the rule for display organism id
        for _line in all_lines:
            if _line.startswith(">" + org_id):
                pair_cnt += 1
                Ensembl_ID = _line.rstrip().split('.')[1]
                veNOG_Ensem.append([veNOG_ID, Ensembl_ID])
    veNOG_Ensem = np.array(veNOG_Ensem)

    os.chdir(path_old)
    print("[PTM-pre] %d files returning %d id pairs: %d veNOG and %d Ensembl." 
        %(file_cnt, pair_cnt, len(np.unique(veNOG_Ensem[:,0])), 
            len(np.unique(veNOG_Ensem[:,1]))))
    return veNOG_Ensem


def map_biomart_ids(biomart_file, veNOG_Ensem):
    """map Ensembl-Uniprot and veNOG-Ensembl ids.
    Note, the Ensembl id in veNOG-Ensembl pair need to be unique.
    """
    fid = open(biomart_file, "r")
    dat = fid.readlines()
    fid.close()

    eff_line = 0
    Ensem_Unip = []
    for i in range(len(dat)):
        if len(dat[i].split()) == 2:
            eff_line += 1
            Ensem_Unip.append(dat[i].split())
    Ensem_Unip = np.array(Ensem_Unip)
    Ensem_Unip = Ensem_Unip[np.argsort(Ensem_Unip[:,0]),:]
    print("[PTM-pre] %d out of %d lines are effective in biomart file." 
        %(eff_line, len(dat)))

    veNOG_Ensem = veNOG_Ensem[np.argsort(veNOG_Ensem[:,1]),:]
    if len(np.unique(veNOG_Ensem[:,1])) < len(veNOG_Ensem[:,1]):
        print("[PTM-pre] Warning: the Ensembl ID in veNOG-Ensembl map is not "
              "unique. Only the first one of the repeated id will be output. "
              "Please check!")
    
    i, j = 0, 0
    veNOG_Ensem_Unip = []
    while i < len(Ensem_Unip) and j < len(veNOG_Ensem):
        if Ensem_Unip[i,0] < veNOG_Ensem[j,1]:
            i += 1
        elif Ensem_Unip[i,0] == veNOG_Ensem[j,1]:
            veNOG_Ensem_Unip.append(np.append(veNOG_Ensem[j,:], Ensem_Unip[i,1]))
            i += 1
        elif Ensem_Unip[i,0] > veNOG_Ensem[j,1]:
            j += 1
    veNOG_Ensem_Unip = np.array(veNOG_Ensem_Unip)
            
    print("[PTM-pre] %d veNOG-Ensembl-Uniprot id pair returned: %d veNOG, %d "
          "Ensemble, %d Unirpot." %(veNOG_Ensem_Unip.shape[0], 
            len(np.unique(veNOG_Ensem_Unip[:,0])), 
            len(np.unique(veNOG_Ensem_Unip[:,1])),
            len(np.unique(veNOG_Ensem_Unip[:,2]))))
    return veNOG_Ensem_Unip

def lite_folder(data_dir, list_to_keep):
    for file_name in os.listdir(data_dir):
        veNOG_ID = file_name.split(".")[1]
        file_name.split(".")[1]
        if veNOG_ID not in list_to_keep:
            os.remove(data_dir + "/" + file_name)
    return None

def main():
    #0. parse command line options
    parser = OptionParser()
    parser.add_option("--veNOG-dir", "-d", dest="data_dir", default=None,
        help="The diroctory of the veNOG data")
    parser.add_option("--biomart-file", "-m", dest="biomart_file", default=None,
        help="The biomart file containing Enesmbl ID and Uniprot ID")
    parser.add_option("--out-file", "-o", dest="out_file", default=None,
        help="The file to save processed veNOG-Ensembl-Uniprot ID pairs")
    parser.add_option("--lite-folder", action="store_true", dest="remove_no_use", 
        default=False, help="remove no used files in veNOG folder.")
    
    (options, args) = parser.parse_args()
    data_dir = options.data_dir
    out_file = options.out_file
    biomart_file = options.biomart_file
    if out_file is None:
        out_file = data_dir + "/../veNOG_Ensem_Unip.txt"

    if os.path.isdir(data_dir) == False:
        print("[PTM-pre] Error: Data directory does not exist. Please check!")
        sys.exit()

    #save mapped veNOG, Ensembl ids
    fid = open(out_file, "w")
    veNOG_Ensem = get_veNOG_id(data_dir)
    veNOG_Ensem_Unip = map_biomart_ids(biomart_file, veNOG_Ensem)

    head_line = "veNOG_ID\tEnsem_ID\tUnip_ID\n"
    fid.writelines(head_line)
    for id_pair in veNOG_Ensem_Unip:
        fid.writelines("\t".join(list(id_pair)) + "\n")
    fid.close()

    #lite folder
    if options.remove_no_use:
        lite_folder(data_dir, veNOG_Ensem_Unip[:,0])


if __name__ == '__main__':
    main()

    