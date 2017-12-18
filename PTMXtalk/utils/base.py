# This file contains some basic object and function as common utilities

import numpy as np

class FastaMSA:
    """This class is to load and to handle fasta file, which
    must be small size, e.g. less than 100M in plain text.
    This because it load the full data; big fasta file will
    take a lot of memenory and a long time to load. For big
    fasta files, you could use pysam.FastaFile. """
    def __init__(self, fasta_file):
        fid = open(fasta_file, "r")
        all_lines = fid.readlines()
        fid.close()
        seq, self.ref, self.seq = "", [], []
        for line in all_lines:
            line = line.split("\n")[0]
            if len(self.ref) == 0 and line[0] != ">": continue
            if line.startswith(">"):
                self.ref.append(line[1:])
                if seq == "": continue
                self.seq.append(seq)
                seq = ""
            else:
                seq = seq + line
        self.seq.append(seq)

    def get_seq(self, qref, start=None, stop=None):
        """get the sequence in a given region, the start is from 1.
        You will get stop-start+1 chars."""
        try:
            idx = self.ref.index(qref)
        except ValueError:
            print "No reference id as the query: %s" %qref
            return None
        if start is None:
            return self.seq[idx]
        try:
            RV = self.seq[idx][start-1 : stop]
        except ValueError:
            print "Wrong start or stop position: %d, %d" %(start, stop)
            return None
        return RV

    def remove_paralog(self, qref, delimiter=".", organism_idx=0):
        """In multiple sequence alignment, multiple paralogs may exist in a 
        file. We could only keep one paralog that has highest consistence with 
        the reference item.
        """
        try:
            idx = self.ref.index(qref)
        except ValueError:
            print "No reference id as the query: %s" %qref
            return None
        qref_seq = []
        qref_seq[:] = self.seq[idx]

        idx_kept = []
        organisms = [x.split(delimiter)[organism_idx] for x in self.ref]
        for _org in np.unique(organisms):
            org_idx = np.where(organisms == _org)[0]
            if len(org_idx) == 1:
                idx_kept.append(org_idx[0])
                continue
            _seq = []
            for ii in range(len(org_idx)):
                _seq[:] = self.seq[ii]
                similarity = np.mean(np.array(qref_seq) == np.array(_seq))
            idx_kept.append(org_idx[np.argmax(similarity)])

        self.ref = [self.ref[i] for i in idx_kept]
        self.seq = [self.seq[i] for i in idx_kept]

    def get_msa_loc(self, qref, original_loc, residue, shift_loc=[-1,1,-2,2], 
        inserts="-", verbose=True):
        """Get the locus in the multiple sequence alignment given the ref name 
        and original locus. Check: inserts="-?X"
        Note, msa_loc starts from 1 rather than 0.
        """
        cnt = -1
        msa_loc = 0
        seq = self.get_seq(qref)
        if seq is None:
            return None
        for i in range(len(seq)):
            if inserts.count(seq[i]) == 0:
                cnt += 1
            if cnt == original_loc:
                msa_loc = i
                break
        if msa_loc == -1:
            if verbose:
                print("%s is too short for %s%d." %(qref, residue, 
                    original_loc))
        elif seq[msa_loc-1] != residue:
            #check locus shift
            for ii in shift_loc:
                near_loc = msa_loc+ii
                if near_loc < 0 or near_loc > len(seq):
                    pass
                elif seq[near_loc-1] == residue:
                    msa_loc = near_loc
                    break
            if verbose and seq[msa_loc-1] != residue:
                idx_min = max(msa_loc-4, 0)
                idx_max = min(msa_loc+5, len(seq))
                ref_motif = seq[idx_min : idx_max]
                print("%s %s%d doesn't match %s" %(qref, residue, original_loc, 
                    ref_motif))
                msa_loc = -1
        return msa_loc


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


def permutation_test(x1, x2, alternative='two-sided', times=1000):
    """
    Classical permutation test implementation, which alows us to test the 
    significant of the distances in two sets of observations. Missing value 
    will omitted.
    
    Parameters
    ----------
    x1, x2 : sequence of 1-D ndarrays
        two arrays of sample observations assumed to be drawn from a continuous 
        distribution, sample sizes can be different
    alternative : {'two-sided', 'less', 'greater'}, optional
        Which alternative hypothesis to the null hypothesis the test uses. 
        Default is 'two-sided'.
    times : int, optional
        Number of permutation times

    Returns
    -------
    p_value : float
        P-value, the probability of obtaining a mean difference after 
        permututation greater than the observed, assuming that the null 
        hypothesis is true.
    """
    x1 = x1[x1==x1]
    x2 = x2[x2==x2]

    L1 = len(x1)
    xx = np.append(x1, x2)
    d0 = np.mean(x1) - np.mean(x2)
    
    N_success = 0
    for i in range(times):
        yy = np.random.permutation(xx)
        d1 = np.mean(yy[:L1]) - np.mean(yy[L1:])
        if alternative == "two-sided" and abs(d1) >= abs(d0):
            N_success += 1
        if alternative == "greater" and d1 >= d0:
            N_success += 1
        if alternative == "less" and d1 <= d0:
            N_success += 1
    p_value = N_success / (times+0.0)
    return p_value 
