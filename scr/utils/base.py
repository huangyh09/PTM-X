# This function file contains several fucntions that can be used widely,
# especially for analyses.
# For use: from func_common import *

import numpy as np

def get_same_distribute(x0, x1, is_log_scale=True):
    """
    This function is generate an identical subset from x1 to x0
    with the sample size as large as possible.
    It supports log scale and alos orignal scale. For the latter,
    the sets x0 ana x1 is required to a suitable scale, it treats
    them in integral format.
    """
    if is_log_scale:
        x0 = np.log(x0).astype("int")
        x1 = np.log(x1).astype("int")
    else:
        x0 = x0.astype("int")
        x1 = x1.astype("int")
    
    x0_unq = np.unique(x0)
    x1_unq = np.unique(x1)
    
    num0 = np.zeros(x0_unq.shape[0]) # for set x0
    num1 = np.zeros(x0_unq.shape[0]) # for set x1
    for i in range(x0_unq.shape[0]):
        _idx0 = np.where(x0_unq[i] == x0)[0]
        _idx1 = np.where(x1_unq[i] == x1)[0]
        num0[i] = _idx0.shape[0]
        num1[i] = _idx1.shape[0]
    min_time = int(min(num1/num0))   # times of the new set from set x0
    
    rv_idx0 = np.array([])
    rv_idx1 = np.array([])
    for i in range(x0_unq.shape[0]):
        _idx0 = np.where(x0_unq[i] == x0)[0]
        _idx1 = np.where(x1_unq[i] == x1)[0]

        rv_idx0 = np.append(rv_idx0, _idx0)
        rv_idx1 = np.append(rv_idx1, np.random.permutation(_idx1)[:min_time*_idx0.shape[0]])
    return rv_idx0.astype("int"), rv_idx1.astype("int")


def permutation_test(feature1, feature2, times=1000):
    """
    This functin is a classical permutation test implementation,
    which alows us to test the significant of the distances in 
    two features.
    It can treat missing features automatically.
    But please bear in mind that it can take a long time when "tiems" is greater than 10^5.
    """
    idx1 = (feature1==feature1)
    idx2 = (feature2==feature2)
    feature1, feature2 = feature1[idx1], feature2[idx2]
    feature_both = np.append(feature1, feature2)
    deta = np.mean(feature1) -  np.mean(feature2)

    cnt = 0
    for i in range(times):
        np.random.shuffle(feature_both)
        deta_tmp = np.mean(feature_both[:sum(idx1)]) - np.mean(feature_both[sum(idx1):])
        if deta_tmp >= deta:
            cnt = cnt + 1
    p_value = cnt / (times+0.0) # this is the p-value of the permutation test
    return p_value 