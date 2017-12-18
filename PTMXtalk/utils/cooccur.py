
import numpy as np
import multiprocessing
from scipy.stats import fisher_exact

def cooccur_test(raw_lines, pro_PTM1, pro_PTM2, min_line=1):
    """raw_lines list all conditions and the PTMs occurring in that condition.
    This function will count the four situations for a pair of PTM:
    1) num_00: pro_PTM1 No and pro_PTM2 No
    2) num_01: pro_PTM1 No and pro_PTM2 Yes
    3) num_10: pro_PTM1 Yes and pro_PTM2 No
    4) num_11: pro_PTM1 Yes and pro_PTM2 Yes
    and return the Fisher exact test results. 
    Null hypothesis: pro_PTM2 is not more likely to occur when pro_PTM1 occurs.

    Note, pro_PTM1 and pro_PTM2 format: A0PJX4+S162
    Adding "+" is in find important to distinguish A0PJX4+S16+ and A0PJX4+S162+
    """
    num_lines = len(raw_lines)
    num_00, num_01, num_10, num_11 = 0, 0, 0, 0
    for _line in raw_lines:
        FIND1 = _line.find(pro_PTM1 + "+")
        FIND2 = _line.find(pro_PTM2 + "+")
        if FIND1 == -1:
            if FIND2 == -1:
                num_00 += 1
            else:
                num_01 += 1
        else:
            if FIND2 == -1:
                num_10 += 1
            else:
                num_11 += 1
    num_PTM1 = num_10 + num_11
    num_PTM2 = num_01 + num_11
    if num_PTM1 < min_line or  num_PTM2 < min_line:
        odd_ratio, pval = None, None
    else:
        odd_ratio, pval = fisher_exact([[num_00, num_01], [num_10, num_11]], 
            alternative='greater')
        pval = -np.log10(pval)
    # print([[num_00, num_01], [num_10, num_11]], num_PTM1, num_PTM2, pval)
    return odd_ratio, pval, min(num_PTM1, num_PTM2)


def fetch_coOccur(samples, cooccur_file, verbose=False, nproc=1):
    f1 = open(cooccur_file, "r")
    occur_lines = f1.readlines()
    f1.close()

    cooccur_val = np.zeros((samples.shape[0], 3)) #odd_ratio, log10_p, min_occur
    if nproc <= 1:
        for i in range(samples.shape[0]):
            pro_PTM1 = "+".join(list(samples[i, 0:2]))
            pro_PTM2 = "+".join(list(samples[i, 3:5]))
            cooccur_val[i, :] = cooccur_test(occur_lines, pro_PTM1, pro_PTM2)
    else:
        pool = multiprocessing.Pool(processes=nproc)
        result = []
        for i in range(samples.shape[0]):
            pro_PTM1 = "+".join(list(samples[i, 0:2]))
            pro_PTM2 = "+".join(list(samples[i, 3:5]))
            result.append(pool.apply_async(cooccur_test, 
                (occur_lines, pro_PTM1, pro_PTM2)))
        pool.close()
        pool.join()
        cooccur_val = np.array([res.get() for res in result], dtype=float)

    idx = cooccur_val[:,1] == cooccur_val[:,1]
    mean_val = np.mean(cooccur_val[idx,1])
    print("[PTM-X] fetched PTM co-occurrence for %d samples. mean: " 
          "%.3f, nan: %.1f%%." %(len(samples), mean_val, 100-np.mean(idx)*100))
    return cooccur_val

