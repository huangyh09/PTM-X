
# Analysis of disordered regions on PTM cross-talk

## Load functions and data


```python
import sys
sys.path.append('../')

from scr import fetch_disorder, count_disorder, load_disorder_file
```


```python
import numpy as np
# data_dir = "/Users/huangh/research/PTM/PTM-X_V2.1/"
data_dir = "/afs/inf.ed.ac.uk/user/s13/s1333321/research/PTM-X/"

sample_file = data_dir + "/interface/crosstalkSamples/positive_samples.txt"
samples_positive = np.loadtxt(sample_file, delimiter='\t', skiprows=0, dtype="str")

sample_file = data_dir + "/interface/crosstalkSamples/negative_samples.txt"
samples_negative = np.loadtxt(sample_file, delimiter='\t', skiprows=0, dtype="str")
```


```python
dis_file = data_dir + "/data/disorder/result.txt"
dis_count_positive = fetch_disorder(samples_positive, dis_file, residue_check=True, verbose=False)
dis_count_negative = fetch_disorder(samples_negative, dis_file, residue_check=True, verbose=False)
```

    [PTM-X] fetched disordered regions on 95 samples.
    [PTM-X] fetched disordered regions on 6001 samples.



```python
dis_count_negative
```




    array([[ 1.,  1.],
           [ 1.,  1.],
           [ 1.,  1.],
           ..., 
           [ 1.,  1.],
           [ 1.,  1.],
           [ 2.,  1.]])



## Preference on disordered regions


```python
idx1 = dis_count_positive == dis_count_positive
print("%.1f%% missing values in positive samples caused by wrong PTM annotation." %(100-np.mean(idx1)*100))

idx2 = dis_count_negative == dis_count_negative
print("%.1f%% missing values in negative samples caused by wrong PTM annotation." %(100-np.mean(idx2)*100))

idx_pair1 = idx1[:,0]*idx1[:,1]
idx_pair2 = idx2[:,0]*idx2[:,1]
print("%d postive and %d negatvie PTM pairs have valid values." %(sum(idx_pair1), sum(idx_pair2)))
```

    1.1% missing values in positive samples caused by wrong PTM annotation.
    0.2% missing values in negative samples caused by wrong PTM annotation.
    93 postive and 5974 negatvie PTM pairs have valid values.



```python
from scipy.stats import fisher_exact

for threshold in [1,2,3]:
    idx_dis_single1 = dis_count_positive[idx1] >= threshold
    idx_dis_single2 = dis_count_negative[idx2] >= threshold
    
    table = [[sum(idx_dis_single1), len(idx_dis_single1)-sum(idx_dis_single1)],
             [sum(idx_dis_single2), len(idx_dis_single2)-sum(idx_dis_single2)]]
    p_val = fisher_exact(table)[1]    
    print("threshold=%d: %.3f positve and %.3f negative PTMs in disordered regions (p=%.3f)." 
          %(threshold, np.mean(idx_dis_single1), np.mean(idx_dis_single2), p_val))
print("")

for threshold in [1,2,3]:  
    idx_dis_both1 = np.sum(dis_count_positive[idx_pair1,:] >= threshold, axis=1) == 2
    idx_dis_both2 = np.sum(dis_count_negative[idx_pair2,:] >= threshold, axis=1) == 2
    
    table = [[sum(idx_dis_both1), len(idx_dis_both1)-sum(idx_dis_both1)],
             [sum(idx_dis_both2), len(idx_dis_both2)-sum(idx_dis_both2)]]
    p_val = fisher_exact(table)[1]  
    print("threshold=%d: %.3f positve and %.3f negative PTM pairs in disordered regions (p=%.3f)." 
          %(threshold, np.mean(idx_dis_both1), np.mean(idx_dis_both2), p_val))
```

    threshold=1: 0.910 positve and 0.852 negative PTMs in disordered regions (p=0.029).
    threshold=2: 0.574 positve and 0.486 negative PTMs in disordered regions (p=0.018).
    threshold=3: 0.207 positve and 0.247 negative PTMs in disordered regions (p=0.233).
    
    threshold=1: 0.817 positve and 0.724 negative PTM pairs in disordered regions (p=0.047).
    threshold=2: 0.312 positve and 0.232 negative PTM pairs in disordered regions (p=0.083).
    threshold=3: 0.032 positve and 0.061 negative PTM pairs in disordered regions (p=0.376).



```python

```

### Check the output of  ``load_disorder_file``


```python
disorder_file = data_dir + "/data/disorder/result.txt"
dis_prot_ID, dis_regions, dis_seq = load_disorder_file(disorder_file)

print(len(dis_prot_ID), len(dis_regions), len(dis_seq))

print(dis_prot_ID[:5])
print(dis_regions[0][:3])
print(disorder_seq[:3])
```

    (19836, 19836, 19836)
    ['P31946', 'P62258', 'Q04917', 'P61981', 'P31947']
    [['105', '116'], ['133', '144'], ['159', '173']]
    ['mtmdkselvqkaklaeqaeryddmaaamkavteqghelsneernllsvayknvvgarrsswrvissieqkternekkqqmgkeyrekieaelqdicndvlelldKYLIPNATQPESkvfylkmkgdyfrylsEVASGDNKQTTVsnsqqayqeafeisKKEMQPTHPIRLGLAlnfsvfyyeilnspekacslaktafdeaiaeldtlneesykdstlimqllrdnlTLWTSENQGDEGDAGEGEN', 'mddredlvyqaklaeqaerydemvesmkkvagmdveltveernllsvayknvigarraswriissieqKEENKGGEDklkmireyrqmvetelkliccdildvLDKHLIPAANTGESkvfyykmkgdyhrylaefATGNDRKEaaenslvaykaasDIAMTELPPTHPIRLGLalnfsvfyyeilnspdracrlakaafddaiaeldtlseesykdstlimqllrdnltLWTSDMQGDGEEqnkealqdvedenq', 'mgdreqllqrarlaeqaeryddmasamkavteLNEPLSNEDrnllsvayknvvgarrsswrvissieqktmadgnekklekvkayrekiekeletvcndvlslldkflikncndfqyeskvfylkmkgdyyrylaevasgekknsvveaseaaykeafeisKEQMQPTHPIRLGLAlnfsvfyyeiqnapeqacllakqafddaiaeldTLNEDSYKdstlimqllrdnlTLWTSDQQDEEAGEGN']



```python

```
