PTM-X: PTM cross-talk predictor
===============================

About PTM-X
-----------

PTM-X is an algorithm to predict post-translational modification(PTM) cross-talk 
both intra- and inter- protein. In current version, the following four features 
are used to predict PTM cross-talk inter proteins, via a 
`random forest <http://scikit-learn.org/stable/modules/ensemble.html#forest>`_ 
classifier.

1. ``sequence residue co-evolution``: The coevolution between two the amino acids across 
multiple vertebrates. The multiple sequence alignment data is downloaded from 
`eggNOG v4.5 <http://eggnogdb.embl.de>`_. Here, the 
1-`Hamming loss <http://scikit-learn.org/stable/modules/model_evaluation.html#hamming-loss>`_ 
is used to measure the co-evolution score.

2. ``sequence motif co-evolution``: The +/-3 surrounding amino acids are used to form 
the motif for a given PTM locus. The fraction of consistence the this motif in 
any vertebrate comparing to its human orthologous protein. Based on the two 
vectors of motif conservation fraction, the co-conservation score is calculated 
by the mean product of the two vectors. Again, the eggNOG v4.5 is used here.

3. ``Co-modification across different species``: The co-conservation of PTM existence across human, 
mouse and rat. The raw PTM data for these three species is downloaded from 
`PhosphoSitePlus <https://www.phosphosite.org>`_.

4. ``Co-modification across different conditions``: The co-occurrence between the two PTMs across 88 
tissue, disease and cellline conditions.

Get started
-----------


Installation
~~~~~~~~~~~~
**Step 1**: you need a Python environment for supporting packages. The easiest way 
might be installing the python platform via 
`Anaconda <https://www.anaconda.com/download/>`_. PTM-X is only compatible with 
Python 2.7. If you are using Anaconda 3, create a 
`conda environment <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_ 
with Python 2.7 and dependent packages as follows,

.. code:: bash

    conda create --name ptmxPy2 python=2.7 numpy scipy scikit-learn==0.17 joblib==0.11
    # activate the environment
    conda activate ptmxPy2

**Step 2**: install PTM-X package by either downloading the source code here or use the 
following command line for an equivalent setting:

.. code:: bash

    pip install -U git+https://github.com/huangyh09/PTM-X

In either case, if you don't have the root permission, add ``--user``.

**Step 3**, to fetch above four featuresd ownload the data from 
http://ufpr.dl.sourceforge.net/project/ptm-crosstalk/PTM-X_data_v2.2.zip
Unzip it and set the directory as the according parameter in the command lines below.


Fetch features
~~~~~~~~~~~~~~

* To see all arguments, try command line ``PTMX-feature -h``
* For the data_dir, set it as the downloaded data (see above).
* The python file is `feature_extractor.py <https://github.com/huangyh09/PTM-X/blob/master/PTMXtalk/feature_extractor.py>`_

Prediction
~~~~~~~~~~

* To see all arguments, try command line ``PTMX-predict -h``
* The python file is `predict.py <https://github.com/huangyh09/PTM-X/blob/master/PTMXtalk/predict.py>`_.

Links
-----
* **demo file**: https://github.com/huangyh09/PTM-X/blob/master/demo.sh
* web server: http://bioinfo.bjmu.edu.cn/ptm-x/
* data repository: http://ptm-crosstalk.sourceforge.net


References
----------
* For version 1 (intra-protein): Huang Y. et al.: `Systematic Characterization and Prediction of Post-Translational Modification Cross-talk <http://www.mcponline.org/content/14/3/761>`_. **Molecular & Cellular Proteomics**, 2015, 14 (3), 761-770.

* For version 2 (inter-protein): Huang R. et al.: `Systematic characterization and prediction of post-translational modification cross-talk between proteins <https://doi.org/10.1093/bioinformatics/bty1033>`_. **Bioinformatics**, 2019, 35 (15), 2626-2633.
