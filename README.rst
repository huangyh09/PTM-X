PTM-X: PTM cross-talk predictor
===============================

About PTM-X
-----------

PTM-X is an algorithm to predict posttranslational modification(PTM) cross-talk 
both intra- and inter proteins. In current version, the following four features 
are used to predict PTM cross-talk inter proteins, via a 
`random forest <http://scikit-learn.org/stable/modules/ensemble.html#forest>`_ 
classifier.

1. ``sequence co-evolution``: The coevolution between two the amino acids across 
multiple vertebrates. The multiple sequence alignment data is downloaded from 
`eggNOG v4.5 <http://eggnogdb.embl.de>`_. Here, the 
1-`Hamming loss <http://scikit-learn.org/stable/modules/model_evaluation.html#hamming-loss>`_ 
is used to measure the co-evolution score.

2. ``motif co-conservation``: The +/-4 surrounding amino acids are used to form 
the motif for a given PTM locus. The fraction of consistence the this motif in 
any vertebrate comparing to its human orthologous protein. Based on the two 
vectors of motif conservation fraction, the co-conservation score is calculated 
by the mean product of the two vectors. Again, the eggNOG v4.5 is used here.

3. ``PTM co-conservation``: The co-conservation of PTM existence across human, 
mouse and rat. The raw PTM data for these three species is downloaded from 
`PhosphoSitePlus <https://www.phosphosite.org>`_.

4. ``PTM co-occurrence``: The co-occurrence between the two PTMs across across 
88 tissue, disease and cellline conditions.

Get started
-----------
Note, you need download the data to fetch above four features. Downloaded the 
data from `here <http://ufpr.dl.sourceforge.net/project/ptm-crosstalk/PTM-X_data_v2.2.zip>`_,
unzip it, and set the directory as the accordding parameter.

Installation
~~~~~~~~~~~~
**Installation**
Download the codes from this github repository [`master.zip <https://github.com/huangyh09/PTM-X/archive/master.zip>`_] 
and then run the following command line:

::

    python setup.py install

If you don't have the root permission, add ``--user``.

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
* demo file: https://github.com/huangyh09/PTM-X/blob/master/demo.sh
* web server: http://bioinfo.bjmu.edu.cn/ptm-x/
* data repository: http://ptm-crosstalk.sourceforge.net


Reference
---------
* For version 1: Huang et al.: `Systematic Characterization and Prediction of Post-Translational Modification Cross-talk <http://www.mcponline.org/content/14/3/761>`_. **Molecular & Cellular Proteomics**, 2015, 14 (3), 761-770.

* For version 2: comming soon.
