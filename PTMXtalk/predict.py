#!/usr/bin/python2.7

# This file is to predict crosstalk of PTMs across proteins via a random
# forest classifier.

import os
import joblib
import numpy as np
from optparse import OptionParser, OptionGroup
from PTMXtalk.utils.classifier import MultiModel
from sklearn.ensemble import RandomForestClassifier


def main():
    #0. parse command line options
    parser = OptionParser()
    
    parser.add_option("-i", "--input-file", dest="test_file", 
        help="The feature file to predict.")
    parser.add_option("-o", "--output-file", dest="predict_file", 
        help="The output file with predicted results", default=None)
    parser.add_option("-m", "--model-file", dest="model_file", 
        help="The file with trained models. If lack, please input positive and "
        "negative feature files", default=None)

    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--positive", dest="positive_file", 
        help="The feature file for positve samples", default=None)
    group.add_option("--negative", dest="negative_file", 
        help="The feature file for negative samples", default=None)
    group.add_option("--model-out-file", dest="model_out", 
        help="The file for newly trained mdoels", default=None)
    group.add_option("--N-model", "-N", type="int", dest="n_model", 
        default="100", help="Number of models [default: %default]")
    group.add_option("-s", "--seed", dest="seed", default=None, type="float",
        help="Seed for randomness. default not set.")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default="1",
        help="Number of subprocesses to fit and predict [default: %default]")
    parser.add_option_group(group)

    # main arguments
    (options, args) = parser.parse_args()
    nproc = options.nproc
    test_file = options.test_file
    model_file = options.model_file
    predict_file = options.predict_file

    if options.seed is not None:
        np.random.seed(int(options.seed))
    if os.access(os.path.dirname(predict_file), os.W_OK) == False:
        print("[PTM-X] Error: output-file is not writable. Please check!")
        exit()

    # load model if exist otherwise fit model
    if model_file is not None:
        (RT_model, att_sets) = joblib.load(model_file)
    else:
        n_model = options.n_model
        model_out = options.model_out
        positive_file = options.positive_file
        negative_file = options.negative_file
        if model_out is None:
            model_out = os.path.dirname(predict_file) + "/trained_model_set.pkl"
        
        X1 = np.genfromtxt(positive_file, delimiter='\t', skip_header=1, 
            dtype="str")[:, 6:10].astype(float)
        X2 = np.genfromtxt(negative_file, delimiter='\t', skip_header=1, 
            dtype="str")[:, 6:10].astype(float)
        
        X = np.append(X1, X2, axis=0)
        Y = np.append(np.ones(X1.shape[0]), np.zeros(X2.shape[0]))

        RT_model = []
        #att_sets = [[0, 1, 2, 3], [0, 1, 2], [0, 1, 3], [0, 1]]
        att_sets = [[0, 1, 3], [0, 1]]
        for k in range(len(att_sets)):
            xx = X[:, att_sets[k]]
            ii = np.min(xx == xx, axis=1)
            RF_model = RandomForestClassifier(n_estimators=100, n_jobs=1)
            MM_model = MultiModel(model=RF_model, n_model=n_model, n_jobs=nproc)
            MM_model.fit(xx[ii, :], Y[ii])
            RT_model.append(MM_model)

        if os.access(os.path.dirname(model_out), os.W_OK) == False:
            print("[PTM-X] Warning: model-out-file is not writable!")
        else:
            joblib.dump((RT_model, att_sets), model_out)
    
    # predict test data
    Xtest = np.genfromtxt(test_file, delimiter='\t', skip_header=1, 
        dtype="str")[:, 6:10].astype(float)

    scores = np.zeros(Xtest.shape[0])
    states = np.zeros(Xtest.shape[0])
    scores[:], states[:] = None, None

    for k in range(len(att_sets)):
        xx = Xtest[:, att_sets[k]]
        ii = np.min(xx == xx, axis=1)
        ii = ii * (~(states == states)) # only unpredicted samples, still None
        if sum(ii) == 0:
            continue
        states[ii] = RT_model[k].predict(xx[ii,:], n_jobs=nproc)
        scores[ii] = RT_model[k].predict_proba(xx[ii,:], n_jobs=nproc)[:,1]

    # save results
    data = np.genfromtxt(test_file, delimiter='\t', skip_header=0, dtype="str")
    
    fid = open(predict_file, "w")
    head_line = "\t".join(list(data[0,:]) + ["predict_score", "predict_state"])
    fid.writelines(head_line + "\n")
    for i in range(len(scores)):
        val_out = list(data[i+1,:]) + ["%.3f" %scores[i],  "%.3f" %states[i]]
        fid.writelines("\t".join(val_out) + "\n")
    fid.close()

    print("[PTM-X] saved prediction results into file!")


if __name__ == '__main__':
    main()
