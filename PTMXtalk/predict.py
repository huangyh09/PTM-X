#!/usr/bin/python2.7

# This file is to predict crosstalk of PTMs across proteins via a random
# forest classifier.

import os
import numpy as np
from optparse import OptionParser
from utils.classifier import MultiModel
from sklearn.ensemble import RandomForestClassifier

def main():
    #0. parse command line options
    parser = OptionParser()
    
    parser.add_option("-i", "--input-file", dest="test_file", 
        help="The feature file to predict.")
    parser.add_option("-o", "--output-file", dest="predict_file", 
        help="The output file with predicted results", default=None)
    parser.add_option("--positive", dest="positive_file", 
        help="The feature file for positve samples", default=None)
    parser.add_option("--negative", dest="negative_file", 
        help="The feature file for negative samples", default=None)
    parser.add_option("-s", "--seed", dest="seed", default=None, type="float",
        help="Seed for randomness. default not set.")
    
    # main arguments
    (options, args) = parser.parse_args()
    test_file = options.test_file
    predict_file = options.predict_file
    positive_file = options.positive_file
    negative_file = options.negative_file

    if options.seed is not None:
        np.random.seed(int(options.seed))
    if os.access(os.path.dirname(predict_file), os.W_OK) == False:
        print("[PTM-X] Error: output-file is not writable. Please check!")
        exit()

    X1 = np.genfromtxt(positive_file, delimiter='\t', skip_header=1, 
        dtype="str")[:, 6:].astype(float)
    X2 = np.genfromtxt(negative_file, delimiter='\t', skip_header=1, 
        dtype="str")[:, 6:].astype(float)
    Xtest = np.genfromtxt(test_file, delimiter='\t', skip_header=1, 
        dtype="str")[:, 6:].astype(float)

    X = np.append(X1, X2, axis=0).astype(float)
    Y = np.append(np.ones(X1.shape[0]), np.zeros(X2.shape[0]))

    scores = np.zeros(Xtest.shape[0])
    states = np.zeros(Xtest.shape[0])
    scores[:], states[:] = None, None

    att_code_test = np.dot(Xtest==Xtest, 2**np.arange(4))
    att_code_train = np.dot(X==X, 2**np.arange(4))

    att_use = [[0,1,2,3], [0,1,2], [0,1]]
    att_code_use = [15, 7, 3]
    for i in range(len(att_use)):
        idx_test = att_code_test == att_code_use[i]
        idx_train = np.mean(X[:,att_use[i]] == X[:,att_use[i]], axis=1) == 1

        _Y = Y[idx_train]
        _X = X[idx_train, :][:, att_use[i]]
        _Xtest = Xtest[idx_test, :][:, att_use[i]]

        if sum(idx_test) == 0 or len(np.unique(_Y)) < 2:
            continue

        # RF_model = RandomForestClassifier(n_estimators=100, n_jobs=-1)
        # RF_model.fit(_X, _Y)

        # states[idx_test] = RF_model.predict(_Xtest)
        # scores[idx_test] = RF_model.predict_proba(_Xtest)[:,1]

        RF_model = RandomForestClassifier(n_estimators=100, n_jobs=-1)
        mModel = MultiModel(model=RF_model, n_model=10)
        mModel.fit(_X, _Y)

        states[idx_test] = mModel.predict(_Xtest)
        scores[idx_test] = mModel.predict_proba(_Xtest)[:,1]

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
