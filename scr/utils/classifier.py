# This is file contains an object for cross-validation and an object of naive 
# Bayes with Gaussian kernel density estimate.

# To import this: 
# import sys
# sys.path.append('../scr')
# from classifier import KernalNB, CrossValidation

import numpy as np
import scipy.stats as st

class KernalNB:
    def __init__(self, priors=None):
        self.priors = priors

    def fit(self, X, Y):
        self.C = np.sort(np.unique(Y))
        self.density_all = []
        for i in range(len(self.C)):
            _density = []
            for j in range(X.shape[1]):
                idx = (Y == self.C[i]) * (X[:,j] == X[:,j])
                _density.append(st.gaussian_kde(X[idx, j]))
            self.density_all.append(_density)

    def predict(self, Xtest):
        proba = self.predict_proba(Xtest)
        state = np.argmax(proba, axis=1)
        return state

    def predict_proba(self, Xtest):
        proba = np.ones((Xtest.shape[0], len(self.C)))

        for j in range(Xtest.shape[1]):
            idx = Xtest[:,j] == Xtest[:,j]
            for i in range(len(self.C)):
                proba[idx,i] *= self.density_all[i][j](Xtest[idx,j])
        proba = proba / proba.sum(axis=1).reshape(-1,1)
        return proba


class CrossValidation:
    def __init__(self, X, Y):
        """General cross-validation for both classification and regression.
        1) check the shape of X and Y are compatible with the input model;
        2) the model doesn't have memory on previous trainings, i.e., all 
        parameters are only based on current fit.
        """
        self.X, self.Y = X, Y
        
    def cv_regression(self, model, folds=3, shuffle=True):
        """run cross-validation for regression.
        For regression, make sure the input model has the following 
        functions: fit and predict
        """
        self.Y_pre = np.zeros(self.Y.shape[0])
        fold_len = int(self.Y.shape[0] / folds)
        idx_all = np.arange(self.Y.shape[0])
        if shuffle:
            np.random.shuffle(idx_all)
            
        for i in range(folds):
            if i < folds - 1:
                _idx = idx_all[i*fold_len : (i+1)*fold_len]
            else:
                _idx = idx_all[i*fold_len : self.Y.shape[0]]
            Xtest = self.X[_idx, :]
            Xtrain = np.delete(self.X, _idx, 0)
            Ytrain = np.delete(self.Y, _idx)
            
            model.fit(Xtrain, Ytrain)
            self.Y_pre[_idx] = model.predict(Xtest)
        return self.Y_pre
    
    def cv_classification(self, model, folds=3, shuffle=True):
        """Run cross-validation for classification.
        For classification, make sure the input model has the following 
        functions: fit, predict and predict_proba.
        """
        cc = np.unique(self.Y)
        self.Ystate = np.zeros(self.Y.shape[0])
        self.Yscore = np.zeros((self.Y.shape[0], len(cc)))
        idx_all = []
        fd_lens = []
        for i in range(len(cc)):
            _idx = np.where(self.Y == cc[i])[0]
            if shuffle: np.random.shuffle(_idx)
            idx_all.append(_idx)
            fd_lens.append(int(len(_idx)/folds))
        
        for i in range(folds):
            idx_use = np.array([], "int")
            for j in range(len(cc)):
                if i < folds-1:
                    _idx = idx_all[j][i*fd_lens[j]: (i+1)*fd_lens[j]]
                else:
                    _idx = idx_all[j][i*fd_lens[j]:]
                idx_use = np.append(idx_use, _idx)
                    
            Xtest = self.X[idx_use, :]
            Xtrain = np.delete(self.X, idx_use, 0)
            Ytrain = np.delete(self.Y, idx_use)
            
            model.fit(Xtrain, Ytrain)
            self.Ystate[idx_use] = model.predict(Xtest)
            model.fit(Xtrain, Ytrain)
            self.Yscore[idx_use,:] = model.predict_proba(Xtest)
        return self.Ystate, self.Yscore

