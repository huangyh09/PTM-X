# This is file contains an object for cross-validation and an object of naive 
# Bayes with Gaussian kernel density estimate.

# To import this: 
# import sys
# sys.path.append('../PTMXtalk')
# from classifier import KernalNB, CrossValidation

import numpy as np
import scipy.stats as st
import multiprocessing

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


def lite_fit(model, X, y):
    model.fit(X, y)
    return model

def lite_predict(model, X):
    return model.predict(X)

def lite_predict_proba(model, X):
    return model.predict_proba(X)


class MultiModel:
    """Bootstraping samples in each class with option to generate balanced 
    subsamples. On each bootstrap sample set, a separate model will be trained. 
    """
    def __init__(self, model, n_model=1, bootstrap=True, balanced=True, 
        n_sample=None, n_jobs=1):
        """Initialize the MultiModel object

        Parameters:
        -----------
        mdoel: classifier object
            It should have these functions: fit, predict, predict_proba
        n_model: int
            The number of models to fit, same as bootstrap sets
        bootstrap: bool
            If True, do bootstrap otherwise no repeat samples
        balanced: bool
            If True, generate the same size for each class
        n_sample: int
            The number of samples in each class, shoud be smaller than min class
        n_jobs: int
            Number of cores for parallel fit and predict on multiple models. 
            When it is -1, it will use half of the total cores.
        """
        self.models = [model] * n_model
        self.n_jobs = n_jobs
        self.n_model = n_model
        self.n_sample = n_sample
        self.balanced = balanced
        self.bootstrap = bootstrap
        
    def fit(self, X, Y, n_jobs=None):
        """If n_jobs is None, than use the initalized n_jobs.
        """
        y_uniq, y_len = np.unique(Y, return_counts=True)
        xx_set = [X[Y==Ytmp, :] for Ytmp in y_uniq]
        if self.balanced:
            if self.n_sample is None or self.n_sample > np.min(y_len):
                y_len[:] = np.min(y_len)
            else:
                y_len[:] = self.n_sample
            
        X_list, Y_list = [], []
        for i in range(self.n_model):
            xx, yy = np.zeros((0, X.shape[1])), np.zeros(0)
            for j in range(len(y_uniq)):
                if self.bootstrap:
                    idx = np.random.choice(sum(Y==y_uniq[j]), y_len[j])
                else:
                    remainder = i % sum(Y==y_uniq[j])
                    idx = np.arange(y_len[j]) + remainder
                yy = np.append(yy, [y_uniq[j]] * len(idx))
                xx = np.append(xx, xx_set[j][idx, :], axis=0)
            X_list.append(xx)
            Y_list.append(yy)

        if n_jobs is None:
            n_jobs = self.n_jobs
        if n_jobs == -1:
            n_jobs = int(multiprocessing.cpu_count() / 2)
        if n_jobs > 1:
            pool = multiprocessing.Pool(processes=n_jobs)
            result = []
            for i in range(self.n_model):
                result.append(pool.apply_async(lite_fit, (self.models[i], 
                    X_list[i], Y_list[i])))
            pool.close()
            pool.join()
            self.models = [res.get() for res in result]
        else:
            for i in range(self.n_model):
                self.models[i].fit(X_list[i], Y_list[i])

            
    def predict(self, Xtest, merge=True, n_jobs=None):
        """If n_jobs is None, than use the initalized n_jobs.
        """
        if n_jobs is None:
            n_jobs = self.n_jobs
        if n_jobs == -1:
            n_jobs = int(multiprocessing.cpu_count() / 2)
        if n_jobs > 1:
            pool = multiprocessing.Pool(processes=n_jobs)
            result = []
            for _model in self.models:
                result.append(pool.apply_async(lite_predict, (_model, Xtest)))
            pool.close()
            pool.join()
            states_list = [res.get() for res in result]
        else:
            states_list = []
            for _model in self.models:
                states_list.append(_model.predict(Xtest))
            
        if merge:
            RT_states = np.zeros(len(states_list[0]))
            for j in range(len(RT_states)):
                _states = [s[j] for s in states_list]
                (values, counts) = np.unique(_states, return_counts=True)
                RT_states[j] = values[np.argmax(counts)]
        else:
            RT_states = states_list

        return RT_states
    
    def predict_proba(self, Xtest, merge=True, n_jobs=None):
        """If n_jobs is None, than use the initalized n_jobs.
        """
        if n_jobs is None:
            n_jobs = self.n_jobs
        if n_jobs == -1:
            n_jobs = int(multiprocessing.cpu_count() / 2)
        if n_jobs <= 1:
            pool = multiprocessing.Pool(processes=n_jobs)
            result = []
            for _model in self.models:
                result.append(pool.apply_async(lite_predict_proba, 
                    (_model, Xtest)))
            pool.close()
            pool.join()
            scores_list = [res.get() for res in result]
        else:
            scores_list = []
            for _model in self.models:
                scores_list.append(_model.predict_proba(Xtest))

        if merge:
            RT_scores = scores_list[0]
            for j in range(1, len(scores_list)):
                RT_scores += scores_list[j]
            RT_scores = RT_scores / len(scores_list)
        else:
            RT_scores = scores_list

        return RT_scores


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

