import sys
import numpy as np
import h5py
from sklearn import preprocessing
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import math
from hmmlearn import hmm

def hmm_gaussian_use_mean_scores(mean_scores, start_prob = np.array([0.95,0.025,0.025]), transmat = 0, means = 0, covars = 0, chrunk_size = 500, n_components = 3, minproba = 0.05 , maxproba = 0.9):
    """ get motif use lsar scores
    Parameters
    -------
    mean_scores : list
        mean scores in ATCG minus ref scores after scale  (>=0)
    n_components : int
        The number of state
    transmat : array, shape (n_components, n_components)
        Matrix of transition probabilities between states. 
    means : list
        Mean and precision of the Normal prior distribtion
    covars : list
        gaussion distribution covariance(scores covariance)
    block_size : int
        use n bp block to analyse
    """
    if len(mean_scores) < chrunk_size:
        chrunk_size = len(mean_scores)
    if transmat == 0:
        transmat = np.array([[0.98, 0.01, 0.01],
                     [0.15, 0.84, 0.01],
                     [0.15, 0.01, 0.84]])
    mean_scores_list = list(divide_chunks(mean_scores, chrunk_size))
    chrunk_num = len(mean_scores_list)
    state_list=[]
    for i in range(chrunk_num):
        mean_scores_tmp = mean_scores_list[i]
        if means == 0:
            means_ = np.array([[0], [np.quantile(mean_scores_tmp, 0.05)], [np.quantile(mean_scores_tmp, 0.95)]])
        else:
            means_ = means
        if covars == 0:
            covars_ = np.array([[np.std(mean_scores_tmp)**2],[np.std(mean_scores_tmp)**2],[np.std(mean_scores_tmp)**2]])
        else:
            covars_ = covars
        model_guassion = hmm.GaussianHMM(n_components=n_components, covariance_type="diag")
        model_guassion.startprob_ = start_prob
        model_guassion.transmat_ = transmat
        model_guassion.means_ = np.array(means_)
        model_guassion.covars_ = np.array(covars_)
        seen = np.array( [[i] for i in mean_scores_tmp])
        EM_proba = model_guassion.predict_proba(seen)
        EM_proba[EM_proba > maxproba] = maxproba
        EM_proba[EM_proba < minproba] = minproba
        model_normal = hmm.MultinomialHMM(n_components=n_components)
        model_normal.startprob_ = start_prob
        model_normal.transmat_ = transmat
        model_normal.emissionprob_ = EM_proba.T
        logprob, state = model_normal.decode(np.array([range(len(seen))]).T, algorithm="viterbi")
        state_list.extend(state)
    return(state_list)
def divide_chunks(l, n):
    for i in range(0, len(l), n):
        if i + 2*n <= len(l):
            yield(l[i:i + n])
        elif i + n > len(l):
            continue
        else:
            yield(l[i:i+2*n])

sample = sys.argv[1]
h5file = "/public/home/jcli/data/orth_compare/Basenji/h5.list"
with open(h5file) as fh:
    for line in fh:
        line = line.strip()
    
        file = "/public/home/jcli/data/orth_compare/Basenji/MH63_score/" + line
        out = "/public/home/jcli/data/orth_compare/Basenji/motif_state/" + sample + "/" + line + ".csv"
        data = h5py.File(file,'r')

        select_index = int(sample)
        chrom = data['chrom'][0].decode()
        start = data['start'][0]
        end = data['end'][0]
        scores = data['scores'][0]
        scores_tmp = scores[:,:,select_index].T-scores[:,:,select_index][data['seqs'][0]]
        mean_scores = preprocessing.scale(pd.DataFrame(scores_tmp).mean().tolist())
        res = hmm_gaussian_use_mean_scores(mean_scores)

        pos = list(range(start,start+len(mean_scores)))
        df = pd.DataFrame(data={"chrom":chrom,"pos":pos,"score":mean_scores,"state":res})
        df.to_csv(out)
