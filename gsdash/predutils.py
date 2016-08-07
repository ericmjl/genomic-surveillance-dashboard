import numpy as np
from sklearn.externals import joblib

def load_model(drug):
    print('loading model for drug {0}'.format(drug))
    mdl = joblib.load("../models/base/{drug}/{drug}.pkl".format(drug=drug))
    return mdl, drug


def pred_range(model, datum):
    """
    Returns the full range of predictions for a given ensemble model.
    """
    preds = np.zeros(len(model.estimators_))
    for i, est in enumerate(model.estimators_):
        preds[i] = est.predict(datum).reshape(1, -1)[0]
    return preds

def intervals(data, percentile=95):
    """
    Given a numpy array of data, return the 0th, lower bound, median,
    upper bound and 100th percentile of the data. Generally useful for drawing
    box-plots.
    """
    low = (100 - percentile) / 2
    upp = 100 - low
    med = 50

    return np.percentile(data,
                         [0, low, med, upp, 100])

def predictions(drugs, models, seq):
    preds = list()  # we will store the data records-style
    # preds['drug'] = list()
    # preds['log10(DR)'] = list()
    # preds['yerr'] = list()
    for drug, mdl in zip(drugs, models):
        print(drug)
        pred = mdl.predict(seq)[0]
        prange = pred_range(mdl, seq)
        # zeroth, low, med, upp, hundreth = intervals(prange)
        # preds.append(dict(drug=drug, pred=pred))
        # preds['drug'].append(drug)
        # preds['log10(DR)'].append(pred)
        # preds['yerr'].append(prange.std())
        # print(preds)
        for p in prange:
            pred = dict()
            pred['drug'] = drug
            pred['log10(DR)'] = p
            preds.append(pred)
            print(pred)
    return preds
