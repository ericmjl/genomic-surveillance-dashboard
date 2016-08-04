import numpy as np
from sklearn.externals import joblib

def pred_range(model, datum):
    """
    Returns the full range of predictions for a given ensemble model.
    """
    preds = np.zeros(len(model.estimators_))
    for i, est in enumerate(model.estimators_):
        preds[i] = est.predict(datum).reshape(1, -1)[0]
    return preds

# def pred_intervals(data, percentile=95):
#     mean = np.mean(data)
#     lowbo = np.percentile(data, (100 - percentile) / 2. )
#     uppbo = np.percentile(data, 100 - (100 - percentile) / 2.)
#
#     return mean, lowbo, uppbo

def intervals(data, percentile=95):
    """
    Given a numpy array of data, return the 0th, lower bound, median,
    upper bound and 100th percentile of the data. Generally useful for drawing
    box-plots.
    """
    low = (100 - percentile) / 2
    upp = 100 - low
    med = 50

    return np.percentile(0, low, med, upp, 100)

def predict_withrange(drugs, seq):
    all_predictions = list()
    for drug in drugs:
        print(drug)
        model = joblib.load("../models/base/{drug}/{drug}.pkl".format(drug=drug))
        data = dict()
        data['drug'] = drug
        data['log10(DR)'] = pred_range(model, seq)
        print(data)
    all_predictions.append(data)
    return all_predictions

def predict_norange(drugs, seq):
    preds = dict()
    preds['drug'] = list()
    preds['log10(DR)'] = list()

    for drug in drugs:
        # print(drug)
        mdl = joblib.load("../models/base/{drug}/{drug}.pkl".format(drug=drug))
        pred = mdl.predict(seq)[0]
        # preds.append(dict(drug=drug, pred=pred))
        preds['drug'].append(drug)
        preds['log10(DR)'].append(pred)
        # print(pred)
    return preds
