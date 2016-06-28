"""
A Python module that automatically selects best parameters for a specified
regression ensemble learning model.
"""

from sklearn.ensemble import (AdaBoostRegressor,
                              BaggingRegressor,
                              ExtraTreesRegressor,
                              GradientBoostingRegressor,
                              RandomForestRegressor)
from sklearn.grid_search import GridSearchCV
import numpy as np

shortnames = dict()
shortnames['abr'] = 'AdaBoost Regressor'
shortnames['bgr'] = 'Bagging Regressor'
shortnames['etr'] = 'ExtraTrees Regressor'
shortnames['gbr'] = 'GradientBoosting Regressor'
shortnames['rfr'] = 'RandomForest Regressor'

models = dict()
models['abr'] = AdaBoostRegressor()
models['bgr'] = BaggingRegressor()
models['etr'] = ExtraTreesRegressor()
models['gbr'] = GradientBoostingRegressor()
models['rfr'] = RandomForestRegressor()

params = dict()
params['abr'] = {'n_estimators': np.arange(50, 501, 50),
                 'learning_rate': np.arange(0.1, 1.1, 0.1),
                 }
params['bgr'] = {'n_estimators': np.arange(10, 100, 20),
                 'bootstrap': [False, True],
                 'bootstrap_features': [False, True],
                 }
params['etr'] = {'max_features': np.arange(0.1, 1.0, 0.2),
                 'bootstrap': [False, True],
                 }
params['gbr'] = {'n_estimators': np.arange(100, 701, 100),
                 # 'max_depth': np.arange(3, 6, 1),
                 'learning_rate': np.arange(0.02, 0.21, 0.04),
                 # 'subsample': np.arange(0.6, 1.01, 0.1)
                 }
params['rfr'] = {'max_features': np.arange(0.1, 1.1, 0.2),
                 'bootstrap': [False, True],
                 }


def find_best_params(mdl, cv, scoring, X, Y):
    """
    Uses scikit-learn's GridSearchCV class to search across reasonable
    parameter range defaults, which are in turn specified above.
    """
    assert mdl in models.keys(), "mdl must be one of {0}".format(models.keys())

    gs = GridSearchCV(models[mdl], params[mdl], n_jobs=-1, verbose=3, cv=cv,
                      scoring=scoring)
    gs.fit(X, Y)

    return gs
