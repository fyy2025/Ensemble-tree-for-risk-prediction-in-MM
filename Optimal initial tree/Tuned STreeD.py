from pystreed import STreeDSurvivalAnalysis
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.metrics import accuracy_score
import numpy as np
import pandas as pd
import time
from sksurv.metrics import concordance_index_censored
from sklearn.base import BaseEstimator, ClassifierMixin

def fit(data):
    def concordance_index_scorer(estimator, X, y):
        # Predict risk scores
        risk_scores = estimator.predict(X)
        
        # Extract event indicators and event times
        event_indicators = y['status'].astype(bool)
        event_times = y['time']
        
        # Compute concordance index
        c_index = concordance_index_censored(event_indicators, event_times, risk_scores)[0]
        
        return c_index


    df = data
    X = df.iloc[:,0:4].values
    y = df[['status',"time"]].to_records(index=False, column_dtypes={"time": np.double, 'status': bool})

    ##################################################################
    ##### 1. Tune using GridSearchCV from sklearn ####################
    ##################################################################

    model = STreeDSurvivalAnalysis(max_depth=4,min_leaf_node_size=data.shape[0]//10)

    params = [{ 'max_num_nodes': list(range(1,16)), 'n_thresholds': list(range(1,10))}]

    gs_knn = GridSearchCV(model,
                        param_grid=params,
                        scoring=concordance_index_scorer,
                        cv=5,
                        n_jobs=1,
                        verbose=0)
    start = time.perf_counter()
    gs_knn.fit(X, y)
    gs_duration = time.perf_counter() - start
    best_params = gs_knn.best_params_
    best_model = STreeDSurvivalAnalysis(max_depth=4,**best_params)

    best_model.fit(X, y)
    return best_model

