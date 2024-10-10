from pystreed import STreeDSurvivalAnalysis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.metrics import concordance_index_censored

def fit(data):   

    # df = pd.read_csv("/Users/fanyiyang/Desktop/Research/summer/Implementation of risk models/merge data continuous.csv")
    df=data
    X = df.iloc[:,0:4].values
    y = df[['status',"time"]].to_records(index=False, column_dtypes={"time": np.double, 'status': bool})
    events = y["status"]
    times = y["time"]

    # Train an optimal survival tree model
    model = STreeDSurvivalAnalysis(min_leaf_node_size=data.shape[0]//10)
    model.fit(X, y)
    return model

# Measure the performance of the model
# prediction = model.predict(X)
# result = concordance_index_censored(events, times, prediction)

# model.export_dot("tree.dot")

# import graphviz
# with open("tree.dot") as f:
#     dot_graph = f.read()
# g = graphviz.Source(dot_graph)
# g.render('tree.dot', outfile="tree.pdf", view=True, cleanup=True)