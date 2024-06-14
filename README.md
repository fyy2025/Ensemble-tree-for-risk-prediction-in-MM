# Ensemble-tree-for-risk-prediction-in-MM

This file describes the code layout of the ensemble tree in MM project.

## Sample code

The `real_data_experiment.Rmd` file shows an example of training the ensemble tree-like model to predict survival outcome. 

- The ranger() function from the Ranger package is used to the initial train decision tree. 

- The pairwise_survdiff() function from the Survminer package is used to perform pairwise log-rank tests to find the nearest leaf nodes.

- The survfit() function from the Survival package is used to fit Nelson-Aalen esimators for each combined leaf node.

- The cIndex() function from the Intsurv package is used to compute the concordence index of our tree-based models.

- The constrOptim.nl() function from the Alabama package is used to find the optimal ensemble weights for the ensemble tree-like model.

## Implementation of risk prediction models

The folder includes Rmd files implementing risk prediction models (EMC-92, UAMS-70) using gene expression levels to predict risk of MM patients using the MMRF data. Gene expression level data and patient survival data of the MMRF study was used. Brainarray platform and biomaRt package were used to convert probe set IDs to gene ensembl IDs.

The brainarray file can be downloaded at http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/25.0.0/ensg.asp.

## Ensemble tree with risk prediction scores

The folder contains csv files for respective risk scores obtained from risk prediction models.

- The `risk score data experiment.Rmd` file contains example code for applying ensemble method on the data frame containing risk scores and survival outcome.

- The `risk score data experiment bootstrap.Rmd` file runs 500 bootstapping simulations to compare model performance on different datasets (continuous score and binary class).

- The `risk score data experiment with pruning.Rmd` file runs 500 bootstapping simulations to explore performance of ensemble model on pruned trees (using rpart package).

## Optimal initial tree

The folder contains exploration on methods to develop an optimal tree to begin the ensemble with.

- The `SurTree.Rmd` file runs 500 bootstapping simulations to explore performance of ensemble model on optimal survival trees (using STreeD package in python). STreeD is developed in python (https://github.com/AlgTUDelft/pystreed), the `SurTree Implementation.py` calls the function in python.
