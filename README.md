# Ensemble-tree-for-risk-prediction-in-MM

This file describes the code layout of the ensemble tree in MM project.

## Sample code

The `real_data_experiment.Rmd` file shows an example of training the ensemble tree-like model to predict survival outcome. 

- The ranger() function from the Ranger package is used to the initial train decision tree. 

- The pairwise_survdiff() function from the Survminer package is used to perform pairwise log-rank tests to find the nearest leaf nodes.

- The survfit() function from the Survival package is used to fit Nelson-Aalen esimators for each combined leaf node.

- The cIndex() function from the Intsurv package is used to compute the concordence index of our tree-based models.

- The constrOptim.nl() function from the Alabama package is used to find the optimal ensemble weights for the ensemble tree-like model.

## Implementation in MMRF/IFM/Gamer/UAMS

The folders include Rmd files implementing risk prediction models (EMC-92, UAMS-70, EI-score, Gene Proliferation Index, UAMS-17, UAMS-80, IFM-15, HM19) using gene expression levels to predict risk of MM patients using each dataset. Gene expression level data and patient survival data of each MM study was used. Brainarray platform and biomaRt package were used to convert probe set IDs to gene ensembl IDs.

The brainarray file can be downloaded at http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/25.0.0/ensg.asp.

## Ensemble tree with risk prediction scores

The folder contains csv files for respective risk scores obtained from risk prediction models.

- The `risk score data experiment.Rmd` file contains example code for applying ensemble method on the data frame containing risk scores and survival outcome.

- The `risk score data experiment bootstrap.Rmd` file runs 500 bootstapping simulations to compare model performance on different datasets (continuous score and binary class).

- The `risk score data experiment with pruning.Rmd` file runs 500 bootstapping simulations to explore performance of ensemble model on pruned trees (using rpart package).

## Optimal initial tree

The folder contains exploration on methods to develop an optimal tree to begin the ensemble with using the MMRF data.

- The files including STreeD in their names explore the performance of ensemble model on optimal survival trees (using STreeD package in python). STreeD is developed in python (https://github.com/AlgTUDelft/pystreed), the `.py` files call the function in python.

- The files including rpart in their names explore the performance of ensemble model with initial tree fit by rpart package in R.

- The files including OST in their names explore the performance of ensemble model on optimal survival trees (using OST algorithm implemented in Julia). OST can be called using the IAI package in R. https://docs.interpretable.ai/stable/IAI-R/quickstart/ot_survival/

## Cross Study Validation

The folder contains exploration on ensemble tree performance while validated on an independent dataset. We train the model using the MMRF dataset and validate the models on the IFM dataset, after implementing all risk prediction models in the IFM dataset. `Cross Study experiment in MMRF and IFM` includes codes to train and validate the methods in MMRF and IFM. 

Later we added the UAMS dataset and GAMER dataset. `Cross Study Validation` included codes to study CSV behavior.
