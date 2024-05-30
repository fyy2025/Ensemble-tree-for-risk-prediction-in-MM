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
