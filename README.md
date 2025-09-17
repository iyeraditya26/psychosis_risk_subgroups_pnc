# Salience Network Segregation and Symptom Profiles in Psychosis Risk Subgroups among Youth and Early Adults
This repository contains the code for the analyses presented in the paper, "Salience Network Segregation and Symptom Profiles in Psychosis Risk Subgroups among Youth and Early Adults". This code does not contain any files from the Philadelphia Neurodevelopment Cohort (PNC) dataset that was analyzed, as we do not own the dataset and are thus prevented from redistributing it. However, the PNC dataset can be obtained from [the database of Genotypes and Phenotypes (dbGaP)](http://www.ncbi.nlm.nih.gov/sites/entrez?db=gap) through dbGaP accession phs000607.v3.p2. The scripts can be run to reproduce our results once the PNC data is obtained from dbGaP, though this will require some very minor edits to the scripts in order to specify the appropriate filepaths.

## Directory Descriptions

### notebooks/
Jupyter and RStudio notebooks used to run the analyses
* : constructs the neuroimaging layer from functional connectivity matrices
* : constructs the symptom layer from self-reported psychosis risk symptom data
* : fits a multiplex Stochastic Block Model (SBM) to the neuroimaging and symptom layers
* : fits a simple SBM to only the symptom layer
* : Computes brain system segregation values for each subject
* : 

### src/
Source code for functions used by the notebooks
* : 
* : 

### results/
Output directory for the results produced by the notebooks
