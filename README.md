# Salience Network Segregation and Symptom Profiles in Psychosis Risk Subgroups among Youth and Early Adults
This repository contains the code for the analyses presented in the paper, "Salience Network Segregation and Symptom Profiles in Psychosis Risk Subgroups among Youth and Early Adults". This code does not contain any files from the Philadelphia Neurodevelopment Cohort (PNC) dataset that was analyzed, as we do not own the dataset and are thus prevented from redistributing it. However, the PNC dataset can be obtained from [the database of Genotypes and Phenotypes (dbGaP)](http://www.ncbi.nlm.nih.gov/sites/entrez?db=gap) through dbGaP accession phs000607.v3.p2. The scripts can be run to reproduce our results once the PNC data is obtained from dbGaP, though this will require some very minor edits to the scripts in order to specify the appropriate filepaths.

## Directory Descriptions

### notebooks/
Jupyter and RStudio notebooks used to run the analyses and generate visualizations presented in the figures. Run the notebooks in the listed order to reproduce the results.
* create_neuroimaging_layers.ipynb: constructs the neuroimaging layer from functional connectivity matrices
* create_symptom_layers.ipynb: constructs the symptom layer from self-reported psychosis risk symptom data
* SBM_fitting.Rmd: fits all Stochastic Block Models (both simple and multiplex) that were used in the analyses.
* symptom_analysis.ipynb: produces tables containing aggregated psychopathology symptom scores for each subject. Also performs permutation tests to assess the significance of between-block differences in SIP score means.
* segregation_analysis.ipynb: Computes, describes, and visualizes brain system segregation values for each subject, breaking down these values by age group and the identified subject communities.
* general_analysis_results.ipynb: performs a wide range of analyses related to demographic information, psychopathology symptoms, and model fitting.
* sensitivity_analysis_results.ipynb: visualizes the parameter estimates and block-wise pychosis risk scores of the simple SBMs fit to only the symptom layer. It also visualizes the model fitting process for the simple SBM and the multiplex SBM fit to a nueroimaging layer constructed using pairwise Euclidean distances.
* community_detection_concordance_analysis.ipynb: performs analyses that evaluate the concordance between the community detection results from the simple and multiplex Stochastic Block Models.

### src/
Source code for functions used by the notebooks or for miscellaneous tasks
* conn_batch_execute.m: MATLAB script that preprocesses each subject's fMRI data using the CONN pipeline. Note: this script was submitted multiple times (once per subject) as a SLURM job on a high performance compute cluster.
* distance_FC.py: contains functions for computing the Pearson dissimilarities and Euclidean distances.
* neuroimaging_layer_constructor.py: contains functions that help to construct the neuroimaging layers.
* symptom_layer_constructor.py: contains functions that help to construct the symptom layers.
* table_generator.py: contains functions that help to produce tables containing aggregated psychopathology symptom scores for each subject.

### results/
Output directory for the results produced by the notebooks

### resources/
Miscellaneous files that are referenced or implemented in the analysis notebooks.
