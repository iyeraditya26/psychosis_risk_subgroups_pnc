# Salience Network Segregation and Symptom Profiles in Psychosis Risk Subgroups among Youth and Early Adults
This repository contains the code for the analyses presented in the paper, "Salience Network Segregation and Symptom Profiles in Psychosis Risk Subgroups among Youth and Early Adults". This code does not contain any files from the Philadelphia Neurodevelopment Cohort (PNC) dataset that was analyzed, as we do not own the dataset and are thus prevented from redistributing it. However, the PNC dataset can be obtained from [the database of Genotypes and Phenotypes (dbGaP)](http://www.ncbi.nlm.nih.gov/sites/entrez?db=gap) through dbGaP accession phs000607.v3.p2. The scripts can be run to reproduce our results once the PNC data is obtained from dbGaP, though this will require some very minor edits to the scripts in order to specify the appropriate filepaths.

## Directory Descriptions

### notebooks/
Jupyter and RStudio notebooks used to run the analyses. Run the notebooks in the listed order to reproduce the results.
* create_neuroimaging_layers.ipynb: constructs the neuroimaging layer from functional connectivity matrices
* create_symptom_layers.ipynb: constructs the symptom layer from self-reported psychosis risk symptom data
* SBM_fitting.Rmd: fits all Stochastic Block Models (both simple and multiplex) that were used in the analyses.
* symptom_analysis.ipynb: produces tables .
* segregation_analysis.ipynb: Computes, describes, and visualizes brain system segregation values for each subject, breaking down these values by age group and the identified subject communities.
* general_analysis_results.ipynb: .
* sensitivity_analysis_results.ipynb: .
* community_detection_concordance_analysis.ipynb: .

### src/
Source code for functions used by the notebooks
* : 
* : 

### results/
Output directory for the results produced by the notebooks

### resources/

