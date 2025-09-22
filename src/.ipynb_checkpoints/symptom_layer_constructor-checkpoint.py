import numpy as np
import pandas as pd
import itertools as it
from numpy import linalg as LA
from os.path import exists
from scipy.spatial import distance

class symptom_layer_constructor(object):
    # define lists for the different categories of variables
    ID_vars = ['SUBJID', 'INT_NUM', 'INT_TYPE', 'Sex']
    ADD = ['ADD011', 'ADD012', 'ADD013', 'ADD014', 'ADD015', 'ADD016', 'ADD020', 'ADD021', 'ADD022']
    AGR = ['AGR001', 'AGR002', 'AGR003', 'AGR004', 'AGR005', 'AGR006', 'AGR007', 'AGR008']
    CDD = ['CDD001', 'CDD002', 'CDD003', 'CDD004', 'CDD005', 'CDD006', 'CDD007', 'CDD008']
    DEP = ['DEP001', 'DEP002', 'DEP004', 'DEP006']
    GAD = ['GAD001', 'GAD002']
    MAN = ['MAN001', 'MAN002', 'MAN003', 'MAN004', 'MAN005', 'MAN006', 'MAN007']
    OCD = ['OCD001', 'OCD002', 'OCD003', 'OCD004', 'OCD005', 'OCD006', 'OCD007', 'OCD008', 'OCD011', 'OCD012', 'OCD013', 'OCD014', 'OCD015', 'OCD016', 'OCD017', 'OCD018', 'OCD019']
    ODD = ['ODD001', 'ODD002', 'ODD003', 'ODD005', 'ODD006']
    PAN = ['PAN001', 'PAN003', 'PAN004']
    PHB = ['PHB001', 'PHB002', 'PHB003', 'PHB004', 'PHB005', 'PHB006', 'PHB007', 'PHB008']
    PSY = ['PSY001', 'PSY029', 'PSY050', 'PSY060', 'PSY070', 'PSY071']
    SIP = ['SIP003', 'SIP004', 'SIP005', 'SIP006', 'SIP007', 'SIP008', 'SIP009', 'SIP010', 'SIP011', 'SIP012', 'SIP013', 'SIP014']
    PTD = ['PTD001', 'PTD002', 'PTD003', 'PTD004', 'PTD006', 'PTD007', 'PTD008', 'PTD009']
    SEP = ['SEP500', 'SEP508', 'SEP509', 'SEP510', 'SEP511']
    SOC = ['SOC001', 'SOC002', 'SOC003', 'SOC004', 'SOC005']
    SUI = ['SUI001', 'SUI002']
    
    # the columns of the phenotypic data are the variables in the list above
    cols = ID_vars + ADD + AGR + CDD + DEP + GAD + MAN + OCD + ODD + PAN + PHB + PSY + SIP + PTD + SEP + SOC + SUI
    
    def __init__(self,
                 phenotype_data_filepath='[Insert filepath to phenotype data]',
                 quality_control_data_filepath='[Insert filepath to the quality control data]',
                 fc_matrix_dir = '[Insert filepath to the directory containing the functional connectivity matrices]',
                 fc_matrix_suffix = '_resultsROI_Condition001.mat'):
        
        # read the phenotype data and store it in a dataframe
        self.phenotype_df = pd.read_csv(phenotype_data_filepath, sep='\t', usecols=symptom_layer_constructor.cols, low_memory=False)
        
        # remove duplicate subjects from phenotype dataframe
        self.phenotype_df = self.phenotype_df[self.phenotype_df['INT_NUM'] == 1]
        
        # read the quality control data and store it in a dataframe
        self.QC_df = pd.read_csv(quality_control_data_filepath)
        
        # remove samples with a valid scan percentage less than or equal to 80 from quality control dataframe
        self.QC_df = self.QC_df[self.QC_df['QC_ValidScanPercentage'] > 80]
        
        # remove all subjects from phenotype dataframe that are not in quality control dataframe
        self.phenotype_df = pd.merge(self.phenotype_df, self.QC_df['SUBJID'], how='inner', on='SUBJID')
        
        # remove all rows with missing data from phenotype dataframe
        for column in symptom_layer_constructor.SIP:
            self.phenotype_df = self.phenotype_df[self.phenotype_df[column].isnull() == False]

        # remove all rows with subjects that do not have a functional connectivity matrix from phenotype dataframe
        for subject in self.phenotype_df['SUBJID']:
            if not exists(fc_matrix_dir + subject + fc_matrix_suffix):
                self.phenotype_df = self.phenotype_df.drop(self.phenotype_df[self.phenotype_df['SUBJID'] == subject].index)
        
        # assign a value of 0.0 to all cells with unknown values in phenotype dataframe
        self.phenotype_df[symptom_layer_constructor.cols[3:]] = self.phenotype_df[symptom_layer_constructor.cols[3:]].applymap(lambda x: 0.0 if x == 9.0 else x)

    # function for returning the fully processed dataframe
    def dataFrame(self):
        return self.phenotype_df

    # function for computing the value of a single element in the matrix
    @classmethod
    def __compute_correlation(self, subj_1, subj_2, VI=None):
        assert VI is not None, 'Inverse of the covariance matrix (VI) is set to None.'
        vec1 = subj_1[symptom_layer_constructor.SIP].to_numpy(dtype=np.float32)
        vec2 = subj_2[symptom_layer_constructor.SIP].to_numpy(dtype=np.float32)
        return distance.mahalanobis(vec1, vec2, VI)

        # adjust value of similarity score so that Fisher transform function does not return inf or -inf
        if similarity == 1:
            similarity = 0.99
        elif similarity == -1:
            similarity = -0.99

        # return Fisher transform of similarity score
        return np.arctanh(similarity)
    
    # function for generating the matrix
    def generate_matrix(self, age_group='youth', SUBJIDs_filepath='../results/youth_SUBJIDs', output_filepath=''):
        
        if age_group == 'youth':
            age_group_df = self.phenotype_df[self.phenotype_df['INT_TYPE'] == 'MP']
            
            # store SUBJID column of DataFrame as a csv file
            age_group_df['SUBJID'].to_csv(SUBJIDs_filepath, index=False)
        
        elif age_group == 'early adults':
            age_group_df = self.phenotype_df[self.phenotype_df['INT_TYPE'] == 'AP']
            
            # store SUBJID column of DataFrame as a csv file
            age_group_df['SUBJID'].to_csv(SUBJIDs_filepath, index=False)
        
        else:
            raise Exception('Invalid value for age_group. Valid values for age_group: "youth" or "early adults".')
        
        n_subjects = len(age_group_df)
        
        # initialize a matrix with a default value of 0.0 for each element
        matrix = np.zeros(n_subjects, n_subjects))
    
        VI = LA.pinv(np.cov(age_group_df[symptom_layer_constructor.SIP].to_numpy().T))
            
        # iterate through each row and column
        for i,j in it.product(range(n_subjects), repeat=2):
            if (i%20 == 0 and j==0):
                print(f'Finished creating {i} out of {n_subjects} rows')
            matrix[i][j] = symptom_layer_constructor.__compute_correlation(age_group_df.iloc[i], age_group_df.iloc[j], VI=VI)
        
        if output_filepath != '':
            matrix_df = pd.DataFrame(matrix)
            matrix_df.to_csv(output_filepath, index=False)
            print(f'Matrix has been saved as a csv file with the following path: {output_filepath}')
        
        return matrix