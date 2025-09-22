import numpy as np
import pandas as pd
from os.path import exists

class table_generator(object):
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
    
    # dictionary for grouping the disorder variables
    disorders = {
        'ADD': ADD,
        'AGR': AGR,
        'CDD': CDD,
        'DEP': DEP,
        'GAD': GAD,
        'MAN': MAN,
        'OCD': OCD,
        'ODD': ODD,
        'PAN': PAN,
        'PHB': PHB,
        'PSY': PSY,
        'SIP': SIP,
        'PTD': PTD,
        'SEP': SEP,
        'SOC': SOC,
        'SUI': SUI
    }
    
    # list of possible values for the race variable
    races = ['AA', 'EA', 'OT', 'Mixed', 'Unspecified']
    
    # columns of the table to be generted
    cols = ['Block', 'n', 'M', 'F'] + races + list(disorders.keys())
    
    def __init__(self,
                 phenotype_data_filepath='[Insert filepath to phenotype data]',
                 SUBJID_filepath='../results/youth_SUBJIDs',
                 membership_list_filepath='../results/youth_multiplex_Pearson_memberships.csv'):
        
        # read the phenotype data and store it in a dataframe
        self.phenotype_df = pd.read_csv(phenotype_data_filepath, sep='\t', low_memory=False)
        
        # remove duplicate subjects from phenotype dataframe
        self.phenotype_df = self.phenotype_df[self.phenotype_df['INT_NUM'] == 1]
        
        # read the subject IDs and store them in a dataframe
        SUBJID_df = pd.read_csv(SUBJIDs_filepath)
        
        # remove subjects from the phenotype dataframe that aren't in the subject ID dataframe
        self.phenotype_df = pd.merge(self.phenotype_df, SUBJID_df['SUBJID'], how='inner', on='SUBJID')
        
        # read the list of block memberships for each subject
        membership_list = pd.read_csv(membership_list_filepath)['memberships'].values.tolist()
        
        # save the number of blocks
        self.Q = max(membership_list)
        
        # add the block memberships to the phenotype dataframe
        self.phenotype_df['Block'] = membership_list

        # create new dataframes for each block
        block_dfs= [self.phenotype_df[self.phenotype_df['Block'] == i+1] for i in range(self.Q)]
        
        # create a new empty table
        self.table = pd.DataFrame(columns=self.cols)
        
        for i in range(self.Q):
            # total number of subjects in a block
            n = block_dfs[i].shape[0]
            
            # create a new entry for the block to be added to the table
            row_dict = {'Block': i+1, 'n': n}
            
            # only proceed to generate results if there are subjects in the block
            if n > 0:
                block_df = block_dfs[i]
                
                # compute and store the proportion of each gender in the block
                for gender in ('M', 'F'):
                    row_dict[gender] = len(block_df[block_df['Sex'] == gender])/n
                
                # compute and store the proportion of each race in the block
                null_race_count = block_df['Race'].isnull().sum()
                for race in self.races:
                    if race == 'OT':
                        row_dict[race] = len(block_df[block_df['Race'].isin({'OT', 'PI', 'AI', 'HI', 'Eskimo/Alaskan'})])/n
                    elif race == 'Mixed':
                        row_dict[race] = (len(block_df[~block_df['Race'].isin({'AA', 'EA', 'OT', 'PI', 'AI', 'HI', 'Eskimo/Alaskan'})]) - null_race_count)/n
                    elif race == 'Unspecified':
                        row_dict[race] = null_race_count/n
                    else:
                        row_dict[race] = len(block_df[block_df['Race'] == race])/n
                
                disorder_list = []
                for disorder in self.disorders.keys():
                    disorder_list += self.disorders[disorder]
                
                block_df['SIP'] = 0
                for var in self.disorders['SIP']:
                    block_df['SIP'] = block_df['SIP'] + block_df[var]
                row_dict['SIP'] = np.mean(block_df['SIP'].to_numpy())
                    
                block_df = block_df.dropna(subset=disorder_list)
                
                row_dict['n_missing'] = row_dict['n'] - len(block_df)
                
                for disorder in self.disorders.keys():
                    if disorder != 'SIP':
                        block_df[disorder] = 0
                        for var in self.disorders[disorder]:
                            block_df[disorder] = block_df[disorder] + block_df[var]
                        row_dict[disorder] = np.mean(block_df[disorder].to_numpy())
                
            # add the new entry to the table
            self.table = pd.concat([self.table, pd.DataFrame([row_dict], index=[0])], ignore_index=True)


    # function for returning the table
    def table(self):
        return self.table