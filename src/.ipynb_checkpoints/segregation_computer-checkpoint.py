import numpy as np
import pandas as pd
from os.path import exists
from scipy.io import loadmat

class segregation_computer(object):
    
    # columns of the table to be generted
    cols = ['Block', 'Total FC Within Network', 'Mean FC Within Network', 'Segregation Between All Other Regions', 'Segregation Between Default Network', 'Segregation Between Salience Network', 'Segregation Between FrontoParietal Network', 'Segregation Between DorsalAttn Network']
    
    # networks that were analyzed
    networks = ['Default', 'Salience', 'FrontoParietal', 'DorsalAttn']
    
    def __init__(self,
                 parcels_filepath='../resources/Parcels.xlsx',
                 fc_matrix_dir = '[Insert filepath to the directory containing the functional connectivity matrices]',
                 fc_matrix_suffix = '_resultsROI_Condition001.mat',
                 SUBJID_filepath='../results/youth_SUBJIDs',
                 membership_list_filepath='../results/youth_multiplex_Pearson_memberships.csv',
                 network='Salience'):
        
        # read the parcels and store them in a dataframe
        self.parcels = pd.read_excel(parcels_filepath)
        
        assert network in set(self.parcels['Community'].tolist()) | set(['']), 'Invalid network. Valid networks: Auditory, CinguloOperc, CinguloParietal, Default, DorsalAttn, FrontoParietal, None, RetrosplenialTemporal, SMhand, SMmouth, Salience, VentralAttn, Visual. Note: you can enter an empty string if you want to include all regions of interest in computation.'
        
        # read the subject IDs and store them in a dataframe
        SUBJID_df = pd.read_csv(SUBJID_filepath)
        
        # read the list of block memberships for each subject
        membership_list = pd.read_csv(membership_list_filepath)['memberships'].values.tolist()
        
        # save the number of blocks
        self.Q = max(membership_list)
        
        # Add the block memberships of each subject to the subject ID dataframe
        SUBJID_df['Block'] = membership_list
        
        # create new dataframes for each block
        block_dfs = [SUBJID_df[SUBJID_df['Block'] == i+1] for i in range(self.Q)]
        
        # create a new empty table
        cols = self.cols.copy()
        cols.remove(f'Segregation Between {network} Network')
        self.table = pd.DataFrame(columns=cols)
        
        # remove the network of interest from the list of networks we want to study
        networks = self.networks.copy()
        networks.remove(network)
        
        for i in range(self.Q):
            # total number of subjects in a block
            n = block_dfs[i].shape[0]
            
            # create a new entry for the block to be added to the table
            row_dict = {'Block': i+1}
            
            # only proceed to generate results if there are subjects in the block
            if n > 0:
                block_df = block_dfs[i]
                
                # obtain all the networks other than the network of interest
                all_networks = self.parcels['Community'].tolist()
                all_networks.remove(network)
                
                # create a list
                segregation_lists = {net:[] for net in networks}
                segregation_lists['All'] = []
                segregation_lists['Within Sum'] = []
                segregation_lists['Within Mean'] = []
                
                for j,subject in enumerate(block_df['SUBJID'].tolist()):
                    
                    # find filepath to functional connectivity matrix
                    if exists(fc_matrix_dir + subject + fc_matrix_suffix):
                        output_path = fc_matrix_dir + subject + fc_matrix_suffix
                    else:
                        raise Exception(f'Functional connectivity matrix for subject {subject} does not exist.')
                    
                    # load functional connectivity matrix
                    mat = loadmat(output_path)['Z'][:333,:333]
                    
                    # add segregation values of subject to the segregation lists
                    segregation_lists['Within Mean'].append(self.FC_mean(mat, [network], [network]))
                    segregation_lists['Within Sum'].append(self.FC_sum(mat, [network], [network]))
                    segregation_lists['All'].append(self.segregation(mat, [network], all_networks))
                    for net in networks:
                        segregation_lists[net].append(self.segregation(mat, [network], [net]))
                
                # find the means of all the segregation values
                row_dict['Mean FC Within Network'] = np.mean(np.array(segregation_lists['Within Mean']))
                row_dict['Total FC Within Network'] = np.mean(np.array(segregation_lists['Within Sum']))
                row_dict['Segregation Between All Other Regions'] = np.mean(np.array(segregation_lists['All']))
                for net in networks:
                    row_dict[f'Segregation Between {net} Network'] = np.mean(np.array(segregation_lists[net]))
            
            # add the new entry to the table
            self.table = self.table.append(row_dict, ignore_index=True)
            
    # function for returning the table
    def table(self):
        return self.table

    # function for computing the sum of functional connectivities between or within networks
    def FC_sum(self, mat, nets1, nets2):

        mat = mat.copy()
        mat = np.nan_to_num(mat, nan=0)
        mat[mat<0] = 0

        nets1_rois = self.parcels.loc[self.parcels['Community'].isin(nets1)].index[:]
        nets2_rois = self.parcels.loc[self.parcels['Community'].isin(nets2)].index[:]

        # To get FC in block diagonal or off diagonal
        conn_mat_x = mat.take(nets1_rois, axis=0)
        conn_mat_xy = conn_mat_x.take(nets2_rois, axis=1)
        
        return conn_mat_xy.sum()
    
    # function for computing the mean of functional connectivities between or within networks
    def FC_mean(self, mat, nets1, nets2):

        mat = mat.copy()
        mat = np.nan_to_num(mat, nan=0)
        mat[mat<0] = 0

        nets1_rois = self.parcels.loc[self.parcels['Community'].isin(nets1)].index[:]
        nets2_rois = self.parcels.loc[self.parcels['Community'].isin(nets2)].index[:]

        # To get FC in block diagonal or off diagonal
        conn_mat_x = mat.take(nets1_rois, axis=0)
        conn_mat_xy = conn_mat_x.take(nets2_rois, axis=1)
        
        return conn_mat_xy.mean()

    # function for computing the segregation between two different sets of regions
    def segregation(self, mat, nets1, nets2):

        mat = mat.copy()
        mat = np.nan_to_num(mat, nan=0)
        mat[mat<0] = 0

        nets1_rois = self.parcels.loc[self.parcels['Community'].isin(nets1)].index[:]
        nets2_rois = self.parcels.loc[self.parcels['Community'].isin(nets2)].index[:]
        
        target_net_x = mat.take(nets1_rois, axis=0)
        target_net_xy = target_net_x.take(nets1_rois, axis=1)

        outside_net_x = mat.take(nets2_rois, axis=0)
        outside_net_xy = outside_net_x.take(nets2_rois, axis=1)

        # (Mean of target network - mean of outside networks) / mean target network
        # Only works if nets1 != nets2
        segregation = (np.mean(np.abs(target_net_xy)) - np.mean(np.abs(outside_net_xy))) / np.mean(np.abs(target_net_xy))

        return segregation