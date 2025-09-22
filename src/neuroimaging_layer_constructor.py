import numpy as np
import pandas as pd
import itertools
from os.path import exists
from scipy.io import loadmat
from distance_FC import distance_FC

class neuroimaging_layer_constructor(object):
    
    def __init__(self,
                 SUBJID_filepath='../results/youth_SUBJIDs',
                 fc_matrix_dir = '[Insert filepath to the directory containing the functional connectivity matrices]',
                 fc_matrix_suffix = '_resultsROI_Condition001.mat'):
        self.subjects = [str(subject) for subject in pd.read_csv(SUBJID_filepath)['SUBJID'].to_list()]
    
        self.fc_matrices = np.zeros((len(self.subjects), 333, 333))
        
        for i,subject in enumerate(self.subjects):
            if exists(fc_matrix_dir + subject + fc_matrix_suffix):
                output_path = fc_matrix_dir + subject + fc_matrix_suffix
            else:
                raise Exception(f'Functional connectivity matrix for subject {subject} does not exist.')
            self.fc_matrices[i] = np.nan_to_num(loadmat(output_path)['Z'][:333,:333])

    # function for retrieving the list of subjects
    def subjects(self):
        return self.subjects

    # function for retrieving the functional connectivity matrices
    def fc_matrices(self):
        return self.fc_matrices

    # function for computing the value of a single element in the matrix
    @classmethod
    def __distance(self, FC1, FC2, dist_type='Pearson'):
        dist = distance_FC(FC1, FC2)
        if (dist_type == 'Euclidean'):
            return dist.euclidean()
        elif (dist_type == 'Pearson'):
            return dist.pearson()
        else:
            raise Exception(f'Invalid argument for dist_type. Expected \'Euclidean\' or \'Pearson\' but received \'{dist_type}\' instead.')

    # function for generating the matrix
    def generate_matrix(self,
                        parcels_filepath='../resources/Parcels.xlsx',
                        network='',
                        dist_type='Pearson',
                        output_filepath=''):
        parcels = pd.read_excel(parcels_filepath)
        
        assert network in set(parcels['Community'].tolist()) | set(['']), 'Invalid network. Valid networks: Auditory, CinguloOperc, CinguloParietal, Default, DorsalAttn, FrontoParietal, None, RetrosplenialTemporal, SMhand, SMmouth, Salience, VentralAttn, Visual. Note: you can enter an empty string if you want to include all regions of interest in computation.'
        
        del_indices = [] if network=='' else parcels.index[parcels['Community'] != network].tolist()

        num_subjects = len(self.subjects)
        matrix = np.zeros((num_subjects, num_subjects))

        for i,j in itertools.product(range(num_subjects), repeat=2):
            if (i%20==0 and j==0):
                print(f'Finished computing {i} out of {num_subjects} rows')

            FC1, FC2 = self.fc_matrices[i].copy(), self.fc_matrices[j].copy()
            FC1 = np.delete(FC1, del_indices, axis=0)
            FC1 = np.delete(FC1, del_indices, axis=1)
            FC2 = np.delete(FC2, del_indices, axis=0)
            FC2 = np.delete(FC2, del_indices, axis=1)

            matrix[i,j] = neuroimaging_data_processor_updated.__distance(FC1, FC2, dist_type=dist_type)

        print(f'Finished computing all {num_subjects} rows')

        if output_filepath != '':
            matrix_df = pd.DataFrame(matrix)
            matrix_df.to_csv(output_filepath, index=False)
            print(f'Matrix has been saved as a csv file with the following path: {output_filepath}')

        return matrix