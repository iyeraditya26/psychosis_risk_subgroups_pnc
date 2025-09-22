import numpy as np
from numpy import linalg as LA
from scipy.spatial import distance

class distance_FC(object):
    def __init__(self, FC1, FC2):
        self.FC1 = FC1
        self.FC2 = FC2

        # ensure symmetric
        self.FC1 = self._ensure_symmetric(self.FC1)
        self.FC2 = self._ensure_symmetric(self.FC2)

    def _info(self, s):
        print('INFO: %s' % s)

    def _ensure_symmetric(self, Q):
        '''
        computation is sometimes not precise (round errors),
        so ensure matrices that are supposed to be
        symmetric are symmetric
        '''
        return (Q + np.transpose(Q))/2

    def _vectorize(self, Q):
        '''
        given a symmetric matrix (FC), return unique
        elements as an array. Ignore diagonal elements
        '''
        # extract lower triangular matrix
        tri = np.tril(Q, -1)

        vec = []
        for ii in range(1, tri.shape[0]):
            for jj in range(ii):
                vec.append(tri[ii, jj])
        
        return np.asarray(vec)

    def pearson(self):
        '''
        conventional Pearson distance between
        two FC matrices. The matrices are vectorized
        '''
        vec1 = self._vectorize(self.FC1)
        vec2 = self._vectorize(self.FC2)

        return (1 - np.corrcoef(vec1, vec2)[0, 1])/2
    
    def euclidean(self):
        '''
        conventional Euclidean distance between
        two FC matrices. The matrices are vectorized
        '''
        vec1 = self._vectorize(self.FC1)
        vec2 = self._vectorize(self.FC2)

        return LA.norm(vec1 - vec2)
