#!/usr/bin/env python
# coding: utf-8

# Functions to add to tensor-cell2cell

from typing import Union
import numpy as np
import pandas as pd
def get_joint_loadings(self, dim1: str, dim2: str, factor: Union[str,int]) -> pd.DataFrame:
    """Create the joint probability distribution between two tensor dimensions for a 
    given factor output from decomposition.
    
    Parameters
    ----------
    dim1 : str
        One of the tensor dimensions (options are self.factors.keys())
    dim2 : str
        A second tensor dimension (options are self.factors.keys())
    factor: Union[str,int]
        one of the factors output from decompositon [1, self.rank]
    
    Returns
    -------
    joint_dist : pd.DataFrame
        joint distribution of factor loadings for the specified dimensions. 
        Rows correspond to dim1 and columns to dim2
    """
    if dim1 not in self.factors:
        raise ValueError('The specified dimension ' + dim1 + ' is not present in the tensor')
    if dim2 not in self.factors:
        raise ValueError('The specified dimension ' + dim2 + ' is not present in the tensor')
    if int(factor) > list(self.factors.values())[0].shape[1]:
        raise ValueError('The specified factor lies outside the computed low-rank tensor')
    
    vec1 = self.factors[dim1]['Factor ' + str(factor)]
    vec2 = self.factors[dim2]['Factor ' + str(factor)]
    
    # calculate the outer product
    joint_dist = pd.DataFrame(np.outer(vec1, vec2), 
                              index = vec1.index, 
                              columns = vec2.index)

    return joint_dist


from typing import Dict
import pandas as pd
import scipy.cluster.hierarchy as hcluster
from scipy.stats import zscore
def get_hierarchical_clusters(X: pd.DataFrame, height: float, 
                              z_rows = False, z_cols = False, 
                             **kwargs) -> Dict[str, int]:
    '''Assign LR pairs to a cluster using hierarchical clusters.
    
    Parameters
    ----------
    X : Union[np.ndarray, pd.DataFrame]
        Index is the LR pairs. Should be the factor loadings or a factor specific 
        joint distribution
    height : float
        t parameter in scipy.cluster.hierarchy.fcluster. Allows adjusting of # of clusters output. 
    z_rows : bool, optional
        z-score the rows, by default False
    z_cols : bool, optional
        z-score the columns, by default False
    **kwargs : 
        Parameters for scipy.cluster.hierarchy.linkage
    
    Returns
    -------
    clusters : Dict[str, int]
        Keys are the LR pair, values are the assigned cluster label
    '''
    if z_rows:
        X = X.T.apply(zscore).T
    if z_cols:
        X = X.apply(zscore)
        
    Z = hcluster.linkage(X, **kwargs)
    cluster_labels = hcluster.fcluster(Z, t=height, criterion='distance')
    clusters = dict(zip(X.index, cluster_labels))
    
    return clusters