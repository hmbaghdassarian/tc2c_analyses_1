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
    X : pd.DataFrame
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

from typing import Dict, List
import pandas as pd
from scipy.stats import hypergeom 
import statsmodels as sm
def loading_ora(lr_loadings: pd.Series, 
                lr_functions: Dict[str, List[str]],
                percentile: float = 0.9) -> pd.DataFrame:
    """Use an LR factor loading vector for over-representation analysis

    Parameters
    ----------
    lr_loadings : pd.Series
        a factors LR loadings; index is the LR pair, value is the loading 
    lr_functions : Dict[str, List[str]]
        value is the meta-term (e.g., signalling pathway)
        key is the list of LR pairs associated with the meta-term
    percentile: float, optional
        the top q-percentile loadings to consider, by default 0.9
    background : str, optional [DEPRECATED]
        the background genes to be considered, by default 'loadings' 
        Either 'loadings' for all LR pairs in the lr_loadings vector or 
        'universe' for all pairs in lr_functions

    Returns
    -------
    ora_res : pd.DataFrame
        summarizes over-representation analysis results for each term
    """
    background = 'universe'
    background = sorted(set([item for sublist in lr_functions.values() for item in sublist]))
#     if background == 'loadings':
#         background = lr_loadings.index.tolist()
#     elif background == 'universe':
#         background = sorted(set([item for sublist in lr_functions.values() for item in sublist]))
#     else: 
#         raise ValueError("Background must be on of ['loadings', 'universe']")
    
    lr_loadings = lr_loadings.sort_values(ascending = False)
    top_n = lr_loadings[lr_loadings >= lr_loadings.quantile(q = percentile)]
    
    
    ora_res = pd.DataFrame(columns = ['Term', 'Gene Ratio', 'BG Ratio', 'p_val'], 
                          index = range(len(lr_functions)))
    N = len(background)
    n = top_n.shape[0]
    for idx, meta_term in enumerate(lr_functions):
        M = len(lr_functions[meta_term]) # K 
        k = len(set(lr_functions[meta_term]).intersection(top_n.index))
        
        p_enrich = hypergeom.pmf(k,N,M,n) # TODO: make sure this is correct
        ora_res.loc[idx, :] = [meta_term, k/n, M/N, p_enrich]
        
    ora_res['BH_FDR'] = sm.stats.multitest.multipletests(pvals=ora_res.p_val, method='fdr_bh')[1]
    return ora_res
    