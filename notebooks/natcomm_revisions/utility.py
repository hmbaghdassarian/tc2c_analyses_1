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


# going from a communication edge list (LR, CC, score) or matrix (sender-receiver columns x ligand^receptor pairs) to a communication tensor
from typing import Dict, List, Optional
from collections import OrderedDict
from numpy.testing import assert_almost_equal
import warnings
import pandas as pd
from tqdm import tqdm
import tensorly as tl
from cell2cell.tensor import PreBuiltTensor

def _lr_to_matrix(df, lr_pair, cell_delim='-', lr_delim='^'):
    '''Convert a row of a communication matrix to a slice of the 3D tensor'''
    slice_df = pd.DataFrame(df.loc[lr_pair, :])
    slice_df = pd.concat([pd.Series(slice_df.index).str.split(cell_delim, expand = True),slice_df.reset_index(drop = True)], 
             axis = 1) 
    slice_df = slice_df.pivot(index = 0, columns=1).values
    
    return slice_df

def _matrix_to_3d_tensor(df, cell_delim='-', lr_delim='^'):
    """Reformat from CC-LR pair to 3D tensor
    
    df: pd.DataFrame
        columns are sender-receiver cell pairs, separated by '-'
        index is ligand-receptor pairs, separated by '^'
    """
    tensor_3d = np.stack([_lr_to_matrix(df=df, lr_pair=lr_pair, cell_delim=cell_delim, lr_delim=lr_delim) for lr_pair in df.index])
    return tensor_3d

# replace lines 756-776 in c2c.tensor.tensor (since used by multiple functions)
def get_cells_and_lrs(df_list: List, lr_how: str = 'outer', cell_how: str ='outer'):
    """Generate a comprehensive list of LR pairs and cells (or cell pairs) from a list of dataframes.
    Deals with inconsistencies in indeces across context-specific matrices. 

    Parameters
    ----------
    df_list : List
        a list of expression or communication matrices, each corresponding to a different 
    lr_how : str, optional
        take the union (outer) or intersection ('inner') of all matrix row names, by default 'outer'
    cell_how : str, optional
        take the union (outer) or intersection ('inner') of all matrix column names, by default 'outer'

    Returns
    -------
    [type]
        [description]
    """
    df_idxs = [list(df.index) for df in df_list]
    df_cols = [list(df.columns) for df in df_list]
    
    if lr_how == 'outer':
        lr_pairs = set.union(*map(set, df_idxs))
    elif lr_how == 'inner':
        lr_pairs = set.intersection(*map(set, df_idxs))
    
    if cell_how == 'outer':
        cells = set.union(*map(set, df_cols))
    elif cell_how == 'inner':
        cells = set.intersection(*map(set, df_cols))
        
    # Preserve order or sort new set (either inner or outer)
    if set(df_idxs[0]) == lr_pairs:
        lr_pairs = df_idxs[0]
    else:
        lr_pairs = sorted(lr_pairs)

    if set(df_cols[0]) == cells:
        cells = df_cols[0]
    else:
        cells = sorted(cells)
    
    return lr_pairs, cells

def matrix_to_interaction_tensor(scores: Dict[str, pd.DataFrame], 
                                 cell_delim: str = '-', lr_delim: str = '^', 
                                lr_how: str = 'outer', lr_fill: float = float('nan'),
                                cell_how: str = 'outer', cell_fill: float = 0,
                                prioritize_lr_fill: bool = True):
    """Combine communication matrices from multiple contexts into an interaction tensor.
    
    *Note on 'outer' filling: When taking the union across context-specific matrices, any missing LR or cell pair indicates that this was measured
    in atleast one other context. By default, we treat missing LR pairs as NaN (e.g., missing due to technical limitations) and missing cell pairs
    as 0 (e.g., missing due to being biologically zero). 

    Parameters
    ----------
    scores : Dict[str, pd.DataFrame]
        values are the context label, keys are the cell-cell communication matrix for that context
        matrix columns are sender-receiver cell pairs, matrix rows are ligand&receptor cell pairs, and entries are the non-negative communication score
    cell_delim : str, optional
       delimiter that separates two cell types in the cell type pair for the matrix columns, by default '-'
    lr_delim : str, optional
        delimiter that separates ligands and receptors in the ligand-receptor pair for the matrix indeces, by default '^'
    lr_how : str, optional
        each matrix will contain the union (outer) or intersection ('inner') of all matrix row names, by default 'outer'
    lr_fill : float, optional
        value to fill in missing LR indices with, by default float('nan')
        decomposition will treat NaN as missing (impute or ignore depending on algorithm) and 0 as true biological 0s
        only applicble when lr_how = 'outer'
    cell_how: str, optional
        each matrix will contain the union (outer) or intersection ('inner') of all matrix column names, by default 'outer'
    cell_fill : float, optional
        value to fill in missing cell pair indices with, by default 0
        decomposition will treat NaN as missing (impute or ignore depending on algorithm) and 0 as true biological 0s
        only applicble when cell_how = 'outer'
    prioritize_lr_fill : bool, optional
        prioritize lr_fill (True) or cell_fill (False) value, by default True
        only applicable when cell_how = 'outer', lr_how = 'outer', lf_fill != cell_fill, and the matrix is missing both the cell-cell and LR pair

    Returns
    -------
    tensor : PreBuiltTensor
        interaction tensor restructured from scores' communication matrices 
    """
    scores=OrderedDict(scores)
    lr_pairs, cell_pairs = get_cells_and_lrs(df_list = scores.values(), lr_how=lr_how, cell_how=cell_how)
    cell_order = pd.Series(cell_pairs).str.split(cell_delim, expand = True).pivot(index = 0, columns=1).index.tolist()
    
    if prioritize_lr_fill:
        for context, df in scores.items():
            scores[context] = df.reindex(lr_pairs, fill_value=lr_fill).reindex(cell_pairs, fill_value=cell_fill, axis='columns')
    else:
        for context, df in scores.items():
            scores[context] = df.reindex(lr_pairs, fill_value=lr_fill).reindex(cell_pairs, fill_value=cell_fill, axis='columns')
    
    context_order = list()
    for context, df in tqdm(scores.items()):
        scores[context] = _matrix_to_3d_tensor(df=df, cell_delim=cell_delim, lr_delim=lr_delim)
        context_order.append(context)
    
    # sender, receiver, lr, context
    tensor_4d = np.stack(scores.values())
                        
    if (lr_how == 'outer' and lr_fill == float('nan')) or (cell_how == 'outer' and cell_fill == float('nan')):
        mask = (~np.isnan(np.asarray(tensor_4d))).astype(int)
    else:
        mask = None
    
    tensor = PreBuiltTensor(tensor = tensor_4d, 
                         order_names=[context_order, lr_pairs, cell_order, cell_order], 
                         order_labels=['Samples/Contexts', 'Ligand-Receptor Pairs', 
                                      'Sender Cells', 'Receiver Cells'], 
                         mask = mask)

    return tensor

def test_matrix_to_interaction_tensor(scores, **kwargs):
    
    tensor = matrix_to_interaction_tensor(scores, **kwargs)
    
    correct_score = False
    dims = tensor.tensor.shape
    while not correct_score: # ensures non Nan/0
        context_index = np.random.randint(0, dims[0])
        lr_index = np.random.randint(0, dims[1])
        sender_index = np.random.randint(0, dims[2])
        receiver_index = np.random.randint(0, dims[3])

        context_name = tensor.order_names[0][context_index]
        lr_name = tensor.order_names[1][lr_index]
        sender_name = tensor.order_names[2][sender_index]
        receiver_name = tensor.order_names[3][receiver_index]

        tensor_score = tensor.tensor[(context_index, lr_index, sender_index, receiver_index)]
        if not (tensor_score == 0 or np.isnan(tensor_score)):
            correct_score = True

    assert_almost_equal(tensor_score, 
                  scores[context_name].loc[lr_name, cell_delim.join([sender_name, receiver_name])])

def edgelist_to_communication_matrix(edge_list: pd.DataFrame, 
                                     score_col: str, 
                                     sender_cell_col: str, receiver_cell_col: str,
                                     ligand_col: str, receptor_col: str,
                                     cell_delim: str = '-', lr_delim: str = '^', 
                                    fillna_val: float = 0.0):
    """Converts an edge list of cell-cell communication scores for LR pairs into a communication matrix with columns as (Sender-Receiver) pairs and rows as (Ligand&Receptor) pairs

    Parameters
    ----------
    edge_list : pd.DataFrame
        edge list should contain a communication score for each cell-cell and each LR pair
    score_col : str
        name of column containing communication scores
    sender_cell_col : str
        name of column containing sender cell ID
    receiver_cell_col : str
        name of column containing receiver cell ID
    ligand_col : str
        name of column containing ligand ID
    receptor_col : str
        name of column containing receptor ID
    cell_delim : str, optional
        delimiter by which to join sender<DELIM>receiver cell IDs, by default '-'
    lr_delim : str, optional
        delimiter by which to join ligand<DELIM>receptor gene IDs, by default '^'
    fillna_val : float, optional
        value to fill NA in the communication matrix, by default 0.0
        set to None or float('nan') to avoid filling
        in their output edge list, many communication scoring tools will exclude edges (LR communication scores between a 
        given cell pair) that are expressed due to various thresholds. If these LR pairs received
        a non-zero score for another cell pair, when restructuring into a matrix, the excluded edges will results in a NaN
        value. One may want to instead treat these as a communication score of 0, rather than missing (i.e., NaN).

    Returns
    -------
    cm : pd.DataFrame
        communication matrix with columns as (Sender-Receiver) cell pairs and rows as (Ligand^Receptor) gene pairs
    """
    
    cm = edge_list.copy()
    cm[cell_delim.join(['Sender', 'Receiver'])] = cm[[sender_cell_col, receiver_cell_col]].agg(cell_delim.join, axis=1)
    cm[lr_delim.join(['Ligand', 'Receptor'])] = cm[[ligand_col, receptor_col]].agg(lr_delim.join, axis=1)
    cm = cm.loc[:, [cell_delim.join(['Sender', 'Receiver']), lr_delim.join(['Ligand', 'Receptor']), score_col]] 
                     
    cm = cm.pivot(index = cell_delim.join(['Sender', 'Receiver']), columns = lr_delim.join(['Ligand', 'Receptor']), values=score_col).T
    if fillna_val is not None:
        if cm.isna().all(axis=0).any():
            msg = 'Dataset contains LR pairs that are scored NA for all cell pairs. '
            msg += 'You may want to consider dropping these LR pairs rather than filling NA'
            warnings.warn(msg)

        cm.fillna(value = fillna_val, inplace = True)
    return cm