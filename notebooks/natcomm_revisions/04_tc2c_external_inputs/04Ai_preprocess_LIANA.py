#!/usr/bin/env python
# coding: utf-8

# Format the expression matrices from the BALF COVID dataset 
# 
# Erick sent the log(1+CPM) as an aggregated h5ad file, we will separate into expression matrices for each sample. 

# In[1]:


import os
from tqdm import tqdm

import pandas as pd
import numpy as np
import scanpy as sc


expr_files = '/data2/eric/Tensor-Revisions/COVID-19-BALF-log1p.h5ad' # the log(1+CPM) files from Erick
rev_path = '/data3/hratch/tc2c_analyses_1/natcomm_revisions/'


# In[2]:


def format_for_liana(adata_sample):
    """Format adata of a sample for input to LIANA"""
    expr = pd.DataFrame(adata_sample.X.toarray()).T # liana takes log(1+CPM)
    expr.index = adata_sample.var.index.tolist()
    expr.columns = adata_sample.obs.index.tolist()

    # get the cell type annotation
    md = pd.DataFrame(adata_sample.obs['celltype'])
    md.reset_index(inplace = True)
    md.rename(columns = {'index': 'Cell', 'celltype': 'Annotation'}, inplace = True)
    
    return expr, md


# In[ ]:


adata = sc.read_h5ad(expr_files) # log(1+CPM) BALF from Erick
samples = adata.obs['sample'].unique()

for sample in tqdm(samples):
    fp = rev_path + 'interim/tc2c_external_inputs/liana/liana_inputs2/' + sample 

    adata_sample = adata[adata.obs['sample'] == sample] # subset to sample of interest
    expr,md = format_for_liana(adata_sample) # format for input to NATMI

    # files must be written because NATMI is command-line only
    expr.to_csv(fp + '_CPM_expr.csv', index=True, header=True)
    md.to_csv(fp + '_metadata.csv', index=False, header=True)

