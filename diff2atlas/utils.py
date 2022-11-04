### --- MISC UTILITY FUNCTIONS --- ###

import os,sys
import scanpy as sc 
import pandas as pd
import numpy as np
import milopy
import scipy
from scipy.sparse import csr_matrix
import anndata


def get_cells_in_nhoods(adata, nhood_ids):
    '''
    Get cells in neighbourhoods of interest, store the number of neighbourhoods for each cell in
    adata.obs['in_nhoods']
    '''
    in_nhoods = np.array(adata.obsm['nhoods'][:,nhood_ids.astype('int')].sum(1))
    adata.obs['in_nhoods'] = in_nhoods
    
def anndata2pseudobulk(adata, group_by, 
                       agg="s", 
                       min_ncells = 10,
                       use_layer = None
                      ):
    '''
    Do pseudo-bulking of raw counts in anndata object
    
    Params:
    ------
    - adata: the anndata object
    - group_by: list of obs columns to use for aggregation
    - agg: "s" for sum (if adata.X are counts), "m" for mean (if adata.X are log-counts)
    - min_ncells: minimum number of cells to keep pseudobulk sample (default: 10)
    - use_layer: which layer to use for pseudobulking (default: None, use adata.X)
    
    Returns:
    -------
    AnnData object of same n_vars as adata, but pseudobulked obs
    '''    
    if use_layer is None:
        X = adata.X.copy()
    else:
        X = adata.layers[use_layer].copy()
        
    ## Make obs for pseudobulk
    pseudobulk_obs = adata.obs[group_by].drop_duplicates()
    pseudobulk_obs = pseudobulk_obs[group_by].astype("str")
    pseudobulk_obs.index = pseudobulk_obs[group_by].agg("-".join, axis=1)
    
    ## Add column to obs assigning cells to pseudobulk samples
    adata.obs[group_by] = adata.obs[group_by].astype("str")
    adata.obs["pseudobulk_sample"] = adata.obs[group_by].agg("-".join, axis=1)
    
    ## Sum counts from same sample
    sample_dummies = pd.get_dummies(adata.obs["pseudobulk_sample"])[pseudobulk_obs.index].values
    sample_dummies = scipy.sparse.csr_matrix(sample_dummies)
    pseudobulk_X = X.T.dot(sample_dummies)
    
    ## Make new anndata object
    pseudobulk_adata = anndata.AnnData(pseudobulk_X.T, obs=pseudobulk_obs, var=adata.var)
    
    ## Add number of cells to obs 
    n_cells = adata.obs.groupby('pseudobulk_sample').count().iloc[:,0]
    n_cells.name = "n_cells"
    pseudobulk_adata.obs = pd.concat([pseudobulk_adata.obs, n_cells], axis=1)
    
    ## Filter obs by number of cells threshold
    pseudobulk_adata = pseudobulk_adata[pseudobulk_adata.obs['n_cells'] >= min_ncells].copy()
    return(pseudobulk_adata)

