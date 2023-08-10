import os
import scanpy as sc 
import pandas as pd
import numpy as np

## r2py setup
import rpy2.rinterface_lib.callbacks
import logging
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

import celltypist

from pynndescent import NNDescent
from scipy.sparse import csc_matrix, identity
import pickle as pkl
from collections import Counter
import diff2atlas
from diff2atlas.utils import anndata2pseudobulk

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("simdir",
                    type=str,
                    help="path to input anndata file")
parser.add_argument("ref_design",
                    type=str,
                    help="ID of reference design")
parser.add_argument("emb_method",
                    type=str,
                    help="ID of embedding method (scArches or scVI)")
args = parser.parse_args()

def _train_weighted_knn(train_adata, outfile: str = None, use_rep: str = "X_scVI", n_neighbors: int = 50):
    """Trains a weighted KNN classifier on ``train_adata.obsm[use_rep]``.

    Parameters
    ----------
    train_adata: AnnData
        Annotated dataset to be used to train KNN classifier with ``label_key`` as the target variable.
    outfile: str
        path to pkl file to save trained model
    use_rep: str
        Name of the obsm layer to be used for calculation of neighbors. If set to "X", anndata.X will be
        used (default: X_scVI)
    n_neighbors: int
        Number of nearest neighbors in KNN classifier.
    """
    print(
        f"Weighted KNN with n_neighbors = {n_neighbors} ... ",
        end="",
    )
    if use_rep == "X":
        train_emb = train_adata.X
    elif use_rep in train_adata.obsm.keys():
        train_emb = train_adata.obsm[use_rep]
    else:
        raise ValueError("use_rep should be set to either 'X' or the name of the obsm layer to be used!")
    k_neighbors_transformer = NNDescent(
        train_emb,
        n_neighbors=n_neighbors,
        metric="euclidean",
        n_jobs=-1,
    )
    k_neighbors_transformer.prepare()

    if outfile is not None:
        with open(outfile, "wb") as f:
            pkl.dump(k_neighbors_transformer, f)
    return k_neighbors_transformer

def _pbulk_filter_genes(pbulk_sdata, min_samples = 10, min_counts = 20):
    '''Remove sporadically expressed genes'''
    n_samples = pbulk_sdata.X.copy()
    n_samples[pbulk_sdata.X.nonzero()] = 1
    pbulk_sdata.var['tot_samples'] =np.array(n_samples.sum(0)).ravel()
    pbulk_sdata.var['tot_counts'] = np.array(pbulk_sdata.X.sum(0)).ravel()
    return(pbulk_sdata[:,(pbulk_sdata.var['tot_counts'] > min_counts) & (pbulk_sdata.var['tot_samples'] > min_samples)])

def run_glmGamPoi_DE(pbulk_adata, 
                     design = '~nhood_groups',
                    ref_level = 'in_nhoods_other',
                    contrast = 'nhood_groupsin_nhoods_critical',
                     n_hvgs = 5000):
    '''
    Run R code for DE analysis with glmGamPoi
    
    Params:
    ------
    - pbulk_adata: anndata object with pseudobulked data
    
    '''
    
    ## Define R function
    glmgampoi_str = f'''
        library(SingleCellExperiment)
        library(glmGamPoi)
        library(scran)
        library(scater)

        run_de <- function(args){{
            pbulk_sdata_X <- args[[1]]
            pbulk_sdata_obs <- args[[2]]
            pbulk_sdata_var <- args[[3]]
            n_hvgs <- args[[4]]

            sce <- SingleCellExperiment(assays = list(counts = t(pbulk_sdata_X)), colData = pbulk_sdata_obs)
            sce <- logNormCounts(sce)

            ## Feature selection w scran (just on test cell types)
            dec <- modelGeneVar(sce)
            hvgs <- getTopHVGs(dec, n = n_hvgs)
            subset_sce <- sce[hvgs,]

            ## Fit
            fit <- glm_gp(subset_sce, design = {design}, reference_level = '{ref_level}')

            ## Test 
            de_res <- test_de(fit, contrast = '{contrast}')    
            de_res[,'gene_name'] <- pbulk_sdata_var[hvgs,]['gene_name']
            return(de_res)
        }}'''
    
    import rpy2.robjects.pandas2ri
    import rpy2.robjects.numpy2ri
    from rpy2.robjects.packages import STAP
    
    ## Get anndata components
    pbulk_sdata_X = pbulk_adata.X.toarray().copy()
    pbulk_sdata_obs = pbulk_adata.obs.copy()
    pbulk_sdata_var = pbulk_adata.var.copy()
    
    r_pkg = STAP(glmgampoi_str, "r_pkg")
    # this was needed for the code to run on jhub
    # if you have a different version of rpy2 you may not need these two lines
    rpy2.robjects.pandas2ri.activate()
    rpy2.robjects.numpy2ri.activate()
    
    # PASS OBJECTS INTO FUNCTION
    args = [pbulk_sdata_X, pbulk_sdata_obs, pbulk_sdata_var, n_hvgs]
    de_res_df = r_pkg.run_de(args).to_csvfile('./DE_results.csv')
    de_res_df = pd.read_csv('./DE_results.csv', index_col=0)
    de_res_df.index = de_res_df['gene_name']
    de_res_df.drop('name', 1, inplace=True)
    os.remove('./DE_results.csv')
    return(de_res_df)

def run_DE_benchmark(simdir, design, 
                    emb_method,
                     ct_model,
                     annotation_col = 'cell_annotation',
                     sample_col = 'sample_id',
                     n_hvgs = 10000
                    ):
    '''
    Run DE benchmark
    '''
    classifier_reference = {
        'CR':'ctrl',
        'ACR':'ctrl',
        'AR':'atlas'
    }
    
    diff_reference = {
        'CR':'ctrl',
        'ACR':'ctrl',
        'AR':'atlas'
    }

    
    # Store top markers for each celltype (from celltypist model)
    pbulk_adata = sc.read_h5ad(simdir + f'pseudobulk_DE.{design}.{emb_method}.h5ad')
    ct_markers_df = pd.DataFrame(index=pbulk_adata.var['gene_name'], columns = celltypist_anno_dict.keys())
    for ct in ct_markers_df.columns:
        ct_markers = ct_model.extract_top_markers(celltypist_anno_dict[ct], top_n=50)
        ct_markers_df[ct] = ct_markers_df.index.isin(ct_markers)
    print('Testing')
    ct_oi = simdir.split("cell_type")[-1].split("_queryBatch")[0]

    ## Find affected annotation (where are the OOR cells?)
    pred_annotation = pd.read_csv(simdir +  f'/DE_pred_annotations.{design}.{emb_method}.csv')
    n_oor_cells = pred_annotation['OOR_state'].sum()
    anno_df = pred_annotation[['OOR_state', f'DE_{annotation_col}']].groupby(f'DE_{annotation_col}').sum('OOR_state')
    anno_df['affected_anno'] = anno_df['OOR_state'] > (n_oor_cells*0.05)

    anno_df['FPR'] = np.nan
    anno_df['TPR'] = np.nan 

    # de_ct = 'classical_monocyte'
    for de_ct in anno_df.index:
        try:
#             print(de_ct)
            pbulk_sdata = pbulk_adata[pbulk_adata.obs[f'DE_{annotation_col}'] == de_ct].copy()
            pbulk_sdata = _pbulk_filter_genes(pbulk_sdata).copy()
            if pbulk_sdata.n_vars < 1:
                print(f"Skipping {ct} - no gene passes criteria")
                raise ValueError()
            de_results = run_glmGamPoi_DE(pbulk_sdata, 
                            design='~dataset_group', ref_level=diff_reference[design], 
                            contrast =f'dataset_groupquery',
                            n_hvgs = n_hvgs)
        except ValueError:
            # print(f"Skipping {de_ct} - no gene passes criteria")
            continue
            
        de_results[f'DE_{annotation_col}'] = de_ct

        if anno_df.loc[de_ct]['affected_anno']:
            de_results['is_OORstate_marker'] = de_results['gene_name'].isin(ct_markers_df.index[ct_markers_df[ct_oi]].tolist())
            TP = sum(de_results[de_results['is_OORstate_marker']].adj_pval < 0.1)
            FN = sum(de_results[de_results['is_OORstate_marker']].adj_pval >= 0.1)
            TN = sum(de_results[~de_results['is_OORstate_marker']].adj_pval >= 0.1)
            FP = sum(de_results[~de_results['is_OORstate_marker']].adj_pval < 0.1)
            if TP + FN != 0:
                anno_df.loc[de_ct, 'TPR'] = TP / (TP+FN)
            else:
                anno_df.loc[de_ct, 'TPR'] = 0
        else:
            FP = sum(de_results['adj_pval'] < 0.1)
            TN = sum(de_results['adj_pval'] >= 0.1)
            TP = 0
            FN = 0

        anno_df.loc[de_ct, 'FPR'] = FP / (FP+TN)

    anno_df['design'] = design 
    anno_df['OOR_state_simulation'] = ct_oi
    return(anno_df)
    # return(pbulk_adata)

# Match cell types btw dataset and celltypist labels
celltypist_anno_dict = {
'natural_killer_cell':'NK cells',
 'memory_B_cell':'Memory B cells',
 'central_memory_CD4_positive_alpha_beta_T_cell':'Tem/Effector helper T cells',
 'naive_thymus_derived_CD4_positive_alpha_beta_T_cell':'Tcm/Naive helper T cells',
 'naive_thymus_derived_CD8_positive_alpha_beta_T_cell':'Tcm/Naive cytotoxic T cells',
 'naive_B_cell':'Naive B cells',
 'classical_monocyte':'Classical monocytes',
 'conventional_dendritic_cell':'DC',
'effector_memory_CD8_positive_alpha_beta_T_cell':'Tem/Effector cytotoxic T cells',
 'mucosal_invariant_T_cell':'MAIT cells',
 'plasmacytoid_dendritic_cell':'pDC',
 'platelet':'Megakaryocytes/platelets',
 'plasmablast':'Plasma cells',
 'erythrocyte':'Erythrocytes',
 'CD14_low_CD16_positive_monocyte':'Non-classical monocytes'
}

# Get celltypist model for immune cells
ct_model = celltypist.Model.load('Immune_All_Low.pkl')

# Run DE analysis
de_comparison_results = run_DE_benchmark(args.simdir, design=args.ref_design, emb_method=args.emb_method, ct_model=ct_model)
de_comparison_results.to_csv(f'{args.simdir}/DE_comparison_results.{args.ref_design}_{args.emb_method}.csv')