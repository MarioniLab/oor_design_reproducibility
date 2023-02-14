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

def pseudobulk4DE(simdir, design, 
                    emb_method,
                     k: int = 50,
                     annotation_col = 'cell_annotation',
                     sample_col = 'sample_id'
                    ):
    '''
    Make pseudobulk object from predicted cell type labels to use for DE analysis
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
    
    # load
    print('Loading')
    if design != 'AR':
        acr_adata = sc.read_h5ad(simdir + f'/{design}_design.{emb_method}_milo.h5ad', backed=False)
    else:
        acr_adata = sc.read_h5ad(simdir + '/ar_design.h5ad', backed=False)
        
    ## Train classifier from latent embedding
    try:
        with open(simdir + f"/weighted_KNN_classifier.{design}_{emb_method}.pkl", "rb") as f:
            knn_classifier = pkl.load(f)
        print('KNN model found')
    except:
        print('Training')
        knn_classifier = _train_weighted_knn(
            acr_adata[acr_adata.obs['dataset_group'] == classifier_reference[design]],
            use_rep = 'X_scVI',
            outfile=simdir + f"/weighted_KNN_classifier.{design}_{emb_method}.pkl",
            n_neighbors = k
        )
    
    print('Predicting')
    y_train_labels = acr_adata[acr_adata.obs['dataset_group'] == classifier_reference[design]].obs[annotation_col]
    query_adata = acr_adata[acr_adata.obs['dataset_group'] != classifier_reference[design]]

    # Predict neighbors
    query_emb = query_adata.obsm['X_scVI'].toarray()
    top_k_indices, top_k_distances = knn_classifier.query(query_emb, k=knn_classifier.n_neighbors)

    stds = np.std(top_k_distances, axis=1)
    stds = (2.0 / stds) ** 2
    stds = stds.reshape(-1, 1)

    top_k_distances_tilda = np.exp(-np.true_divide(top_k_distances, stds))

    # Compute probs and labels
    weights = top_k_distances_tilda / np.sum(top_k_distances_tilda, axis=1, keepdims=True)
    uncertainties = pd.DataFrame(columns=["pred_uncertainty"], index=query_adata.obs_names)
    pred_labels = pd.DataFrame(columns=["pred_label"], index=query_adata.obs_names)

    for i in range(len(weights)):
        counter = Counter(y_train_labels[top_k_indices[i]])
        # Here I assume the highest no of neighbors also has the highest probability
        best_label = max(counter, key=counter.get)
        best_prob = weights[i, y_train_labels[top_k_indices[i]] == best_label].sum()

        uncertainties.iloc[i] = max(1 - best_prob, 0)
        pred_labels.iloc[i] = best_label

    # Save predicted labels for DE analysis
    ct_oi = simdir.split("cell_type")[-1].split("_queryBatch")[0]
    acr_adata.obs['OOR_state'] = acr_adata.obs[annotation_col] == ct_oi
    acr_adata.obs[f'DE_{annotation_col}'] = acr_adata.obs[annotation_col].copy()
    acr_adata.obs.loc[pred_labels.index, f'DE_{annotation_col}'] = pred_labels.values 
    acr_adata.obs[f'DE_{annotation_col}'] = acr_adata.obs[f'DE_{annotation_col}'].cat.remove_unused_categories()

    pred_annotations = acr_adata.obs[[f'DE_{annotation_col}', 'dataset_group', annotation_col, "OOR_state"]].copy()
    
    print("Pseudobulking")
    if 'Site' in acr_adata.obs:
        group_cols = [sample_col, 'dataset_group', "Site", f'DE_{annotation_col}']
    else:
        group_cols = [sample_col, 'dataset_group', f'DE_{annotation_col}']
    
    pbulk_adata = diff2atlas.utils.anndata2pseudobulk(
            acr_adata, 
            group_by=group_cols, 
            min_ncells=1, 
            use_layer=None)
    return(pbulk_adata, pred_annotations)

# Pseudobulk for DE analysis
pbulk_de, pred_annotations = pseudobulk4DE(args.simdir, design=args.ref_design, emb_method=args.emb_method)
pbulk_de.write_h5ad(args.simdir + f'/pseudobulk_DE.{args.ref_design}.{args.emb_method}.h5ad')
pred_annotations.to_csv(args.simdir + f'/DE_pred_annotations.{args.ref_design}.{args.emb_method}.csv')