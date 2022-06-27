# import sys
# import os
# import scvi
# import anndata
# import matplotlib
# import seaborn as sns
# from matplotlib import colors
# import matplotlib.pyplot as plt
# from matplotlib import rcParams
# from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
# import scanpy as sc
# import numpy.random as random
# import scipy
# import sklearn
# import time
# from scipy.spatial import cKDTree
# import scvelo
# from scipy.sparse import csr_matrix, issparse
# from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple
from anndata import AnnData
# import torch

from collections import Counter
from sklearn.neighbors import KNeighborsTransformer


## --- KNN classifier uncertainty --- ##
# As implemented in https://github.com/LungCellAtlas/mapping_data_to_the_HLCA/blob/main/scripts/scarches_label_transfer.py


def weighted_knn_trainer(
        train_adata: AnnData,
        train_adata_emb: str,
        n_neighbors: int = 50):
    """Trains a weighted KNN classifier on ``train_adata``.
    Parameters
    ----------
    train_adata: :class:`~anndata.AnnData`
        Annotated dataset to be used to train KNN classifier with ``label_key`` as the target variable.
    train_adata_emb: str
        Name of the obsm layer to be used for calculation of neighbors. If set to "X", anndata.X will be
        used
    n_neighbors: int
        Number of nearest neighbors in KNN classifier.
    """
    print(
        f"Weighted KNN with n_neighbors = {n_neighbors} ... ",
        end="",
    )
    k_neighbors_transformer = KNeighborsTransformer(
        n_neighbors=n_neighbors,
        mode="distance",
        algorithm="brute",
        metric="euclidean",
        n_jobs=-1,
    )
    if train_adata_emb == "X":
        train_emb = train_adata.X
    elif train_adata_emb in train_adata.obsm.keys():
        train_emb = train_adata.obsm[train_adata_emb]
    else:
        raise ValueError(
            "train_adata_emb should be set to either 'X' or the name of the obsm layer to be used!"
        )
    k_neighbors_transformer.fit(train_emb)
    return k_neighbors_transformer


def weighted_knn_transfer(
    query_adata,
    query_adata_emb,
    ref_adata_obs,
    label_keys,
    knn_model,
    threshold=1,
    pred_unknown=False,
    mode="package",
):
    """Annotates ``query_adata`` cells with an input trained weighted KNN classifier.
    Parameters
    ----------
    query_adata: :class:`~anndata.AnnData`
        Annotated dataset to be used to queryate KNN classifier. Embedding to be used
    query_adata_emb: str
        Name of the obsm layer to be used for label transfer. If set to "X",
        query_adata.X will be used
    ref_adata_obs: :class:`pd.DataFrame`
        obs of ref Anndata
    label_keys: str
        Names of the columns to be used as target variables (e.g. cell_type) in ``query_adata``.
    knn_model: :class:`~sklearn.neighbors._graph.KNeighborsTransformer`
        knn model trained on reference adata with weighted_knn_trainer function
    threshold: float
        Threshold of uncertainty used to annotating cells as "Unknown". cells with
        uncertainties higher than this value will be annotated as "Unknown".
        Set to 1 to keep all predictions. This enables one to later on play
        with thresholds.
    pred_unknown: bool
        ``False`` by default. Whether to annotate any cell as "unknown" or not.
        If `False`, ``threshold`` will not be used and each cell will be annotated
        with the label which is the most common in its ``n_neighbors`` nearest cells.
    mode: str
        Has to be one of "paper" or "package". If mode is set to "package",
        uncertainties will be 1 - P(pred_label), otherwise it will be 1 - P(true_label).
    """
    if not type(knn_model) == KNeighborsTransformer:
        raise ValueError(
            "knn_model should be of type sklearn.neighbors._graph.KNeighborsTransformer!"
        )

    if query_adata_emb == "X":
        query_emb = query_adata.X
    elif query_adata_emb in query_adata.obsm.keys():
        query_emb = query_adata.obsm[query_adata_emb]
    else:
        raise ValueError(
            "query_adata_emb should be set to either 'X' or the name of the obsm layer to be used!"
        )
    top_k_distances, top_k_indices = knn_model.kneighbors(X=query_emb)

    stds = np.std(top_k_distances, axis=1)
    stds = (2.0 / stds) ** 2
    stds = stds.reshape(-1, 1)

    top_k_distances_tilda = np.exp(-np.true_divide(top_k_distances, stds))

    weights = top_k_distances_tilda / np.sum(
        top_k_distances_tilda, axis=1, keepdims=True
    )
    cols = ref_adata_obs.columns[ref_adata_obs.columns.str.startswith(
        label_keys)]
    uncertainties = pd.DataFrame(columns=cols, index=query_adata.obs_names)
    pred_labels = pd.DataFrame(columns=cols, index=query_adata.obs_names)
    for i in range(len(weights)):
        for j in cols:
            y_train_labels = ref_adata_obs[j].values
            counter = Counter(y_train_labels[top_k_indices[i]])
            best_label = max(counter, key=counter.get)
            best_prob = weights[i, y_train_labels[top_k_indices[i]]
                                == best_label].sum()

            # Annotating low probability as unknown
            if pred_unknown:
                if best_prob >= threshold:
                    pred_label = best_label
                else:
                    pred_label = "Unknown"
            else:
                pred_label = best_label

            uncertainties.iloc[i][j] = (max(1 - best_prob, 0))
            pred_labels.iloc[i][j] = (pred_label)

    print("finished!")

    return pred_labels, uncertainties


# ## --- KNN conservation --- ##


# def _knn_from_embedding(X_emb, k):
#     if int(''.join(scipy.__version__.split("."))) >= 160:
#         KNN_pred = cKDTree(X_emb).query(x=X_emb, k=k, workers=5)
#     else:
#         KNN_pred = cKDTree(X_emb).query(x=X_emb, k=k, n_jobs=5)
#     KNN_pred_mat = np.zeros(shape=[KNN_pred[1].shape[0], KNN_pred[1].shape[0]])
#     for i in range(KNN_pred[0].shape[0]):
#         KNN_pred_mat[i, KNN_pred[1][i, :]] = 1
#     return(KNN_pred_mat)


# def KNNconservation(adata_query, q2r_emb='X_scVI_train', q_emb='X_scVI_fit', k=50):
#     # Find NNs in query
#     X_emb_true = adata_query.obsm[q2r_emb].copy()
#     KNN_true_mat = _knn_from_embedding(X_emb_true, k=k)

#     # Find NNs in common embedding
#     X_emb_pred = adata_query.obsm[q_emb].copy()
#     KNN_pred_mat = _knn_from_embedding(X_emb_pred, k=k)

#     KNN_conservation = ((KNN_pred_mat == KNN_true_mat)
#                         & (KNN_pred_mat == 1)).sum(1)/k
#     return(KNN_conservation)


# ## --- Model comparison metrics --- ##

# def trueVSpred_gex_cosine(
#     adata: AnnData,
#     pred_layer: str,
#     true_layer: str = 'logcounts',
#     scale: bool = False
# ):
#     '''
#     Compute cosine distance between true and predicted gene expression profile

#     Params:
#     ------
#     - adata: AnnData object
#     - pred_layer: element in adata.layers storing predicted gene expression profile
#     - true_layer: element in adata.layers storing true gene expression profile
#     - scale: should gex profiles be scaled before computing the distance (== calculating correlation)

#     Returns:
#     -------
#     None, modifies adata in place adding `adata.obs['trueVSpred_gex_cosine']`
#     '''
#     X_true = adata.layers[true_layer].copy()
#     if scipy.sparse.issparse(X_true):
#         X_true = X_true.toarray()

#     X_pred = adata.layers[pred_layer].copy()
#     if scipy.sparse.issparse(X_pred):
#         X_pred = X_pred.toarray()

#     if scale:
#         X_pred = sc.pp.scale(X_pred, zero_center=False)
#         X_true = sc.pp.scale(X_true, zero_center=False)

#     cosine_all = sklearn.metrics.pairwise.cosine_distances(X_true, X_pred)
#     adata.obs['trueVSpred_gex_cosine'] = np.diag(cosine_all)


# def trueVSpred_gex_euclidean(
#     adata: AnnData,
#     pred_layer: str,
#     true_layer: str = 'logcounts',
#     scale: bool = False
# ):
#     '''
#     Compute cosine distance between true and predicted gene expression profile

#     Params:
#     ------
#     - adata: AnnData object
#     - pred_layer: element in adata.layers storing predicted gene expression profile
#     - true_layer: element in adata.layers storing true gene expression profile
#     - scale: should gex profiles be scaled before computing the distance

#     Returns:
#     -------
#     None, modifies adata in place adding `adata.obs['trueVSpred_gex_cosine']`
#     '''
#     X_true = adata.layers[true_layer].copy()
#     if scipy.sparse.issparse(X_true):
#         X_true = X_true.toarray()

#     X_pred = adata.layers[pred_layer].copy()
#     if scipy.sparse.issparse(X_pred):
#         X_pred = X_pred.toarray()

#     if scale:
#         X_pred = sc.pp.scale(X_pred, zero_center=False)
#         X_true = sc.pp.scale(X_true, zero_center=False)

#     cosine_all = sklearn.metrics.pairwise.euclidean_distances(X_true, X_pred)
#     adata.obs['trueVSpred_gex_eucl'] = np.diag(cosine_all)


# def trueVSpred_gex_mse(
#     adata: AnnData,
#     pred_layer: str,
#     true_layer: str = 'logcounts',
#     scale: bool = False
# ):
#     '''
#     Compute cosine distance between true and predicted gene expression profile

#     Params:
#     ------
#     - adata: AnnData object
#     - pred_layer: element in adata.layers storing predicted gene expression profile
#     - true_layer: element in adata.layers storing true gene expression profile
#     - scale: should gex profiles be scaled before computing the distance

#     Returns:
#     -------
#     None, modifies adata in place adding `adata.obs['trueVSpred_gex_cosine']`
#     '''
#     X_true = adata.layers[true_layer].copy()
#     if scipy.sparse.issparse(X_true):
#         X_true = X_true.toarray()

#     X_pred = adata.layers[pred_layer].copy()
#     if scipy.sparse.issparse(X_pred):
#         X_pred = X_pred.toarray()

#     if scale:
#         X_pred = sc.pp.scale(X_pred, zero_center=False)
#         X_true = sc.pp.scale(X_true, zero_center=False)

#     mse = np.array([sklearn.metrics.mean_squared_error(
#         X_true[i, :], X_pred[i, :]) for i in range(X_pred.shape[0])])
#     adata.obs['trueVSpred_gex_mse'] = mse


# # def trueVSpred_gex_weigthed_euclidean(
# #     adata: AnnData,
# #     pred_layer: str,
# #     true_layer: str = 'logcounts'
# # ):
# #     '''
# #     Compute weighted euclidean distance between true and predicted gene expression profile

# #     Params:
# #     ------
# #     - adata: AnnData object
# #     - pred_layer: element in adata.layers storing predicted gene expression profile
# #     - true_layer: element in adata.layers storing true gene expression profile

# #     Returns:
# #     -------
# #     None, modifies adata in place adding `adata.obs['trueVSpred_gex_WNN']`
# #     '''
# #     X_true = adata.layers[true_layer].copy()
# #     if scipy.sparse.issparse(X_true):
# #         X_true = X_true.toarray()

# #     X_pred = adata.layers[pred_layer].copy()
# #     if scipy.sparse.issparse(X_pred):
# #         X_pred = X_pred.toarray()

# #     # Distance as seen in Hao et al. 2020
# #     # ...
# #     adata.obs['trueVSpred_gex_weigthed_euclidean'] = trueVSpred_gex_weigthed_euclidean


# def fitVStrain_error(adata_query,
#                      pred_layer_prefix,
#                      dist_metric='cosine'
#                      ):
#     '''
#     Compute difference in reconstruction error between fit and trained models on query

#     Params:
#     ------
#     - adata_query: anndata object of query dataset
#     - pred_layer_prefix: prefix of type of gene expression prediction to use
#     (i.e. from KNN imputation or posterior sampling)
#     - dist_metric: distance metric to use to compare true and predicted gex profiles
#     (one of 'cosine' or 'weigthed_euclidean', default: cosine)

#     Returns:
#     -------
#     None, modifies adata_query in place adding fitVStrain_error in adata_query.obs
#     '''
#     if dist_metric == 'cosine':
#         unc_function = trueVSpred_gex_cosine
#     if dist_metric == 'weigthed_euclidean':
#         unc_function = trueVSpred_gex_WNN
#     else:
#         raise ValueError(
#             f"{dist_metric} is not a recognized distance function.\n Use one of 'cosine' or 'WNN'")
#     true_layer = pred_layer_prefix.split("_")[1]
#     unc_function(
#         adata_query, pred_layer=f'{pred_layer_prefix}_fit', true_layer=true_layer)
#     error_fit = adata_query.obs[f'trueVSpred_gex_{dist_metric}'].copy()
#     unc_function(
#         adata_query, pred_layer=f'{pred_layer_prefix}_train', true_layer=true_layer)
#     error_train = adata_query.obs[f'trueVSpred_gex_{dist_metric}'].copy()
#     adata_query.obs[f'fitVStrain_error_{pred_layer_prefix}_{dist_metric}'] = error_fit - error_train


# def inference_posterior_distance(vae: scvi.model.SCVI,
#                                  adata: AnnData,
#                                  n_samples: int = 50,
#                                  seed: int = 42):
#     '''
#     Compute distance of sampled z positions to inferred (mean) z position
#     (~ uncertainty around the position in the latent space)

#     Params:
#     -------
#     -

#     Returns:
#     --------
#     -
#     '''

#     scvi.settings.seed = seed

#     # Get inferred latent dims
#     Z_dims = vae.get_latent_representation()

#     # Sample latent dims
#     scdl = vae._make_data_loader(adata=adata)

#     for tensors in scdl:
#         inference_kwargs = dict(n_samples=n_samples)
#         inference_outputs, _ = vae.module.forward(
#             tensors=tensors,
#             inference_kwargs=inference_kwargs,
#             compute_loss=False,
#         )

#     # Compute distance 2 samples
#     dmat_all = torch.cdist(Z_dims, inference_outputs['z'], p=2).cpu().numpy()
#     dmat = np.diagonal(dmat_all, axis1=1, axis2=2)

#     adata.obsm['dist_z_samples'] = dmat.T
#     adata.obs['mean_dist_z_samples'] = dmat.mean(0)


## --- old code --- ##

# def _trueVSimputed_cosine(
#     adata_test,
#     use_rep = 'X_scVI_train',
#     k = 50
#     ):
#     '''
#     Compute cosine distance between true gex profile and profile imputed from nearest
#     neighbors in embedding.

#     Uses scvelo.utils.moments for KNN imputation (does this use the real cell profile too?)
#     '''

#     sc.pp.neighbors(adata_test, n_neighbors=k, use_rep=use_rep)

#     emb = use_rep.lstrip("X_")
#     X_true = sc.pp.normalize_total(adata_test, layers=['counts'], inplace=False)['counts']
#     X_true = np.log1p(X_true)
#     adata_test.layers['logcounts'] = X_true.copy()

#     X_knn_mean = get_moments_binary(adata_test, layer='logcounts')
#     X_knn_var = get_moments_binary(adata_test, layer='logcounts', second_order=True)
#     adata_test.layers['pred_logcounts_' + emb] = X_knn_mean.copy()
#     adata_test.layers['pred_logcounts_var_' + emb] = X_knn_var.copy()

#     _compute_gex_distance(adata_test, pred_layer='pred_logcounts_' + emb)

#     adata_test.obs['cosine_' + emb] = adata_test.obs['cosine'].copy()
#     adata_test.obs = adata_test.obs.drop(['cosine'], 1)
