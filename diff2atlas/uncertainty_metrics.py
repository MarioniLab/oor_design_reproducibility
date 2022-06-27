import numpy as np
import pandas as pd
import scanpy as sc

import scipy
from scipy.sparse import csr_matrix, issparse
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple
from anndata import AnnData
import scvi
import torch
import sklearn.metrics
from collections import Counter
from sklearn.neighbors import KNeighborsTransformer
from pynndescent import NNDescent
import pickle as pkl


## --- KNN classifier uncertainty --- ##
# Adapting implementation in https://github.com/LungCellAtlas/mapping_data_to_the_HLCA/blob/main/scripts/scarches_label_transfer.py


def weighted_knn_transfer_uncertainty(
    model: Union[NNDescent, str],
    query_adata: AnnData,
    train_labels,
    use_rep: str = 'X_scVI',
    return_labels: bool = False,
):
    """Annotates ``query_adata`` cells with an input trained weighted KNN classifier.
    Parameters
    ----------
    model: :class:`~pynndescent._neighbors.NNDescent`
        knn model trained on reference adata with models.weighted_knn_trainer function
    query_adata: :class:`~anndata.AnnData`
        Annotated dataset to be used to queryate KNN classifier. Embedding to be used
    train_labels: list or series or array of labels from training data 
    use_rep: str
        Name of the obsm layer to be used for label transfer. If set to "X",
        query_adata.X will be used
    ref_adata_obs: :class:`pd.DataFrame`
        obs of ref Anndata
    label_keys: str
        Names of the columns to be used as target variables (e.g. cell_type) in ``query_adata``.
    """
    if type(model) == NNDescent:
        knn_model = model
    elif type(model) == str:
        try:
            with open(model, 'rb') as f:
                knn_model = pkl.load(f)
        except:
            raise FileNotFoundError(
                f'{model} should be either a trained NNDescent object or a path to a pickle file')

    if isinstance(train_labels, pd.Series):
        y_train_labels = train_labels.values
    else:
        y_train_labels = train_labels

    if use_rep == "X":
        query_emb = query_adata.X
    elif use_rep in query_adata.obsm.keys():
        query_emb = query_adata.obsm[use_rep]
    else:
        raise ValueError(
            "use_rep should be set to either 'X' or the name of the obsm layer to be used!"
        )
    top_k_indices, top_k_distances = knn_model.query(
        query_emb, k=knn_model.n_neighbors)

    stds = np.std(top_k_distances, axis=1)
    stds = (2.0 / stds) ** 2
    stds = stds.reshape(-1, 1)

    top_k_distances_tilda = np.exp(-np.true_divide(top_k_distances, stds))

    weights = top_k_distances_tilda / np.sum(
        top_k_distances_tilda, axis=1, keepdims=True
    )
    uncertainties = pd.DataFrame(
        columns=['pred_uncertainty'], index=query_adata.obs_names)
    pred_labels = pd.DataFrame(
        columns=['pred_label'], index=query_adata.obs_names)

    for i in range(len(weights)):
        counter = Counter(y_train_labels[top_k_indices[i]])
        # Here I assume the highest no of neighbors also has the highest probability
        best_label = max(counter, key=counter.get)
        best_prob = weights[i, y_train_labels[top_k_indices[i]]
                            == best_label].sum()

        uncertainties.iloc[i] = (max(1 - best_prob, 0))
        pred_labels.iloc[i] = (best_label)

    if return_labels:
        return(uncertainties, pred_labels)
    else:
        return(uncertainties)


def trueVSpred_gex_cosine(
    model: Union[scvi.model.SCVI, str],
    query_adata: AnnData,
    n_samples=50,
    seed: int = 42,
    scale: bool = False
):
    '''
    Compute cosine distance between true and predicted gene expression profile

    Params:
    ------
    - adata: AnnData object
    - pred_layer: element in adata.layers storing predicted gene expression profile
    - true_layer: element in adata.layers storing true gene expression profile
    - scale: should gex profiles be scaled before computing the distance (== calculating correlation)    

    Returns:
    -------
    None, modifies adata in place adding `adata.obs['trueVSpred_gex_cosine']`
    '''
    if type(model) == scvi.model.SCVI:
        vae = model
    elif type(model) == str:
        try:
            vae = scvi.model.SCVI.load(model)
        except:
            raise FileNotFoundError(
                f'{model} should be either a trained scvi.model.SCVI object or a path to a model dir')

    if 'log1p' not in query_adata.uns.keys():
        sc.pp.normalize_per_cell(query_adata)
        sc.pp.log1p(query_adata)
    X_true = query_adata[:, vae.adata.var_names].X.copy()

    scvi.settings.seed = seed
    post_sample = vae.posterior_predictive_sample(n_samples=n_samples)
    post_sample = np.log1p(post_sample)
    X_pred = post_sample.mean(2)

    if scipy.sparse.issparse(X_true):
        X_true = X_true.toarray()
    if scipy.sparse.issparse(X_pred):
        X_pred = X_pred.toarray()

    if scale:
        X_pred = sc.pp.scale(X_pred, zero_center=False)
        X_true = sc.pp.scale(X_true, zero_center=False)

    cosine_all = sklearn.metrics.pairwise.cosine_distances(X_true, X_pred)
    return(np.diag(cosine_all))


def inference_posterior_distance(vae: scvi.model.SCVI,
                                 adata: AnnData,
                                 n_samples: int = 50,
                                 seed: int = 42):
    '''
    Compute distance of sampled z positions to inferred (mean) z position
    (~ uncertainty around the position in the latent space)

    Params:
    -------
    -

    Returns:
    --------
    -
    '''

    scvi.settings.seed = seed

    # Get inferred latent dims
    Z_dims = vae.get_latent_representation()

    # Sample latent dims
    scdl = vae._make_data_loader(adata=adata)

    for tensors in scdl:
        inference_kwargs = dict(n_samples=n_samples)
        inference_outputs, _ = vae.module.forward(
            tensors=tensors,
            inference_kwargs=inference_kwargs,
            compute_loss=False,
        )

    # Compute distance 2 samples
    dmat_all = torch.cdist(Z_dims, inference_outputs['z'], p=2).cpu().numpy()
    dmat = np.diagonal(dmat_all, axis1=1, axis2=2)

    adata.obsm['dist_z_samples'] = dmat.T
    adata.obs['mean_dist_z_samples'] = dmat.mean(0)


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
