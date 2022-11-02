import numpy as np
import pandas as pd
from anndata import AnnData
import scanpy as sc
import torch
import scvi
# import celltypist

from pynndescent import NNDescent
import pynndescent
import pickle as pkl

## --- scVI --- ##


def train_scVI(
        train_adata: AnnData,
        outfile: str = None,
        **kwargs):
    '''
    Train scVI model
    Params:
    ------
    - train_adata: AnnData object of training data (already subset to highly variable genes)
        counts should be stored in train_adata.layers['counts']
    - outfile: path to dir to save trained model
    - **kwargs: arguments for scvi.model.SCVI.setup_anndata (specifying batch etc)
    '''

    scvi.model.SCVI.setup_anndata(train_adata, layer='counts', **kwargs)

    arches_params = dict(
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
        n_layers=2,
    )

    model_train = scvi.model.SCVI(
        train_adata,
        **arches_params
    )

    model_train.train()
    if outfile is not None:
        model_train.save(outfile, save_anndata=True, overwrite=True)
    return(model_train)


def fit_scVI(
    model_train: scvi.model.SCVI,
    query_adata: AnnData,
    outfile: str
) -> scvi.model.SCVI:
    '''
    Fit query data to scVI model with scArches 
    '''
    vae_q = scvi.model.SCVI.load_query_data(
        query_adata,
        model_train,
        inplace_subset_query_vars=True
    )

    vae_q.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0))
    if outfile is not None:
        vae_q.save(outfile, save_anndata=True, overwrite=True)
    return(vae_q)


# ## --- CellTypist --- ##

# def train_celltypist(
#         train_adata: AnnData,
#         labels_obs: str,
#         outfile: str = None,
#         **kwargs):
#     '''
#     Train celltypist model
#     Params:
#     ------
#     - train_adata: AnnData object of training data (already subset to highly variable genes)
#         counts should be stored in train_adata.layers['counts']
#     - labels_obs: column in obs storing cell type labels for classification
#     - outfile: path to pkl file to save trained model
#     - **kwargs: arguments for celltypist.train
#     '''

#     # Preprocessing
#     train_adata.X = train_adata.layers['counts'].copy()
#     sc.pp.normalize_total(train_adata, target_sum=10000)
#     sc.pp.log1p(train_adata)

#     # Train
#     ct_model = celltypist.train(
#         train_adata, labels=labels_obs, n_jobs=10, feature_selection=True, **kwargs)

#     # Save
#     if outfile is not None:
#         ct_model.write(outfile)
#     return(ct_model)

## --- KNN classifier --- ##


def train_weighted_knn(
        train_adata: AnnData,
        outfile: str = None,
        use_rep: str = 'X_scVI',
        n_neighbors: int = 50):
    """Trains a weighted KNN classifier on ``train_adata.obsm[use_rep]``.
    Parameters
    ----------
    train_adata: :
        Annotated dataset to be used to train KNN classifier with ``label_key`` as the target variable.
    outfile: path to pkl file to save trained model 
    use_rep: 
        Name of the obsm layer to be used for calculation of neighbors. If set to "X", anndata.X will be
        used (default: X_scVI)
    n_neighbors: 
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
        raise ValueError(
            "use_rep should be set to either 'X' or the name of the obsm layer to be used!"
        )
    k_neighbors_transformer = pynndescent.NNDescent(
        train_emb,
        n_neighbors=n_neighbors,
        metric="euclidean",
        n_jobs=-1,
    )
    k_neighbors_transformer.prepare()

    if outfile is not None:
        with open(outfile, 'wb') as f:
            pkl.dump(k_neighbors_transformer, f)
    return(k_neighbors_transformer)
