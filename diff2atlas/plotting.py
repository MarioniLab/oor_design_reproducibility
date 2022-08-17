from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple
from anndata import AnnData
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt

from . import core


def plot_nhood_confidence_boxplot(adata: AnnData,
                                  obs_names: List[str],
                                  n_cols: int = None,
                                  ylim: List[float] = None) -> None:
    '''
    Plot boxplot to show results of differential confidence analysis.

    Params:
    -------
    - adata: AnnData object (after running test_confidence)
    - obs_names: list of nhood names to plot
    '''
    try:
        nhood_adata = adata.uns['nhood_adata'].copy()
    except:
        raise ValueError(
            "adata.uns['nhood_adata'] not found -- please run nhood_confidence first")

    # Get covariates used for testing
    test_covariate = nhood_adata.uns['test_confidence']['test_covariate']
    ref_level = nhood_adata.uns['test_confidence']['ref_level']
    X_cond = core._get_test_column(
        adata, test_covariate, ref_level)

    if n_cols is None:
        n_cols = len(obs_names)
    n_rows = round(len(obs_names)/n_cols)

    for i, o in enumerate(obs_names):
        plt.subplot(n_rows, n_cols, i+1)
        #  Get confidence statistic
        Y_conf = nhood_adata[o].layers['confidence'].flatten()

        # Get test statistic
        test_stat = nhood_adata[o].obs['confidence_test_statistic'][0]
        stat_type = nhood_adata.uns['test_confidence']['method']

        pl_df = pd.DataFrame([X_cond, Y_conf]).T
        pl_df.columns = ['condition', 'confidence']
        test_covariate_levels = nhood_adata.var[test_covariate].cat.categories
        if len(test_covariate_levels) == 2:
            ctrl_level = test_covariate_levels[test_covariate_levels != ref_level].tolist()[
                0]
        else:
            ctrl_level = 'other'
        pl_df['condition'] = [ref_level if x ==
                              1 else ctrl_level for x in pl_df.condition]
        pl_df['condition'] = pl_df.condition.astype(
            "category").cat.reorder_categories([ctrl_level, ref_level])

        # Plot
        sns.boxplot(data=pl_df, x='condition', y='confidence')
        sns.stripplot(data=pl_df, x='condition', y='confidence', color='black')
        plt.title(f'{o}\n{stat_type} = {np.round(test_stat, 3)}')
        plt.xlabel('')
        if ylim is not None:
            plt.ylim(ylim[0], ylim[1])

    plt.tight_layout()


def plot_nhood_graph(
    adata,
    basis: str = "X_umap",
    use_nhood_index: bool = False,
    min_size=10,
    plot_edges=False,
    nhoods_key: str = None,
    **kwargs
):
    '''
    Visualize DA results on abstracted graph (wrapper around sc.pl.embedding)

    - adata: AnnData object
    - alpha: significance threshold
    - min_logFC: minimum absolute log-Fold Change to show results (default: 0, show all significant neighbourhoods)
    - min_size: minimum size of nodes in visualization (default: 10)
    - plot_edges: boolean indicating if edges for neighbourhood overlaps whould be plotted (default: False)
    - title: plot title (default: 'DA log-Fold Change')
    - **kwargs: other arguments to pass to scanpy.pl.embedding
    '''
    _build_nhood_graph(adata, basis=basis, use_nhood_index=use_nhood_index)
    nhood_adata = adata.uns["nhood_adata"].copy()

    sc.pl.embedding(nhood_adata, "X_nhood_graph",
                    size=adata.uns["nhood_adata"].obs["Nhood_size"]*min_size,
                    edges=plot_edges, neighbors_key="nhood",
                    frameon=False,
                    **kwargs
                    )


def _build_nhood_graph(adata: AnnData,
                       basis: str = "X_umap",
                       nhoods_key: str = None,
                       use_nhood_index: bool = False):
    '''
    Build graph of neighbourhoods used for visualization on embedding
    Params:
    -------
    - adata: AnnData object
    - basis: string indicating the name of the obsm basis to use to use for layout of neighbourhoods (key in `adata.obsm`)
    '''
    if nhoods_key is None:
        nhoods_key = 'nhoods'
    # Add embedding positions
    if use_nhood_index:
        adata.uns["nhood_adata"].obsm["X_nhood_graph"] = adata[adata.obs["nhood_ixs_refined"] == 1].obsm[basis]
    else:
        adata.uns["nhood_adata"].obsm["X_nhood_graph"] = adata.obsm[nhoods_key].T.dot(
            adata.obsm[basis])/adata.obsm[nhoods_key].T.sum(1)
        adata.uns["nhood_adata"].obsm["X_nhood_graph"] = np.array(
            adata.uns["nhood_adata"].obsm["X_nhood_graph"])
    # Add nhood size
    adata.uns["nhood_adata"].obs["Nhood_size"] = np.array(
        adata.obsm[nhoods_key].sum(0)).flatten()
    # Add adjacency graph
    adata.uns["nhood_adata"].obsp["nhood_connectivities"] = adata.obsm[nhoods_key].T.dot(
        adata.obsm[nhoods_key])
    adata.uns["nhood_adata"].uns["nhood"] = {
        "connectivities_key": "nhood_connectivities", "distances_key": ""}


def highly_variable_nhoods(adata: AnnData):
    try:
        nhood_adata = adata.uns['nhood_adata'].copy()
    except:
        raise ValueError(
            "adata.uns['nhood_adata'] not found -- please run nhood_confidence first")

    plt.hist2d(nhood_adata.obs['confidence_mean'], nhood_adata.obs['confidence_var'],
               norm=matplotlib.colors.LogNorm(), bins=100)
    plt.plot(nhood_adata.obsm['confidence_var_lowess'][:, 0],
             nhood_adata.obsm['confidence_var_lowess'][:, 1], color='red')
    plt.xlabel('nhood mean confidence')
    plt.ylabel('nhood confidence var')
    plt.title(f"{sum(nhood_adata.obs['highly_variable'])} HV nhoods")
