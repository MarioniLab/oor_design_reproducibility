### Differential analysis wrappers ###
import graphtools as gt
import meld
import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Union
from numpy.typing import ArrayLike
from anndata import AnnData

import milopy

import sccoda.util.cell_composition_data as scc_dat
import sccoda.util.comp_ana as scc_ana
import sccoda.util.data_visualization as scc_viz
import tensorflow as tf

## ---- DA analysis ---- ##


def run_milo(adata_design: AnnData,
             query_group: str, reference_group: str,
             sample_col: str = 'sample_id',
             annotation_col: str = 'cell_type',
             design: str = '~ is_query'
             ):
    try:
        adata_design.obs['dataset_group']
    except KeyError:
        raise KeyError(
            'adata_design must have an obs column called dataset_group')
    milopy.core.make_nhoods(adata_design, prop=0.1)
    milopy.core.count_nhoods(adata_design, sample_col=sample_col)
    milopy.utils.annotate_nhoods(
        adata_design[adata_design.obs['dataset_group'] == reference_group], annotation_col)
    adata_design.obs['is_query'] = adata_design.obs['dataset_group'] == query_group
    milopy.core.DA_nhoods(adata_design, design=design)


def run_sccoda(adata_design: AnnData,
               condition_col: str = 'dataset_group',
               sample_col: str = 'sample_id',
               annotation_col: str = 'cell_type',
               ref_cell_type: str = '1_CD4_T',
               seed: int = 1234
               ):

    # Prepare fraction dataframe
    frac_by_condition = (
        adata_design.obs
        .groupby([condition_col, sample_col])
        .apply(lambda x: x.value_counts(annotation_col, normalize=False))
        .reset_index(name="n_cells")
        .assign(condition=lambda x: x[condition_col].astype(str))
    )
    frac_pivot = frac_by_condition.pivot(
        index=[sample_col, "condition"],
        columns=annotation_col,
        values="n_cells"
    )
    frac_pivot = frac_pivot.fillna(0).reset_index()

    scc_df = scc_dat.from_pandas(frac_pivot, covariate_columns=[
                                 sample_col, "condition"])

    # Run test
    np.random.seed(seed)
    tf.random.set_seed(seed)
    sccoda_mod = scc_ana.CompositionalAnalysis(
        scc_df,
        formula="condition",
        reference_cell_type=ref_cell_type,
    )
    sccoda_res = sccoda_mod.sample_hmc(num_results=20000)
    return(sccoda_res)


def run_meld(X_red_dim: ArrayLike,
             sample_labels: List[str],
             conditions: List[str],
             k: int = 15):
    '''
    Run MELD
    - X_red_dim: c x d matrix of dimensionality reduction to use for graph construction
    - sample_labels: assignment of cells to samples
    - conditions: vector of condition names
    '''
    # Make graph
    graph = gt.Graph(X_red_dim, knn=int(k))
    # Make MELD object
    meld_op = meld.MELD()
    meld_op.graph = graph
    # Compute density
    meld_fit = meld_op.transform(sample_labels=np.array(sample_labels))

    # Mean density per replicates
    mean_density = pd.DataFrame(
        np.zeros(shape=(meld_fit.shape[0], len(conditions))),
        index=meld_fit.index,
        columns=conditions,
    )

    for c in conditions:
        c_mean = meld_fit.loc[:, [c in x for x in meld_fit.columns]].mean(1)
        mean_density[c] = c_mean

    # From density to likelihood per condition
    likelihoods = meld.utils.normalize_densities(mean_density)
    likelihoods.columns = [col.split("_")[0] for col in likelihoods.columns]
    return(likelihoods)
