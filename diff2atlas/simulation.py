import numpy as np
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple
from anndata import AnnData


def _split_train_test(
    adata: AnnData,
    population_obs: str = 'leiden',
    test_frac: float = 0.2
):

    test_cells = np.empty(shape=[0, 0]).ravel()
    for p in np.unique(adata.obs[population_obs]):
        p_cells = adata.obs_names[adata.obs[population_obs] == p].values
        p_test_cells = np.random.choice(
            p_cells, size=int(np.round(len(p_cells)*test_frac)))
        test_cells = np.hstack([test_cells, p_test_cells])

    train_cells = adata.obs_names[~adata.obs_names.isin(test_cells)]
    return(train_cells, test_cells)


def simulate_query_reference(
        adata: AnnData,
        batch_obs: str = None,
        query_batch: List[str] = None,
        ctrl_batch: List[str] = None,
        population_obs: str = 'leiden',
        query_population: List[str] = None,
        perturbation_type: str = 'remove',
        test_frac: float = 0.2,
        DA_frac: float = 0.2,
        seed=42):
    '''
    Split single-cell dataset in a training and test set (reference/query) 
    where one population will be either absent, depleted or enriched from the 
    reference and control datasets.
    Params:
    ------
    - adata: AnnData
    - batch_obs: column in adata.obs containing sample identity (should correspond to samples used in differential analysis)
    (default: None, cells are split at random)
    - query_batch: list of samples to assign to query group
    - ctrl_batch: list of samples to assign to control group (default: None, no control group is specified)
    - population_obs: column in adata.obs containing cell type population identity
    - query_population: which cell type population should be perturbed (defaults to None, pick one population at random)
    - perturbation_type: one of 'remove', 'expansion' or 'depletion'. If equal to 'remove' (default) the population specified in query_population will be removed
    from the reference and control, 
    if equal to 'expansion' a fraction of the cells in population specified in query_population 
    will be removed from the samples in ctrl_batch (the fraction specified by DA_test)
    if equal to 'depletion' a fraction of the cells in population specified in query_population 
    will be removed from the samples in query_batch (the fraction specified by DA_test)
    - test_frac: fraction of cells in each population to be included in the query group (only used if batch_obs is None)
    - DA_frac: the fraction of cells of query_population to keep in control if perturbation_type is 'expansion', or in query if perturbation_type is 'depletion'
    - seed: random seed for sampling
    Returns:
    --------
    - None, updates adata.obs in place adding columns 
        - is_train, storing assignment of cells to reference (atlas) group
        - is_test, storing assignment of cells to query group
        - (if ctrl_batch is not None) is_ctrl, storing assignment of cells to ctrl group
    '''
    np.random.seed(seed)

    if not isinstance(query_population, list):
        raise TypeError(
            'A list of strings should be passed to query_population')
    if not isinstance(query_batch, list):
        raise TypeError('A list of strings should be passed to query_batch')
    if not isinstance(ctrl_batch, list):
        raise TypeError('A list of strings should be passed to ctrl_batch')

    # Split in query-control-reference
    if batch_obs is not None:
        # split by defined batches
        query = np.array([s in query_batch for s in adata.obs[batch_obs]])
        adata.obs["is_train"] = (~query).astype(int)
        adata.obs["is_test"] = query.astype('int')
        if ctrl_batch is not None:
            ctrl = np.array([s in ctrl_batch for s in adata.obs[batch_obs]])
            adata.obs["is_ctrl"] = ctrl.astype('int')
            adata.obs['is_train'] = adata.obs['is_train'] - \
                adata.obs['is_ctrl']
    else:
        # random split
        train_cells, test_cells = _split_train_test(
            adata, population_obs=population_obs, test_frac=test_frac)
        adata.obs["is_train"] = adata.obs_names.isin(train_cells).astype('int')
        adata.obs["is_test"] = adata.obs_names.isin(test_cells).astype('int')

    # Pick cell population to perturb
    if query_population is None:
        query_population = np.random.choice(
            adata.obs[population_obs].unique(), size=1)

    #  Apply perturbation
    if perturbation_type == 'remove':
        adata.obs.loc[(adata.obs[population_obs].isin(
            query_population)), 'is_train'] = 0
        if ctrl_batch is not None:
            adata.obs.loc[(adata.obs[population_obs].isin(
                query_population)), 'is_ctrl'] = 0

    elif perturbation_type == 'expansion':
        for b in ctrl_batch:
            query_pop_cells = adata.obs_names[(adata.obs[batch_obs] == b) & (
                adata.obs[population_obs].isin(query_population))]
            cells2remove = np.random.choice(query_pop_cells, size=int(
                np.round(len(query_pop_cells)*(1 - DA_frac))))
            adata.obs.loc[cells2remove, 'is_ctrl'] = 0

    elif perturbation_type == 'depletion':
        for b in query_batch:
            query_pop_cells = adata.obs_names[(adata.obs[batch_obs] == b) & (
                adata.obs[population_obs].isin(query_population))]
            cells2remove = np.random.choice(query_pop_cells, size=int(
                np.round(len(query_pop_cells)*(1 - DA_frac))))
            adata.obs.loc[cells2remove, 'is_query'] = 0

    else:
        raise ValueError(
            "perturbation type should be one of 'remove' or 'perturb_pc'")
    adata.uns["perturbation"] = {
        'population_obs': population_obs,
        'batch_obs': batch_obs,
        'query_population': query_population,
        'query_batch': query_batch,
        'ctrl_batch': ctrl_batch,
        'perturbation_type': perturbation_type}
    return(None)
