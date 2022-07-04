import pytest
import scanpy as sc
import numpy as np
from scvi.data import heart_cell_atlas_subsampled
from diff2atlas.simulation import *

adata = heart_cell_atlas_subsampled()


def test_query_specific_pop():
    '''
    Test that perturbed population is correctly removed when perturbation_type == 'remove
    '''
    population_obs = 'cell_type'
    perturb_pop = ['Myeloid']
    batch_obs = 'donor'
    query_batch = ['D4', 'D6']
    ctrl_batch = ['D2', 'D3']

    simulate_query_reference(
        adata, query_population=perturb_pop, population_obs=population_obs,
        batch_obs=batch_obs, query_batch=query_batch, ctrl_batch=ctrl_batch,
        perturbation_type='remove'
    )

    # Checks
    assert not any(
        adata.obs[adata.obs[population_obs].isin(perturb_pop)]['is_train'] == 1)
    assert not any(
        adata.obs[adata.obs[population_obs].isin(perturb_pop)]['is_ctrl'] == 1)


def test_ctrl_batch():
    '''Test that no cell belongs to 2 data groups
    '''
    population_obs = 'cell_type'
    perturb_pop = ['Myeloid']
    batch_obs = 'donor'
    query_batch = ['D4', 'D6']
    ctrl_batch = ['D2', 'D3']

    simulate_query_reference(
        adata, query_population=perturb_pop, population_obs=population_obs,
        batch_obs=batch_obs, query_batch=query_batch, ctrl_batch=ctrl_batch
    )

    # Check that there is no cell in 2 splits
    assert adata.obs[['is_train', 'is_ctrl', 'is_test']].sum(1).max() == 1


def test_depletion_expansion():
    '''
    Test that, when perturbation_type 'expansion' or 'depletion', 
    there are cells from perturb population in both query and control
    '''
    population_obs = 'cell_type'
    perturb_pop = ['Myeloid']
    batch_obs = 'donor'
    query_batch = ['D4', 'D6']
    ctrl_batch = ['D2', 'D3']

    simulate_query_reference(
        adata, query_population=perturb_pop, population_obs=population_obs,
        batch_obs=batch_obs, query_batch=query_batch, ctrl_batch=ctrl_batch,
        perturbation_type='expansion'
    )

    assert any(adata.obs[adata.obs[population_obs].isin(
        perturb_pop)]['is_ctrl'] == 1)
    assert any(adata.obs[adata.obs[population_obs].isin(
        perturb_pop)]['is_test'] == 1)


def test_input_errors():
    '''Test that simulation function gives an error if '''
    population_obs = 'cell_type'
    perturb_pop = ['Myeloid']
    batch_obs = 'donor'
    query_batch = 'D4'
    ctrl_batch = ['D2', 'D3']

    with pytest.raises(TypeError):
        simulate_query_reference(
            adata, query_population=perturb_pop, population_obs=population_obs,
            batch_obs=batch_obs, query_batch=query_batch, ctrl_batch=ctrl_batch,
            perturbation_type='expansion'
        )
