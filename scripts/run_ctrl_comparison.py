import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import scvi

import diff2atlas

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("h5ad_file",
                    type=str,
                    help="path to input anndata file")
parser.add_argument("ctrl_class",
                    type=str,
                    help="ID of ctrl class")
parser.add_argument("perturb_pop",
                    type=str,
                    help="ID of perturbed population")
parser.add_argument("--population_obs",
                    default="cell_type",
                    help="column with population info used for simulation")
parser.add_argument("--batch_obs",
                    default="dataset_id",
                    help="column with batch info used for simulation")
parser.add_argument("--outpath",
                    default='/nfs/team205/ed6/data/PBMC_CZI_integration_filtered/',
                    help="path to working directory")
args = parser.parse_args()


def clean_pop_name(string):
    return(''.join(e if e.isalnum() else '_' for e in string))


def _make_design(adata,
                 ctrl_comparison_design,
                 ctrl_class,
                 batch_obs):
    '''
    Assign 2 case/ctrl/atlas
    '''
    query_samples = ctrl_comparison_design.index[ctrl_comparison_design['is_test'] == 1]
    ctrl_samples = ctrl_comparison_design.index[
        ctrl_comparison_design[f'{ctrl_class}_ctrl'] == 1]

    adata.obs['is_test'] = adata.obs['donor_id'].isin(
        query_samples).astype(int)
    adata.obs['is_ctrl'] = adata.obs['donor_id'].isin(ctrl_samples).astype(int)
    adata.obs['is_train'] = (adata.obs['is_ctrl'] +
                             adata.obs['is_test'] == 0).astype('int')
    # Remove matched dataset from atlas
    if ctrl_class != 'matched':
        query_dataset = adata.obs.loc[adata.obs['is_test'] == 1, batch_obs].unique()[
            0]
        adata.obs.loc[(adata.obs[batch_obs] == query_dataset),
                      'is_train'] = 0
    assert all(adata.obs[adata.obs['is_train'] == 1]
               [batch_obs].unique() != query_dataset)


def _train_design(adata,
                  ctrl_class,
                  perturb_pop,
                  population_obs='cell_type',
                  data_dir='/nfs/team205/ed6/data/PBMC_CZI_integration_filtered/',
                  n_hvgs=5000
                  ):

    # Remove query-specific pop from ctrl and atlas
    adata.obs.loc[(adata.obs[population_obs].isin(perturb_pop)),
                  'is_train'] = 0
    adata.obs.loc[(adata.obs[population_obs].isin(
        perturb_pop)), 'is_ctrl'] = 0

    # test that no cell is assigned to multiple splits
    assert adata.obs[['is_train', 'is_test', 'is_ctrl']].sum(1).max() == 1

    # test that perturbed population is in condition dataset only
    assert adata[adata.obs['is_test'] ==
                 1].obs[population_obs].isin(perturb_pop).sum() > 0
    assert adata[adata.obs['is_train'] ==
                 1].obs[population_obs].isin(perturb_pop).sum() == 0
    assert adata[adata.obs['is_ctrl'] == 1].obs[population_obs].isin(
        perturb_pop).sum() == 0

    adata_atlas = adata[adata.obs['is_train'] == 1].copy()
    adata_ctrl = adata[adata.obs['is_ctrl'] == 1].copy()
    adata_query = adata[adata.obs['is_test'] == 1].copy()

    #  Save anndata objects

    #  Make folder to store results
    sim_id = f"ctrl_comparison_PBMC_500cells_demo_perturb_{population_obs}{clean_pop_name('-'.join(perturb_pop))}_ctrl_{ctrl_class}"
    if not os.path.exists(data_dir + sim_id):
        os.mkdir(data_dir + sim_id)
    adata_atlas.write_h5ad(data_dir+sim_id+'/atlas.h5ad')
    adata_query.write_h5ad(data_dir+sim_id+'/query.h5ad')
    adata_ctrl.write_h5ad(data_dir+sim_id+'/ctrl.h5ad')

    adata_atlas.layers['counts'] = adata_atlas.X.copy()
    adata_query.layers['counts'] = adata_query.X.copy()
    adata_ctrl.layers['counts'] = adata_ctrl.X.copy()

    # Train model on atlas
    if 'log1p' not in adata_atlas.uns.keys():
        sc.pp.normalize_per_cell(adata_atlas)
        sc.pp.log1p(adata_atlas)

    sc.pp.highly_variable_genes(
        adata_atlas,
        n_top_genes=n_hvgs,
        subset=True
    )

    vae_ref = diff2atlas.model_wrappers.train_scVI(
        adata_atlas, data_dir + sim_id + "/model_reference", batch_key='sample_id')

    # Map query datasets
    adata_query_fit = adata_query.copy()
    diff2atlas.model_wrappers.fit_scVI(
        vae_ref,
        adata_query_fit,
        outfile=data_dir + sim_id + "/model_fit_query2atlas/")

    adata_ctrl_fit = adata_ctrl.copy()
    diff2atlas.model_wrappers.fit_scVI(
        vae_ref,
        adata_ctrl_fit,
        outfile=data_dir + sim_id + "/model_fit_ctrl2atlas/")

    adata_atlas = sc.read_h5ad(data_dir + sim_id + '/atlas.h5ad')
    adata_query = sc.read_h5ad(data_dir + sim_id + '/query.h5ad')
    adata_ctrl = sc.read_h5ad(data_dir + sim_id + '/ctrl.h5ad')

    # perturbation-atlas-control
    sim_dir = data_dir + sim_id
    pac_design_adata = anndata.concat([adata_atlas, adata_ctrl, adata_query],
                                      label='dataset_group',
                                      keys=['atlas', 'ctrl', 'query'])
    vae_atlas = scvi.model.SCVI.load(sim_dir + "/model_reference/")
    vae_fit_query2atlas = scvi.model.SCVI.load(
        sim_dir + "/model_fit_query2atlas/")
    vae_fit_ctrl2atlas = scvi.model.SCVI.load(
        sim_dir + "/model_fit_ctrl2atlas/")

    pac_design_adata.obsm['X_scVI'] = np.vstack([
        vae_atlas.get_latent_representation(),
        vae_fit_ctrl2atlas.get_latent_representation(),
        vae_fit_query2atlas.get_latent_representation()
    ])

    pac_design_adata = pac_design_adata[pac_design_adata.obs['dataset_group'] != 'atlas'].copy(
    )

    #  Make KNN graph with the same K for all
    n_neighbors = 100
    sc.pp.neighbors(pac_design_adata, use_rep='X_scVI',
                    n_neighbors=n_neighbors)
    #  Save anndata objects
    pac_design_adata.write_h5ad(sim_dir+'/pac_design.h5ad')


adata = sc.read_h5ad(args.outpath + args.h5ad_file)
ctrl_comparison_design = pd.read_csv(
    args.outpath + 'ctrl_comparison_design.csv', index_col=0)

_make_design(adata,
             ctrl_comparison_design,
             ctrl_class=args.ctrl_class,
             batch_obs=args.batch_obs)

_train_design(
    adata,
    ctrl_class=args.ctrl_class,
    perturb_pop=[args.perturb_pop],
    population_obs=args.population_obs,
    data_dir=args.outpath,
)
