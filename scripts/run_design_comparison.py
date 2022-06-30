import os
import numpy as np
import scanpy as sc

import diff2atlas

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("h5ad_file",
                    type=str,
                    help="path to input anndata file")
parser.add_argument("perturb_pop",
                    type=str,
                    help="ID of perturbed population")
parser.add_argument("--population_obs",
                    default="cell_type",
                    help="column with population info used for simulation")
parser.add_argument("--batch_col",
                    default='dataset_id',
                    help="column with batch info used for simulation")
parser.add_argument("--query_batch",
                    default='10_1038_s41591_021_01329_2',
                    help="ID of query batch")
parser.add_argument("--split_seed",
                    default=2022,
                    help="ID of query batch")
parser.add_argument("--outpath",
                    default='/nfs/team205/ed6/data/PBMC_CZI_integration_filtered/',
                    help="path to working directory")
args = parser.parse_args()


def clean_pop_name(string):
    return(''.join(e if e.isalnum() else '_' for e in string))


def _train_design(adata,
                  batch_obs,
                  query_dataset,
                  perturb_pop,
                  population_obs='cell_type',
                  split_seed=2022,
                  data_dir='/nfs/team205/ed6/data/PBMC_CZI_integration_filtered/',
                  n_hvgs=5000
                  ):

    #  Select query batch
    query = np.array([s in query_dataset for s in adata.obs[batch_obs]])
    adata.obs["is_train"] = (~query).astype(int)
    adata.obs["is_test"] = query.astype('int')

    ## Split in case and ctrl
    np.random.seed(split_seed)
    query_samples = adata.obs['sample_id'][adata.obs[batch_obs].isin(
        query_dataset)].unique()
    samples_design = np.random.choice((0, 1), len(query_samples))
    adata.obs['is_ctrl'] = adata.obs['sample_id'].isin(
        query_samples[samples_design == 1]).astype(int)
    adata.obs.loc[adata.obs['is_ctrl'] == 1, 'is_test'] = 0

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
    sim_id = f"qPBMC_500cells_demo_perturb_{population_obs}{clean_pop_name('-'.join(perturb_pop))}_queryBatch{batch_obs}{'-'.join(query_dataset)}_seed{split_seed}"
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

    # Train model on ctrl
    if 'log1p' not in adata_ctrl.uns.keys():
        sc.pp.normalize_per_cell(adata_ctrl)
        sc.pp.log1p(adata_ctrl)

    sc.pp.highly_variable_genes(
        adata_ctrl,
        n_top_genes=n_hvgs,
        subset=True
    )

    vae_ctrl = diff2atlas.model_wrappers.train_scVI(
        adata_ctrl, data_dir + sim_id + "/model_ctrl", batch_key='sample_id')

    # Map query data to ctrl
    adata_query_fit = adata_query.copy()
    diff2atlas.model_wrappers.fit_scVI(
        vae_ctrl,
        adata_query_fit,
        outfile=data_dir + sim_id + "/model_fit_query2ctrl/")


adata = sc.read_h5ad(args.outpath + args.h5ad_file)

_train_design(
    adata,
    perturb_pop=[args.perturb_pop],
    population_obs=args.population_obs,
    batch_obs=args.batch_col,
    query_dataset=[args.query_batch],
    data_dir=args.outpath,
    split_seed=args.split_seed
)
