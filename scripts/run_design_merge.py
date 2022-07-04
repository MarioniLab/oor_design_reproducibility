import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import anndata

import diff2atlas

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("sim_dir",
                    type=str,
                    help="path to simulation directory")
args = parser.parse_args()


def compute_uncertainties(pac_design_adata, d, ref, sim_dir):
    is_d = pac_design_adata.obs['dataset_group'] == d
    unc_labels = diff2atlas.uncertainty_metrics.weighted_knn_transfer_uncertainty(
        model=sim_dir + f'{ref}_model.weighted_KNN_classifier.pkl',
        query_adata=pac_design_adata[is_d],
        train_labels=pac_design_adata[pac_design_adata.obs['dataset_group']
                                      == ref].obs['cell_type']
    )
    pac_design_adata.obs.loc[is_d, 'uncertainty_labels'] = unc_labels.values

    unc_gex = diff2atlas.uncertainty_metrics.trueVSpred_gex_cosine(
        model=sim_dir + f"model_fit_{d}2{ref}/",
        query_adata=pac_design_adata[is_d],
        n_samples=50
    )
    pac_design_adata.obs.loc[is_d, 'uncertainty_gex'] = unc_gex

    unc_z = diff2atlas.uncertainty_metrics.inference_posterior_distance(
        model=sim_dir + f"model_fit_{d}2{ref}/",
        query_adata=pac_design_adata[is_d],
        n_samples=50
    )
    pac_design_adata.obs.loc[is_d, 'uncertainty_z'] = unc_z


def _merge_design(sim_dir):
    adata_atlas = sc.read_h5ad(sim_dir + 'atlas.h5ad')
    adata_query = sc.read_h5ad(sim_dir + 'query.h5ad')
    adata_ctrl = sc.read_h5ad(sim_dir + 'ctrl.h5ad')

    # perturbation-atlas-control
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

    # perturbation-atlas
    pa_design_adata = pac_design_adata[pac_design_adata.obs['dataset_group'] != 'ctrl'].copy(
    )

    # perturbation-control
    pc_design_adata = anndata.concat([adata_ctrl, adata_query],
                                     label='dataset_group',
                                     keys=['ctrl', 'query'])

    pc_design_adata.layers['counts'] = pc_design_adata.X.copy()

    vae_ctrl = scvi.model.SCVI.load(sim_dir + "/model_ctrl/")
    vae_fit_query2ctrl = scvi.model.SCVI.load(
        sim_dir + "/model_fit_query2ctrl/")

    pc_design_adata.obsm['X_scVI'] = np.vstack([
        vae_ctrl.get_latent_representation(),
        vae_fit_query2ctrl.get_latent_representation()
    ])

    #  Quantify uncertainties

    # train KNN classifiers
    diff2atlas.model_wrappers.train_weighted_knn(
        pac_design_adata[pac_design_adata.obs['dataset_group'] == 'atlas'],
        outfile=sim_dir + 'atlas_model.weighted_KNN_classifier.pkl',
        n_neighbors=100
    )
    diff2atlas.model_wrappers.train_weighted_knn(
        pc_design_adata[pc_design_adata.obs['dataset_group'] == 'ctrl'],
        outfile=sim_dir + 'ctrl_model.weighted_KNN_classifier.pkl',
        n_neighbors=100
    )

    compute_uncertainties(pac_design_adata, d='query',
                          ref='atlas', sim_dir=sim_dir)
    compute_uncertainties(pac_design_adata, d='ctrl',
                          ref='atlas', sim_dir=sim_dir)
    compute_uncertainties(pc_design_adata, d='query',
                          ref='ctrl', sim_dir=sim_dir)
    pa_design_adata.obs = pd.concat([pa_design_adata.obs, pac_design_adata.obs.loc[pac_design_adata.obs['dataset_group'] != 'ctrl', [
                                    'uncertainty_labels', 'uncertainty_gex', 'uncertainty_z']]], 1)

    pac_design_adata = pac_design_adata[pac_design_adata.obs['dataset_group'] != 'atlas'].copy(
    )

    #  Make KNN graph with the same K for all
    n_neighbors = 100
    sc.pp.neighbors(pac_design_adata, use_rep='X_scVI',
                    n_neighbors=n_neighbors)
    sc.pp.neighbors(pa_design_adata, use_rep='X_scVI', n_neighbors=n_neighbors)
    sc.pp.neighbors(pc_design_adata, use_rep='X_scVI', n_neighbors=n_neighbors)
    #  Save anndata objects
    pac_design_adata.write_h5ad(sim_dir+'/pac_design.h5ad')
    pa_design_adata.write_h5ad(sim_dir+'/pa_design.h5ad')
    pc_design_adata.write_h5ad(sim_dir+'/pc_design.h5ad')


_merge_design(args.sim_dir)
