import scanpy as sc 
import pandas as pd
import milopy

from oor_benchmark.api import check_dataset
from oor_benchmark.methods._latent_embedding import _train_scVI, _fit_scVI
from oor_benchmark.methods._latent_embedding import embedding_scvi, embedding_scArches
from oor_benchmark.methods.scArches_milo import run_milo

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("simdir", type=str,
                    help="path to simulation directory")
parser.add_argument("design", type=str,
                    help="ID of reference design")
parser.add_argument("--embedding_method", default="scVI",
                    help="Method for latent embedding")
parser.add_argument("--n_hvgs", type=int,default=5000)
parser.add_argument("--fixed_hvgs", action='store_true')
parser.add_argument("--hvg_data", default=None, type=str)
parser.add_argument("--milo_design", default='~Site+is_query')
args = parser.parse_args()

def _get_model_vars(model_dir):
    vars = sc.read(model_dir + '/adata.h5ad', backed=True).var
    return(vars)

def embedding_fixed_hvgs(adata_merge, train_vars, emb_method, ref_dataset=None, outdir=None, train_params=None, **kwargs):
    # Filter genes
    adata_merge.layers["counts"] = adata_merge.X.copy()
    if emb_method == 'scArches':
        assert ref_dataset in adata_merge.obs["dataset_group"].unique().tolist()
        adata_ref_train = adata_merge[adata_merge.obs["dataset_group"] == ref_dataset][:, train_vars].copy()
    elif emb_method == 'scVI':
        dataset_groups = adata_merge.obs["dataset_group"].unique().tolist()
        dataset_groups.sort()
        ref_dataset = "".join(dataset_groups)
        adata_ref_train = adata_merge[:, train_vars].copy()
    else:
        raise NotImplementedError(f'Embedding {emb_method} not implemented')

    # Train scVI model
    if outdir is not None:
        ref_outdir = outdir + f"/model_{ref_dataset}/"
    else:
        ref_outdir = None
    vae_ref = _train_scVI(adata_ref_train, outfile=ref_outdir, train_params=train_params, **kwargs)

    # Fit query data to scVI model
    if emb_method == 'scArches':
        adata_query_fit = adata_merge[adata_merge.obs["dataset_group"] != ref_dataset].copy()
        if outdir is not None:
            q_outdir = outdir + f"/model_fit_query2{ref_dataset}/"
        else:
            q_outdir = None
        vae_q = _fit_scVI(vae_ref, adata_query_fit, train_params=train_params, outfile=q_outdir)

    # Get latent embeddings
        X_scVI_ref = pd.DataFrame(vae_ref.get_latent_representation(), index=vae_ref.adata.obs_names)
        X_scVI_q = pd.DataFrame(vae_q.get_latent_representation(), index=vae_q.adata.obs_names)
        X_scVI = pd.concat([X_scVI_q, X_scVI_ref], axis=0)
    else:
        X_scVI = pd.DataFrame(vae_ref.get_latent_representation(), index=vae_ref.adata.obs_names)
    adata_merge.obsm["X_scVI"] = X_scVI.loc[adata_merge.obs_names].values

def run_hvg_comparison(
    simdir,
    design = 'CR',
    emb_method = 'scVI',
    n_hvgs = 5000,
    fixed_hvgs = True,
    hvg_data=None,
    milo_design = '~Site+is_query'):

    emb_reference_dict = {
        'ACR': 'atlas',
        'CR': 'ctrl',
        'AR':'atlas'
    }

    diff_reference_dict = {
        'ACR': 'ctrl',
        'CR': 'ctrl',
        'AR':'atlas'
        }

    # Load data
    adata = sc.read_h5ad('/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/tmp/PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad')
    adata_sim_obs = sc.read(simdir + '/model_atlasctrlquery/adata.h5ad', backed=False).obs
    adata = adata[adata_sim_obs.index].copy()
    adata.obs = pd.concat([adata.obs, adata_sim_obs[['dataset_group',	'OOR_state',	'cell_annotation','Site']]], 1)

    adata = adata[adata.obs["dataset_group"].isin([emb_reference_dict[design], diff_reference_dict[design], "query"])]
    adata = adata[adata.obs.sort_values("dataset_group").index].copy()

    assert check_dataset(adata)
    if 'X_scVI' in adata.obsm:
        del adata.obsm['X_scVI']

    if fixed_hvgs:
        train_vars = _get_model_vars(simdir + f'/HVG_comparison/model_{n_hvgs}HVG_{hvg_data}').index
        embedding_fixed_hvgs(adata, train_vars, emb_method, ref_dataset=emb_reference_dict[design],  train_params={'max_epochs': 400})
    else:
        if emb_method == 'scVI':
            embedding_scvi(adata, n_hvgs=n_hvgs, batch_key="sample_id", train_params={'max_epochs': 400})
        elif emb_method == 'scArches':
            embedding_scArches(adata, ref_dataset=emb_reference_dict[design], n_hvgs=n_hvgs, batch_key="sample_id", train_params={'max_epochs': 400})

    # remove embedding_reference from anndata if not needed anymore
    if diff_reference_dict[design] != emb_reference_dict[design]:
        adata = adata[adata.obs["dataset_group"] != emb_reference_dict[design]].copy()

    diff_reference = diff_reference_dict[design]
    # Make KNN graph for Milo neigbourhoods
    n_controls = adata[adata.obs["dataset_group"] == diff_reference].obs['sample_id'].unique().shape[0]
    n_querys = adata[adata.obs["dataset_group"] == "query"].obs['sample_id'].unique().shape[0]
    #  Set max to 200 or memory explodes for large datasets
    k = min([(n_controls + n_querys) * 5, 200])
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=k)

    run_milo(adata, "query", diff_reference, sample_col='sample_id', annotation_col='cell_annotation', design=milo_design)

    # Save results
    if not fixed_hvgs:
        hvg_data = 'new'
    milopy.utils.write_milo_adata(adata, f'{simdir}/HVG_comparison/{design}_design.{n_hvgs}HVG_{hvg_data}.{emb_method}_milo.h5ad')


run_hvg_comparison(args.simdir, args.design, args.embedding_method, args.n_hvgs, args.fixed_hvgs, args.hvg_data, args.milo_design)