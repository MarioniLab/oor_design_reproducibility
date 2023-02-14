import warnings

import pandas as pd
import scanpy as sc
import scvi
from scvi.hub import HubModel
from scvi.model.utils import mde
from anndata import AnnData

from oor_benchmark.methods.scArches_milo import run_milo
from oor_benchmark.methods._latent_embedding import embedding_scvi, _fit_scVI
from milopy.utils import write_milo_adata
# from ._latent_embedding import embedding_scArches

# logger = logging.getLogger(__name__)
#  Turn off deprecation warnings in scvi-tools
warnings.filterwarnings("ignore", category=DeprecationWarning)

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("h5ad_path",
                    type=str,
                    help="path to input anndata file")
parser.add_argument("design",
                    type=str,
                    help="ID of reference design (either CR, HLCA, of TabulaSapiens)")
args = parser.parse_args()

# unpack arguments
h5ad_path = args.h5ad_path
design = args.design


def run(
    adata: AnnData,
    reference: str = "atlas",
    sample_col: str = "sample_id",
    annotation_col: str = "cell_type_clean",
    signif_alpha: float = 0.1,
    outdir: str = None,
    harmonize_output: bool = False,
    milo_design: str = "~is_query",
    **kwargs,
    ):

    surgery_epochs = 400
    train_kwargs_surgery = {
            "early_stopping": True,
            "early_stopping_monitor": "elbo_train",
            "early_stopping_patience": 10,
            "early_stopping_min_delta": 0.001,
            "plan_kwargs": {"weight_decay": 0.0},
        }
        
    diff_reference = "ctrl"

    # Subset to datasets of interest
    try:
        assert diff_reference in adata.obs["dataset_group"].unique()
    except AssertionError:
        raise ValueError(f"Differential analysis reference '{diff_reference}' not found in adata.obs['dataset_group']")

    adata = adata[adata.obs.sort_values("dataset_group").index].copy()
    adata.obs['dataset'] = adata.obs['Subject_Identity'].copy()

    # Latent embedding
    if reference == 'HLCA':
        model = scvi.model.SCANVI.load('/lustre/scratch117/cellgen/team205/ed6/OOR_benchmark_lung/HLCA_model_huggingface')

        adata_query = adata.copy()
        scvi.model.SCANVI.prepare_query_anndata(adata_query, model)
        adata_query.obs["scanvi_label"] = "unlabeled"
        query_model = scvi.model.SCANVI.load_query_data(adata_query, model)

        query_model.train(max_epochs=surgery_epochs, **train_kwargs_surgery)
        query_model.save('/lustre/scratch117/cellgen/team205/ed6/OOR_benchmark_lung/model_fit_query2hlca', overwrite=True)
        adata.obsm['X_scVI'] = query_model.get_latent_representation()
    elif reference == 'TabulaSapiens':
        model = scvi.model.SCVI.load('/lustre/scratch117/cellgen/team205/ed6/OOR_benchmark_lung/model_atlas')
        model.adata.var['ensID_train'] = [x[0] for x in model.adata.var_names.str.split('.')]
        model.adata.var_names = model.adata.var['ensID_train'].copy()
        
        adata_query = adata.copy()
        scvi.model.SCVI.prepare_query_anndata(adata_query, model)
        adata_query.obs["scanvi_label"] = "unlabeled"
        adata_query.obs["donor"] = adata_query.obs["Subject_Identity"].copy()
        query_model = scvi.model.SCVI.load_query_data(adata_query, model)

        query_model.train(max_epochs=surgery_epochs, **train_kwargs_surgery)
        query_model.save('/lustre/scratch117/cellgen/team205/ed6/OOR_benchmark_lung/model_fit_query2ts', overwrite=True)
        adata.obsm['X_scVI'] = query_model.get_latent_representation()
    else:
        embedding_scvi(adata, n_hvgs = 2000, outdir=outdir, batch_key="dataset", train_params={'max_epochs':400})

    # Make KNN graph for Milo neigbourhoods
    n_controls = adata[adata.obs["dataset_group"] == diff_reference].obs[sample_col].unique().shape[0]
    n_querys = adata[adata.obs["dataset_group"] == "query"].obs[sample_col].unique().shape[0]
    #  Set max to 200 or memory explodes for large datasets
    k = min([(n_controls + n_querys) * 5, 200])
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=k)

    run_milo(adata, "query", diff_reference, sample_col=sample_col, annotation_col=annotation_col, design=milo_design)

    # Harmonize output
    if harmonize_output:
        sample_adata = adata.uns["nhood_adata"].T.copy()
        sample_adata.var["OOR_score"] = sample_adata.var["logFC"].copy()
        sample_adata.var["OOR_signif"] = (
            ((sample_adata.var["SpatialFDR"] < signif_alpha) & (sample_adata.var["logFC"] > 0)).astype(int).copy()
        )
        sample_adata.varm["groups"] = adata.obsm["nhoods"].T
        adata.uns["sample_adata"] = sample_adata.copy()
    return adata

# emb_reference_assignment = {
#     "ACR":'atlas',
#     "CR":'ctrl'
# }

adata = sc.read_h5ad(h5ad_path)
outdir = '/'.join(h5ad_path.split('/')[:-1]) + '/'
adata.layers['counts'] = adata.X.copy()
adata = run(adata, 
             reference = design, 
             outdir = outdir)

write_milo_adata(adata, outdir + f'/{design}_design.h5ad')
