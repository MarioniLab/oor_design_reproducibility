import warnings

import pandas as pd
import scanpy as sc
import scvi
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
                    help="ID of reference design")
args = parser.parse_args()

# unpack arguments
h5ad_path = args.h5ad_path
design = args.design


def run(
    adata: AnnData,
    embedding_reference: str = "atlas",
    sample_col: str = "sample_id",
    annotation_col: str = "cell_type_clean",
    signif_alpha: float = 0.1,
    outdir: str = None,
    harmonize_output: bool = False,
    milo_design: str = "~ is_query",
    **kwargs,
    ):
    
    diff_reference = "ctrl"

    # Subset to datasets of interest
    try:
        assert diff_reference in adata.obs["dataset_group"].unique()
    except AssertionError:
        raise ValueError(f"Differential analysis reference '{diff_reference}' not found in adata.obs['dataset_group']")

    adata = adata[adata.obs.sort_values("dataset_group").index].copy()

    # Latent embedding
    if embedding_reference == 'atlas':
        vae_ref = scvi.model.SCVI.load("/lustre/scratch117/cellgen/team205/ed6/OOR_benchmark_fetal/model_atlas/")
        # Fit query data to scVI model
        adata_query_fit = adata.copy()
        if outdir is not None:
            q_outdir = outdir + f"/model_fit_query2atlas/"
        else:
            q_outdir = None
        vae_q = _fit_scVI(vae_ref, adata_query_fit, outfile=q_outdir, train_params={'max_epochs':400})

        # Get latent embeddings
        X_scVI_q = pd.DataFrame(vae_q.get_latent_representation(), index=vae_q.adata.obs_names)
        adata.obsm["X_scVI"] = X_scVI_q.loc[adata.obs_names].values
    else:
        embedding_scvi(adata, n_hvgs = 5000, outdir=outdir, batch_key="sample_id", train_params={'max_epochs':400})

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

emb_reference_assignment = {
    "ACR":'atlas',
    "CR":'ctrl'
}

diff_reference_assignment = {
    "ACR":'ctrl',
    "CR":'ctrl'   
}

adata = sc.read_h5ad(h5ad_path)
outdir = '/'.join(h5ad_path.split('/')[:-1]) + '/'
adata = run(adata, 
             embedding_reference = emb_reference_assignment[design], 
             outdir = outdir)

write_milo_adata(outdir + '{design}_design.h5ad')
