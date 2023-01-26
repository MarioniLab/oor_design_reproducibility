import argparse
from multiprocessing.sharedctypes import Value
import scanpy as sc

from oor_benchmark.methods._latent_embedding import embedding_scvi

parser = argparse.ArgumentParser()
parser.add_argument("h5ad_file",
                    type=str,
                    help="path to input anndata file")
parser.add_argument("--outdir",
                    default='/lustre/scratch117/cellgen/team205/ed6/OOR_benchmark_fetal/',
                    help="Output directory")
parser.add_argument("--testing", action='store_true',
                    help="Run in testing mode")
args = parser.parse_args()

# unpack arguments
h5ad_file = args.h5ad_file
outdir = args.outdir

adata = sc.read_h5ad(outdir + h5ad_file)
adata.obs['dataset_group'] = 'atlas'

adata.obs['sample_id'] = adata.obs['donor'].astype("str") + "_" + adata.obs['method'].astype('str')
embedding_scvi(adata, n_hvgs = 5000, outdir=outdir, batch_key="sample_id", train_params={'max_epochs':400})
adata.write_h5ad(outdir + h5ad_file.replace('.h5ad', '_scvi.h5ad'))