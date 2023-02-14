import argparse
from multiprocessing.sharedctypes import Value
import scanpy as sc

from oor_benchmark.methods._latent_embedding import embedding_scvi

parser = argparse.ArgumentParser()
parser.add_argument("h5ad_file",
                    type=str,
                    help="path to input anndata file")
parser.add_argument("--batch_key",
                    default='sample_id',
                    help="Batch key for scVI")
parser.add_argument("--n_hvgs", type=int,
                    default=5000,
                    help="no of HVGs")
parser.add_argument("--outdir",
                    default='/lustre/scratch117/cellgen/team205/ed6/OOR_benchmark_fetal/',
                    help="Output directory")
parser.add_argument("--testing", action='store_true',
                    help="Run in testing mode")
args = parser.parse_args()

# unpack arguments
h5ad_file = args.h5ad_file
batch_key = args.batch_key
n_hvgs = args.n_hvgs
outdir = args.outdir

adata = sc.read_h5ad(outdir + h5ad_file)
adata.obs['dataset_group'] = 'atlas'

embedding_scvi(adata, n_hvgs = n_hvgs, outdir=outdir, batch_key=batch_key, train_params={'max_epochs':400})
adata.write_h5ad(outdir + h5ad_file.replace('.h5ad', '_scvi.h5ad'))