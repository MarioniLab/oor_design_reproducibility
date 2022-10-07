import argparse
from multiprocessing.sharedctypes import Value
import os
import numpy as np
import scanpy as sc


from oor_benchmark.api import check_dataset
from oor_benchmark.methods import scArches_mappingQC

parser = argparse.ArgumentParser()
parser.add_argument("h5ad_file",
                    type=str,
                    help="path to input anndata file")
args = parser.parse_args()

h5ad_file = args.h5ad_file

# Load data
adata = sc.read_h5ad(h5ad_file)
simdir = os.path.dirname(h5ad_file) + '/'

# Extract design
if h5ad_file.split('/')[-1].startswith('cr'):
    embedding_reference = 'ctrl'
if h5ad_file.split('/')[-1].startswith('ar'):
    embedding_reference = 'atlas'
if h5ad_file.split('/')[-1].startswith('acr'):
    embedding_reference = 'atlas'

adata.obs.loc[adata.obs['dataset_group'] !=
              embedding_reference, 'dataset_group'] = 'query'

# Compute label uncertainty
adata = scArches_mappingQC.scArches_mappingQClabels(
    adata,
    embedding_reference=embedding_reference, outdir=simdir)

# Save results
file_prefix = h5ad_file.split('/')[-1].split('.')[0]
adata[adata.obs['dataset_group'] != embedding_reference].obs[["mappingQC_labels", "OOR_state"]].to_csv(
    simdir + file_prefix + '.mappingQC_labels.csv')
