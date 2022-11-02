import argparse
from multiprocessing.sharedctypes import Value
import os
import numpy as np
import scanpy as sc
import anndata

from oor_benchmark.api import check_dataset
from oor_benchmark.methods import scArches_milo
from milopy.utils import write_milo_adata


def clean_pop_name(string):
    return(''.join(e if e.isalnum() else '_' for e in string))


parser = argparse.ArgumentParser()
parser.add_argument("ref_design",
                    type=str,
                    help="reference design")
parser.add_argument("--datadir",
                    default='/nfs/team205/ed6/data/PBMC_COVID/tmp/',
                    help="path to working directory")
args = parser.parse_args()

# unpack arguments
ref_design = args.ref_design
data_dir = args.datadir

# Load data
adata_A = sc.read_h5ad(data_dir + 'PBMC_COVID.subsample500cells.atlas.h5ad')
adata_D = sc.read_h5ad(data_dir + 'PBMC_COVID.subsample500cells.covid.h5ad')
adata_C = sc.read_h5ad(data_dir + 'PBMC_COVID.subsample500cells.ctrl.h5ad')

adata = anndata.concat([adata_dict[d] for d in ref_design])
if ref_design == 'ACR':
    adata = anndata.concat([adata_A, adata_C, adata_D])
elif ref_design == 'CR':
    adata = anndata.concat([adata_C, adata_D])
elif ref_design == 'AR':
    adata = anndata.concat([adata_A, adata_D])

adata.layers['counts'] = adata.X.copy()
assert check_dataset(adata)

# Save intermediate files
simdir = data_dir + f"/model_reference_{ref_design}"
if not os.path.exists(simdir):
    os.mkdir(simdir)

# Run embedding and differential abundance testing
if 'X_scVI' in adata.obsm:
    del adata.obsm['X_scVI']

adata.obs['Site'] = adata.obs['donor_id'].str[0:3]  # Only for Stephenson data

if ref_design == 'ACR':
    acr_adata = scArches_milo.scArches_atlas_milo_ctrl(
        adata,
        train_params={'max_epochs': 400},
        outdir=simdir,
        harmonize_output=False,
        milo_design="~Site+is_query"
    )
    write_milo_adata(acr_adata, simdir + '/acr_design.h5ad')

elif ref_design == 'AR':
    ar_adata = scArches_milo.scArches_atlas_milo_atlas(
        adata,
        train_params={'max_epochs': 400},
        outdir=simdir,
        harmonize_output=False
    )
    write_milo_adata(ar_adata, simdir + '/ar_design.h5ad')
elif ref_design == 'CR':
    cr_adata = scArches_milo.scArches_ctrl_milo_ctrl(
        adata,
        train_params={'max_epochs': 400},
        outdir=simdir,
        harmonize_output=False,
        milo_design="~Site+is_query"
    )
    write_milo_adata(cr_adata, simdir + '/cr_design.h5ad')
else:
    raise ValueError("design should be one of ACR, AR, CR")
