import argparse
from multiprocessing.sharedctypes import Value
import anndata
import numpy as np
import scanpy as sc

from oor_benchmark.api import check_dataset
from oor_benchmark.methods import scArches_milo, scVI_milo
from milopy.utils import write_milo_adata

def clean_pop_name(string):
    return(''.join(e if e.isalnum() else '_' for e in string))

parser = argparse.ArgumentParser()
parser.add_argument("design",
                    type=str,
                    help="ID of reference design (CR, ACR, AR)")
parser.add_argument("method",
                    type=str,
                    help="method for latent embedding (scVI or scArches)")                
parser.add_argument("--outpath",
                    default='/lustre/scratch117/cellgen/team205/ed6/PBMC_COVID_full/',
                    help="path to save files and models")
parser.add_argument("--testing", action='store_true',
                    help="Run in testing mode")
args = parser.parse_args()

ref_design = args.design
method = args.method
outdir = args.outpath
if args.testing:
    m_epochs = 10
else:
    m_epochs = 400

# Load COVID data
adata_dict = {}
adata_dict['query'] = sc.read_h5ad(outdir + 'PBMC_COVID.full.covid.h5ad')
adata_dict['query'].obs['dataset_group'] = 'query'

if 'A' in ref_design: 
    adata_dict['atlas'] = sc.read_h5ad(outdir + 'PBMC_COVID.subsample500cells.atlas.h5ad')
    adata_dict['atlas'].obs['dataset_group'] = 'atlas'
if 'C' in ref_design:
    adata_dict['ctrl'] = sc.read_h5ad(outdir + 'PBMC_COVID.full.ctrl.h5ad')
    adata_dict['ctrl'].obs['dataset_group'] = 'ctrl'

adata = anndata.concat(adata_dict.values(), join='outer') # join outer to keep info on Site
if 'atlas' in adata_dict.keys():
    keep_vars = np.intersect1d(adata_dict['atlas'].var_names, adata_dict['query'].var_names)
    adata = adata[:,keep_vars].copy()

# Run embedding and differential abundance testing
if 'X_scVI' in adata.obsm:
    del adata.obsm['X_scVI']

run_params = {
    'adata': adata,
    'annotation_col' : 'cell_type',
    'train_params' : {'max_epochs': m_epochs},
    'outdir' : outdir,
    'harmonize_output' : False,
    'milo_design' : "~Site+is_query"
}

if args.method == 'scArches':
    if ref_design == 'ACR':
        out_adata = scArches_milo.scArches_atlas_milo_ctrl(**run_params)
    elif ref_design == 'AR':
        out_adata = scArches_milo.scArches_atlas_milo_atlas(**run_params)
    elif ref_design == 'CR':
        out_adata = scArches_milo.scArches_ctrl_milo_ctrl(**run_params)
    else:
        raise ValueError("design should be one of ACR, AR, CR")
elif args.method == 'scVI':
    if ref_design == 'ACR':
        out_adata = scVI_milo.scVI_atlas_milo_ctrl(**run_params)
    elif ref_design == 'AR':
        out_adata = scVI_milo.scVI_atlas_milo_atlas(**run_params)
    elif ref_design == 'CR':
        out_adata = scVI_milo.scVI_ctrl_milo_ctrl(**run_params)
    else:
        raise ValueError("design should be one of ACR, AR, CR")
else:
    raise ValueError("Method should be one of scVI or scArches")

write_milo_adata(out_adata, outdir + f'/PBMC_COVID.full.{ref_design}design_{method}.h5ad')
