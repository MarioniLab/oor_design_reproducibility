### Train reference models for COVID analysis ### 
import scanpy as sc
import anndata

import diff2atlas

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("reference_string",
                    type=str,
                    help="string ID for datasets to include in the reference (A, C, PA or PC)")
parser.add_argument("--datadir",
                    default='/nfs/team205/ed6/data/PBMC_COVID/',
                    help="path to working directory")
args = parser.parse_args()


ref_string = args.reference_string
data_dir = args.datadir

adata_A = sc.read_h5ad(data_dir + 'PBMC_COVID.subsample500cells.atlas.h5ad')
adata_P = sc.read_h5ad(data_dir + 'PBMC_COVID.subsample500cells.covid.h5ad')
adata_C = sc.read_h5ad(data_dir + 'PBMC_COVID.subsample500cells.ctrl.h5ad')
adata_dict = {'A': adata_A, 'P': adata_P, 'C': adata_C}

adata_ref = anndata.concat([adata_dict[d] for d in ref_string])

adata_ref.layers['counts'] = adata_ref.X.copy()

# Filter genes not expressed anywhere 
sc.pp.filter_genes(adata_ref, min_cells=1)

# Select HVGs
if 'log1p' not in adata_ref.uns.keys():
    sc.pp.normalize_per_cell(adata_ref)
    sc.pp.log1p(adata_ref)

sc.pp.highly_variable_genes(
    adata_ref,
    n_top_genes=5000,
    subset=True
)

# Train model
vae_ref = diff2atlas.model_wrappers.train_scVI(
    adata_ref, data_dir + f"/model_reference_{ref_string}", batch_key='sample_id')
