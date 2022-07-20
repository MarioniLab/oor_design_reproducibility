### Train reference models for COVID analysis ### 
import scanpy as sc
import anndata
import numpy as np

import diff2atlas

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("ref_model",
                    type=str,
                    help="path to trained reference model")
parser.add_argument("query_string",
                    type=str,
                    help="string ID for datasets to include in the query (P or PC)")
parser.add_argument("--datadir",
                    default='/lustre/scratch117/cellgen/team205/ed6/PBMC_COVID/',
                    help="path to working directory")
args = parser.parse_args()


ref_model = args.ref_model
query_string = args.query_string
data_dir = args.datadir

#  just for filling missing genes
adata_A = sc.read_h5ad(
    data_dir + 'PBMC_COVID.subsample500cells.atlas.h5ad', backed='r')
adata_P = sc.read_h5ad(data_dir + 'PBMC_COVID.subsample500cells.covid.h5ad')
adata_C = sc.read_h5ad(data_dir + 'PBMC_COVID.subsample500cells.ctrl.h5ad')
adata_dict = {'P': adata_P, 'C': adata_C}

adata_query = anndata.concat([adata_dict[d] for d in query_string])
ref_id = ref_model.split('/')[-1].split('_')[-1]

# Fill to 0 missing genes
missing_genes = adata_A.var_names[~adata_A.var_names.isin(
    adata_query.var_names)]
n_genes = len(missing_genes)
if n_genes > 0:
    empty_X = np.zeros(shape=[adata_query.n_obs, n_genes])
    empty_adata_query = anndata.AnnData(X=empty_X, obs=adata_query.obs)
    empty_adata_query.var_names = missing_genes
    empty_adata_query.var_names.names = ["index"]
    adata_query_filled = anndata.concat(
        [adata_query, empty_adata_query], axis=1)
    # adata_query_filled = adata_query_filled[:,var_names_model].copy()
    adata_query_filled.obs = adata_query.obs.copy()
    adata_query = adata_query_filled.copy()

# Map query data to reference model
adata_query.layers['counts'] = adata_query.X.copy()
diff2atlas.model_wrappers.fit_scVI(
    data_dir + ref_model,
    adata_query,
    outfile=data_dir + f"/model_query_{query_string}_ref{ref_id}")
