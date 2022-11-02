import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import anndata

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("dataset_ID",
                    type=str,
                    help="dataset ID")
parser.add_argument("disease_condition",
                    type=str,
                    help="normal or COVID or lupus")
parser.add_argument("--n_cells_sample", type=int, default=500,
                    help="path to data directory")
parser.add_argument("--data_dir", default='/nfs/team205/ed6/data/PBMC_CZI_integration_filtered/',
                    help="path to data directory")
args = parser.parse_args()

d = args.dataset_ID
disease_condition = args.disease_condition
n_cells_sample = args.n_cells_sample

# Load dataset


def _read_file(normal_sample_obs, d,
               n_cells_sample=500,
               sample_obs_columns=['sex', 'tissue', 'ethnicity', 'disease', 'assay',
                                   'assay_ontology_term_id', 'sample_id', 'donor_id', 'dataset_id', 'development_stage'],
               cell_obs_columns=['cell_type']
               ):
    adata = sc.read_h5ad(metadata_df.loc[d]['file_path'])

    #  Harmonize obs columns
    sample_id_col = metadata_df.loc[d]['sample identifier column']
    if ' + ' in sample_id_col:
        adata.obs['sample_id'] = adata.obs[sample_id_col.split(
            " + ")].astype("str").agg('-'.join, axis=1)
    else:
        adata.obs['sample_id'] = adata.obs[sample_id_col].values
    adata.obs['donor_id'] = adata.obs[metadata_df.loc[d]
                                      ['donor identifier column']]
    adata.obs['dataset_id'] = metadata_df.loc[d]['Dataset ID']

    #  Subsample condition of interest
    s_obs = adata.obs_names[adata.obs['sample_id'].isin(
        normal_sample_obs.index)]
    if n_cells_sample != None:
        obs_df = adata.obs.loc[s_obs]
        obs_df['sample_id'] = obs_df['sample_id'].astype('str')
        s_obs = obs_df[['sample_id']].groupby(
            'sample_id').sample(n_cells_sample).index

    # Store raw counts in adata.X
    if not adata.raw is None:
        adata.X = adata.raw.X.copy()
        adata.var = adata.raw.var.copy()

    # Save essential anndata components
    adata = anndata.AnnData(
        adata[s_obs].X, obs=adata[s_obs].obs[sample_obs_columns + cell_obs_columns], var=adata[s_obs].var)
    return(adata)


# Read metadata
data_dir = '/nfs/team205/ed6/data/PBMC_CZI_integration_filtered/'
metadata_df = pd.read_csv(
    data_dir + 'CZIintegration_workshop_datasets_filtered.csv', index_col=0)
metadata_df.index = metadata_df["Dataset ID"]
sample_obs = pd.read_csv(
    data_dir + f'PBMC_sample_metadata.{disease_condition}.csv', index_col=0)

outdir = data_dir + "tmp/"
if not os.path.exists(outdir):
    os.mkdir(outdir)

## 
adata = _read_file(sample_obs, d, n_cells_sample=n_cells_sample)
adata.write_h5ad(
    outdir + f"{d}.{disease_condition}.subsample{n_cells_sample}cells.h5ad")
