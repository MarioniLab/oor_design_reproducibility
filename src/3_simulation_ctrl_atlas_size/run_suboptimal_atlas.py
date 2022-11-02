import argparse
from multiprocessing.sharedctypes import Value
import os
import numpy as np
import scanpy as sc

from oor_benchmark.api import check_dataset
from oor_benchmark.datasets.simulation import simulate_query_reference
from oor_benchmark.methods import scArches_milo
from milopy.utils import write_milo_adata


def clean_pop_name(string):
    return(''.join(e if e.isalnum() else '_' for e in string))


parser = argparse.ArgumentParser()
parser.add_argument("ct_oi",
                    type=str,
                    help="ID of cell type of interest")
parser.add_argument("design",
                    type=str,
                    help="ID of reference design")
parser.add_argument("keep_n_studies",
                    type=int,
                    help="Number of studies to keep (from 1 to 12)")
parser.add_argument("--annotation_col",
                    default="cell_type",
                    help="column with population info used for simulation")
parser.add_argument("--batch_col",
                    default='dataset_id',
                    help="column with batch info used for simulation")
parser.add_argument("--query_dataset",
                    default='10_1038_s41591_021_01329_2',
                    help="ID of query batch")
parser.add_argument("--split_seed",
                    default=2022,
                    help="ID of query batch")
parser.add_argument("--outpath",
                    default='/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/suboptimal_atlas_analysis/',
                    help="path to working directory")
parser.add_argument("--testing", action='store_true',
                    help="Run in testing mode")
parser.add_argument("--oor_ct_oi", action='store_true',
                    help="Run with OOR population corresponding to cell type of interest")
args = parser.parse_args()

ct_oi = [args.ct_oi]
ref_design = args.design
keep_n_studies = args.keep_n_studies
annotation_col = args.annotation_col
batch_col = args.batch_col
query_dataset = args.query_dataset
split_seed = args.split_seed
outdir = args.outpath
if args.testing:
    m_epochs = 10
else:
    m_epochs = 400

adata = sc.read_h5ad(
    '/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/tmp/PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad')

# Select ctrl and disease samples
np.random.seed(split_seed)
query_samples = adata.obs['sample_id'][adata.obs['dataset_id']
                                       == query_dataset].unique()
samples_design = np.random.choice((0, 1), len(query_samples))
ctrl_samples = query_samples[samples_design == 1]
query_samples = query_samples[samples_design == 0]

adata.obs['dataset_group'] = 'atlas'
adata.obs.loc[adata.obs.sample_id.isin(ctrl_samples), 'dataset_group'] = 'ctrl'
adata.obs.loc[adata.obs.sample_id.isin(
    query_samples), 'dataset_group'] = 'query'

# Count number of cells per celltype and dataset in atlas
atlas_adata = adata[adata.obs['dataset_group'] == 'atlas']
ct_size = atlas_adata.obs.groupby(
    ['dataset_id', 'cell_type']).size().reset_index()
ct_size.columns = ['dataset_id', 'cell_type', 'n_cells']
ct_size = ct_size[ct_size.n_cells > 0]

#  Sort study IDs by abundance of cell type of interest
mean_n_cells_ct = ct_size[ct_size.cell_type.isin(ct_oi)].groupby(
    'dataset_id').mean('n_cells').fillna(0).sort_values('n_cells', ascending=False)
order_datasets_exclude = mean_n_cells_ct.index.tolist()

# Exclude studies where CT oi is most abundant
exclude_studies = order_datasets_exclude[:len(
    order_datasets_exclude) - keep_n_studies]
adata = adata[~adata.obs['dataset_id'].isin(exclude_studies)].copy()
assert adata.obs['dataset_id'].nunique() == keep_n_studies + 1
assert adata.obs[adata.obs['dataset_group'] ==
                 'atlas']['dataset_id'].nunique() == keep_n_studies
assert adata.obs[adata.obs['dataset_group'] !=
                 'atlas']['dataset_id'].unique() == query_dataset

if args.oor_ct_oi:
    # Remove OOR populations corresponding to cell type of interest
    adata = simulate_query_reference(
        adata, query_annotation=ct_oi, annotation_col=annotation_col,
        batch_col='sample_id',
        query_batch=query_samples.tolist(),
        ctrl_batch=ctrl_samples.tolist(),
        perturbation_type='remove'
    )
    assert check_dataset(adata)

# Save intermediate files
if args.oor_ct_oi:
    sim_id = f"suboptimal_atlas_queryBatchdataset_id{query_dataset}_seed{split_seed}_keep{keep_n_studies}studies_removeby{ct_oi[0]}_withOOR/"
else:
    sim_id = f"suboptimal_atlas_queryBatchdataset_id{query_dataset}_seed{split_seed}_keep{keep_n_studies}studies_removeby{ct_oi[0]}/"    
if not os.path.exists(outdir + sim_id):
    os.mkdir(outdir + sim_id)

# Run embedding and differential abundance testing
if 'X_scVI' in adata.obsm:
    del adata.obsm['X_scVI']

adata.obs['Site'] = adata.obs['donor_id'].str[0:3]  # Only for Stephenson data
adata.obs['cell_annotation'] = adata.obs['cell_type'].copy()

if ref_design == 'ACR':
    acr_adata = scArches_milo.scArches_atlas_milo_ctrl(
        adata,
        train_params={'max_epochs': m_epochs},
        outdir=outdir + sim_id,
        harmonize_output=False,
        milo_design="~Site+is_query"
    )
    write_milo_adata(acr_adata, outdir + sim_id + '/acr_design.h5ad')

elif ref_design == 'AR':
    ar_adata = scArches_milo.scArches_atlas_milo_atlas(
        adata,
        train_params={'max_epochs': m_epochs},
        outdir=outdir + sim_id,
        harmonize_output=False
    )
    write_milo_adata(ar_adata, outdir + sim_id + '/ar_design.h5ad')

elif ref_design == 'CR':
    acr_adata = scArches_milo.scArches_ctrl_milo_ctrl(
        adata,
        train_params={'max_epochs': m_epochs},
        outdir=outdir + sim_id,
        harmonize_output=False,
        milo_design="~Site+is_query"
    )
    write_milo_adata(acr_adata, outdir + sim_id + '/cr_design.h5ad')

else:
    raise ValueError("design should be one of ACR, AR, CR")
