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
parser.add_argument("h5ad_file",
                    type=str,
                    help="path to input anndata file")
parser.add_argument("perturb_pop",
                    type=str,
                    help="ID of perturbed population")
parser.add_argument("design",
                    type=str,
                    help="ID of reference design")
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
                    default='/nfs/team205/ed6/data/PBMC_CZI_integration_filtered/tmp/',
                    help="path to working directory")
parser.add_argument("--testing", action='store_true',
                    help="Run in testing mode")
args = parser.parse_args()

h5ad_file = args.h5ad_file
perturb_pop = [args.perturb_pop]
ref_design = args.design
annotation_col = args.annotation_col
batch_col = args.batch_col
query_dataset = args.query_dataset
split_seed = args.split_seed
outdir = args.outpath
if args.testing:
    m_epochs = 10
else:
    m_epochs = 400

# Load data
adata = sc.read_h5ad(outdir + h5ad_file)

# Select ctrl and disease samples
np.random.seed(split_seed)
query_samples = adata.obs['sample_id'][adata.obs['dataset_id']
                                       == query_dataset].unique()
samples_design = np.random.choice((0, 1), len(query_samples))
ctrl_samples = query_samples[samples_design == 1]
query_samples = query_samples[samples_design == 0]

# Simulate dataset groups
adata = simulate_query_reference(
    adata, query_annotation=perturb_pop, annotation_col=annotation_col,
    batch_col='sample_id',
    query_batch=query_samples.tolist(),
    ctrl_batch=ctrl_samples.tolist(),
    perturbation_type='remove'
)
assert check_dataset(adata)

# Save intermediate files
sim_id = f"qPBMC_500cells_demo_perturb_{annotation_col}{clean_pop_name('-'.join(perturb_pop))}_queryBatchdataset_id{query_dataset}_seed{split_seed}"
if not os.path.exists(outdir + sim_id):
    os.mkdir(outdir + sim_id)

# Run embedding and differential abundance testing
if 'X_scVI' in adata.obsm:
    del adata.obsm['X_scVI']

adata.obs['Site'] = adata.obs['donor_id'].str[0:3] ## Only for Stephenson data

if ref_design == 'ACR':
    acr_adata = scArches_milo.scArches_atlas_milo_ctrl(
        adata,
        train_params={'max_epochs': m_epochs},
        outdir=outdir + sim_id,
        harmonize_output=False,
        milo_design = "~Site+is_query"
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
    cr_adata = scArches_milo.scArches_ctrl_milo_ctrl(
        adata,
        train_params={'max_epochs': m_epochs},
        outdir=outdir + sim_id,
        harmonize_output=False,
        milo_design = "~Site+is_query"
    )
    write_milo_adata(cr_adata, outdir + sim_id + '/cr_design.h5ad')
else:
    raise ValueError("design should be one of ACR, AR, CR")
