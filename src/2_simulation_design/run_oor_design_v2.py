import argparse
from multiprocessing.sharedctypes import Value
import os
import numpy as np
import scanpy as sc


from oor_benchmark.api import check_dataset
from oor_benchmark.datasets.simulation import simulate_query_reference
from oor_benchmark.methods import scArches_milo, scArches_meld, scVI_milo, scArches_cna
from milopy.utils import write_milo_adata

# --- Utils ---
def _write_adata(adata, h5ad_path):
    adata.write_h5ad(h5ad_path)

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
parser.add_argument("--embedding_method",
                    default="scArches",
                    help="Method for latent embedding")
parser.add_argument("--diff_method",
                    default="milo",
                    help="Method for differential analysis")
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

embedding_method = args.embedding_method
diff_method = args.diff_method

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

## Define params based on methods
if diff_method == 'milo':
    if ref_design == 'ACR' or ref_design == 'CR':
        workflow_params = {
            'harmonize_output':False,
            'milo_design':"~Site+is_query"
        }
    elif ref_design == 'AR':
        workflow_params = {
            'harmonize_output':False
        }
    save_output = write_milo_adata
    if embedding_method == 'scArches':    
        run_workflow = scArches_milo.scArches_milo
    elif embedding_method == 'scVI':
        run_workflow = scVI_milo.scVI_milo
    else:
        raise ValueError("Unknown embedding method (must be one of scVI/scArches)")
elif diff_method == 'meld':
    workflow_params = {
        'harmonize_output':True,
    }
    save_output = _write_adata
    run_workflow = scArches_meld.scArches_meld
elif diff_method == 'cna':
    workflow_params = {
        'harmonize_output':True,
    }
    save_output = _write_adata
    run_workflow = scArches_cna.scArches_cna
else:
    raise ValueError("Unknown differential analysis method (must be one of milo/meld/cna)")

emb_reference_assignment = {
    "ACR":'atlas',
    "AR":'atlas',
    "CR":'ctrl',
}

diff_reference_assignment = {
    "ACR":'ctrl',
    "AR":'atlas',
    "CR":'ctrl',    
}

adata = run_workflow(adata, 
             embedding_reference = emb_reference_assignment[ref_design], 
             diff_reference = diff_reference_assignment[ref_design], 
             train_params={'max_epochs': m_epochs}, 
             outdir = outdir+sim_id, 
             **workflow_params)

save_output(adata, outdir + sim_id + f'/{ref_design}_design.{embedding_method}_{diff_method}.h5ad')
