import os
import argparse
from multiprocessing.sharedctypes import Value
import numpy as np
import scanpy as sc
import scipy.spatial

from oor_benchmark.api import check_dataset
from oor_benchmark.methods import scArches_milo
from milopy.utils import read_milo_adata, write_milo_adata

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("sim_dir",
                    type=str,
                    help="path to simulation directory")
parser.add_argument("subsampling_method",
                    type=str,
                    help="ID of subsampling method (either random or closest)")
parser.add_argument("--random_seed",
                    type=int,
                    default=42,
                    help="random sampling seed")
args = parser.parse_args()


def clean_pop_name(string):
    return(''.join(e if e.isalnum() else '_' for e in string))


def get_sample_adata(adata, n_pcs=30):
    if adata.isbacked:
        sample_adata = adata.uns['nhood_adata'].to_memory().T
    else:
        sample_adata = adata.uns["nhood_adata"].T.copy()
    sample_adata.layers['counts'] = sample_adata.X.copy()

    # Normalize counts x donor
    sc.pp.normalize_total(sample_adata, target_sum=10000)
    sc.pp.log1p(sample_adata)

    #  Dim reduction
    sc.pp.pca(sample_adata)

    #  Compute similarity between donors
    X_pca = sample_adata.obsm['X_pca']
    sample_adata.obsp['distance_full'] = scipy.spatial.distance.cdist(
        X_pca[:, 0:n_pcs], X_pca[:, 0:n_pcs])
    return(sample_adata)

# Find set of closest samples to subset atlas


def find_ctrl_samples(
    sample_adata,
    ctrl_samples_all,
    query_samples,
    method='closest',  #  one of closest, random
    n_ctrl: int = 15,
    random_seed=42,
    n_pcs=10
):
    '''
    Params:
    ------
    - sample_adata: sample-level anndata object
    - ctrl_samples_all: all possible ctrl samples
    - query_samples: list or collection of samples from query dataset
    - n_ctrl: number of control samples to pick
    '''

    # Get distances to query samples
    query_distances = scipy.spatial.distance.cdist(
        sample_adata[query_samples, :].obsm['X_pca'][:, 0:n_pcs], sample_adata[ctrl_samples_all, :].obsm['X_pca'][:, 0:n_pcs])

    # Get n_ctrl closest samples
    if method == 'closest':
        n_closest = 0
        k = 1
        ctrl_samples = np.array([])
        while n_closest < n_ctrl:
            closest_samples_k = np.argsort(query_distances, axis=1)[:, 0:k]
            closest_samples = np.setdiff1d(
                np.unique(closest_samples_k), ctrl_samples)
            closest_dist = query_distances[:, closest_samples].min(
                axis=0).argsort()
            ctrl_samples = np.hstack(
                [ctrl_samples, closest_samples[closest_dist][0:(n_ctrl - n_closest)]])
            n_closest = len(ctrl_samples)
            k += 1
        ctrl_samples = np.array(ctrl_samples_all)[ctrl_samples.astype('int')]

    if method == 'random':
        np.random.seed(random_seed)
        ctrl_samples = np.random.choice(
            ctrl_samples_all, n_ctrl, replace=False)

    return(ctrl_samples)


sim_dir = args.sim_dir
subsampling_method = args.subsampling_method
random_seed = args.random_seed

# Load data (Query mapped to atlas)
adata = read_milo_adata(
    sim_dir + '/ar_design.h5ad')
sample_adata = get_sample_adata(adata)

# Get controls from subsampling atlas samples
query_samples = adata.obs.loc[adata.obs['dataset_group']
                              == 'query', 'sample_id'].unique().tolist()
ctrl_samples_all = adata.obs.loc[adata.obs['dataset_group']
                                 != 'query', 'sample_id'].unique().tolist()
ctrl_samples = find_ctrl_samples(sample_adata, ctrl_samples_all=ctrl_samples_all, query_samples=query_samples,
                                 n_ctrl=len(query_samples), method=subsampling_method, random_seed=random_seed)
adata.obs['dataset_group'] = np.where(adata.obs['sample_id'].isin(query_samples), 'query',
                                      np.where(adata.obs['sample_id'].isin(ctrl_samples), 'ctrl',
                                      'atlas'))
assert adata.obs['dataset_group'].value_counts().shape[0] == 3
assert check_dataset(adata)

outdir = sim_dir + f'/ctrl_by_subsampling_{subsampling_method}{random_seed}/'
if not os.path.exists(outdir):
    os.mkdir(outdir)

# Run embedding and differential abundance testing
if 'X_scVI' in adata.obsm:
    del adata.obsm['X_scVI']

# adata.obs['Site'] = adata.obs['donor_id'].str[0:3]  # Only for Stephenson data

acr_adata = scArches_milo.scArches_atlas_milo_ctrl(
    adata,
    train_params={'max_epochs': 300},
    outdir=outdir,
    harmonize_output=False,
    milo_design="~is_query"
)
write_milo_adata(acr_adata, outdir + 'acr_design.h5ad')
