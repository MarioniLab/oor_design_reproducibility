import os
import numpy as np
import scanpy as sc
from oor_benchmark.api import check_dataset
from oor_benchmark.methods import scArches_milo, scVI_milo
from oor_benchmark.metrics.utils import make_OOR_per_group
from oor_benchmark.metrics.FDR_TPR_FPR import FDR_TPR_FPR
from oor_benchmark.metrics.auprc import auprc
import pandas as pd

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("sim_dir",
                    type=str,
                    help="path to simulation directory")
parser.add_argument("n_controls",
                    type=int,
                    help="number of control donors to use")
parser.add_argument("n_querys",
                    type=int,
                    help="number of query donors to use")
parser.add_argument("random_seed",
                    type=int,
                    help="seed for sampling of controls")
parser.add_argument("--annotation_col",
                    default="cell_type",
                    type=str,
                    help="obs column for annotation")
args = parser.parse_args()


def clean_pop_name(string):
    return(''.join(e if e.isalnum() else '_' for e in string))


def main(sim_dir, n_controls, n_querys, random_seed, annotation_col):

    # Make new dir for outputs
    outdir = os.path.join(sim_dir, f'ctrl_size_analysis_nquery{n_querys}/')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    sim_id = 'ctrl_size{n_controls}_seed{random_seed}'.format(**locals())
    outdir = outdir + sim_id + '/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Load query and control data
    adata_merge = sc.read_h5ad(sim_dir + 'ACR_design.scArches_milo.h5ad')

    # Control for collection site when subsampling controls in Stephenson
    adata_merge.obs["Site"] = adata_merge.obs.donor_id.str[0:4]

    # Subset query donors (with highest frac of cells in OOR population)
    np.random.seed(random_seed)
    oor_ct = adata_merge.obs[annotation_col][adata_merge.obs["OOR_state"] == 1][0]
    n_oor_cells_donor = adata_merge.obs[adata_merge.obs[annotation_col] == oor_ct].value_counts(
        'donor_id')
    n_cells_donor = adata_merge.obs.value_counts('donor_id')
    disease_samples = (n_oor_cells_donor[n_oor_cells_donor > 0] /
                       n_cells_donor[n_oor_cells_donor > 0]).sort_values(ascending=False)[0:n_querys].index

    # Subset control donors
    np.random.seed(random_seed)
    ctrl_samples = adata_merge.obs[adata_merge.obs['is_ctrl']
                                   == 1]['donor_id'].unique()
    ctrl_samples = np.random.choice(ctrl_samples, n_controls, replace=False)

    keep_samples = disease_samples.tolist() + ctrl_samples.tolist()
    adata_merge = adata_merge[adata_merge.obs.donor_id.isin(
        keep_samples)].copy()

    assert check_dataset(adata_merge)

    ## --- CR design joint --- ##
    cr_scvi_adata = adata_merge.copy()
    if 'X_scVI' in cr_scvi_adata.obsm:
        del cr_scvi_adata.obsm['X_scVI']
    cr_scvi_adata = scVI_milo.scVI_ctrl_milo_ctrl(
        cr_scvi_adata,
        train_params={'max_epochs': 400},
        outdir=outdir,
        harmonize_output=True,
        milo_design='~Site+is_query'
    )
    make_OOR_per_group(cr_scvi_adata)
    cr_res_df = pd.concat([FDR_TPR_FPR(cr_scvi_adata), auprc(cr_scvi_adata)], 1)
    cr_res_df['OOR_state'] = oor_ct
    cr_res_df['n_ctrls'] = n_controls
    cr_res_df['random_seed'] = random_seed
    cr_res_df['ref_design'] = "CR"
    cr_res_df['emb_method'] = "scVI"
    cr_res_df.to_csv(outdir + 'CR_scVI_results.csv')

    ## --- CR design with mapping --- ##
    cr_adata = adata_merge.copy()
    if 'X_scVI' in cr_adata.obsm:
        del cr_adata.obsm['X_scVI']
    cr_adata = scArches_milo.scArches_ctrl_milo_ctrl(
        cr_adata,
        train_params={'max_epochs': 400},
        outdir=outdir,
        harmonize_output=True,
        milo_design='~Site+is_query'
    )

    make_OOR_per_group(cr_adata)
    cr_res_df = pd.concat([FDR_TPR_FPR(cr_adata), auprc(cr_adata)], 1)
    cr_res_df['OOR_state'] = oor_ct
    cr_res_df['n_ctrls'] = n_controls
    cr_res_df['random_seed'] = random_seed
    cr_res_df['ref_design'] = "CR"
    cr_res_df['emb_method'] = "scArches"
    cr_res_df.to_csv(outdir + 'CR_scArches_results.csv')

    ## -- ACR design -- ##
    acr_adata = adata_merge.copy()
    k = min([(n_controls + n_querys) * 5, 200])
    sc.pp.neighbors(acr_adata, use_rep="X_scVI", n_neighbors=k)
    scArches_milo.run_milo(acr_adata, "query", 'ctrl',
                           annotation_col=annotation_col,
                           design='~Site+is_query'
                           )

    # Harmonize output
    sample_adata = acr_adata.uns["nhood_adata"].T.copy()
    sample_adata.var["OOR_score"] = sample_adata.var["logFC"].copy()
    sample_adata.var["OOR_signif"] = (
        ((sample_adata.var["SpatialFDR"] < 0.1) &
         (sample_adata.var["logFC"] > 0)).astype(int).copy()
    )
    sample_adata.varm["groups"] = acr_adata.obsm["nhoods"].T
    acr_adata.uns["sample_adata"] = sample_adata.copy()

    make_OOR_per_group(acr_adata)
    acr_res_df = pd.concat([FDR_TPR_FPR(acr_adata), auprc(acr_adata)], 1)
    acr_res_df['OOR_state'] = oor_ct
    acr_res_df['n_ctrls'] = n_controls
    acr_res_df['random_seed'] = random_seed
    acr_res_df['ref_design'] = "ACR"
    acr_res_df['emb_method'] = "scArches"
    acr_res_df.to_csv(outdir + 'ACR_scArches_results.csv')


main(args.sim_dir, args.n_controls, args.n_querys,
     args.random_seed, args.annotation_col)
