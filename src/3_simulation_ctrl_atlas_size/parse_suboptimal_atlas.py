from oor_benchmark.metrics import FDR_TPR_FPR
from oor_benchmark.metrics import auprc
import milopy
import pandas as pd

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("simdir",
                    type=str,
                    help="path to simulation directory")
parser.add_argument("--outdir",
                    default='/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/suboptimal_atlas_analysis/',
                    help="path to working directory")
args = parser.parse_args()

d = args.simdir
outdir = args.outdir

def harmonize_output(adata, signif_alpha=0.1):
    if adata.isbacked:
        sample_adata = adata.uns['nhood_adata'].to_memory().T
    else:
        sample_adata = adata.uns["nhood_adata"].T.copy()
    sample_adata.var["OOR_score"] = sample_adata.var["logFC"].copy()
    sample_adata.var["OOR_signif"] = (
        ((sample_adata.var["SpatialFDR"] < signif_alpha) &
         (sample_adata.var["logFC"] > 0)).astype(int).copy()
    )
    sample_adata.varm["groups"] = adata.obsm["nhoods"].T
    adata.uns["sample_adata"] = sample_adata.copy()


def parse_design(adata, ct, ref_design):
    harmonize_output(adata)
    tpr_df = FDR_TPR_FPR.FDR_TPR_FPR(adata)
    auprc_df = auprc.auprc(adata, return_curve=True)
    nhoods_df = adata.uns['sample_adata'].var.copy()
    tpr_df['design'] = ref_design
    tpr_df['OOR_state_name'] = ct
    auprc_df['design'] = ref_design
    auprc_df['OOR_state_name'] = ct
    nhoods_df['design'] = ref_design
    nhoods_df['OOR_state_name'] = ct
    return(nhoods_df, tpr_df, auprc_df)


if d.startswith("suboptimal_atlas"):
    n_studies = int(d.split('keep')[1].split('studies')[0])
    ct_oi = d.split('removeby')[-1]
    ar_adata = milopy.utils.read_milo_adata(outdir + d + '/ar_design.h5ad', backed=True)
    atlas_size = sum(ar_adata.obs['dataset_group'] == 'atlas')
    nhoods_df_ar, tpr_df_ar, auprc_df_ar = parse_design(ar_adata, ct_oi, 'AR')
else:
    remove_dataset, ct_oi = d.split('seed2022_remove')[-1].split('_removeby')
    ct_oi = ct_oi.split("_withOOR")[0]

acr_adata = milopy.utils.read_milo_adata(outdir + d + '/acr_design.h5ad', backed=True)
nhoods_df_acr, tpr_df_acr, auprc_df_acr = parse_design(acr_adata, ct_oi, 'ACR')

## store results
if d.startswith("suboptimal_atlas"):
    nhoods_res_df = pd.concat([nhoods_df_acr, nhoods_df_ar])
    tpr_res_df = pd.concat([tpr_df_acr, tpr_df_ar])
    auprc_df = pd.concat([auprc_df_acr, auprc_df_ar])
    nhoods_res_df['n_studies'] = n_studies
    tpr_res_df['n_studies'] = n_studies
    auprc_df['n_studies'] = n_studies
    tpr_res_df['atlas_size'] = atlas_size
    auprc_df['atlas_size'] = atlas_size
else:
    nhoods_res_df = nhoods_df_acr
    tpr_res_df = tpr_df_acr
    auprc_df = auprc_df_acr
    nhoods_res_df['remove_dataset'] = remove_dataset
    tpr_res_df['remove_dataset'] = remove_dataset
    auprc_df['remove_dataset'] = remove_dataset

nhoods_res_df['CT_oi'] = ct_oi
tpr_res_df['CT_oi'] = ct_oi
auprc_df['CT_oi'] = ct_oi


# Save outputs
nhoods_res_df.to_csv(outdir + d + '/nhoods_res_df.csv')
tpr_res_df.to_csv(outdir + d + '/tpr_res_df.csv')
auprc_df.to_csv(outdir + d + '/auprc_df.csv')
