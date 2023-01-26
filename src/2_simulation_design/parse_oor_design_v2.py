import pandas as pd
import scanpy as sc
import milopy
from oor_benchmark.metrics import FDR_TPR_FPR
from oor_benchmark.metrics import auprc


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("sim_dir",
                    type=str,
                    help="path to simulation directory")
parser.add_argument("--embedding_method",
                    default="scArches",
                    help="Method for latent embedding")
parser.add_argument("--diff_method",
                    default="milo",
                    help="Method for differential analysis")
parser.add_argument("--annotation_col",
                    default="cell_type",
                    type=str,
                    help="obs column for annotation")
args = parser.parse_args()

# Parse args
simdir = args.sim_dir
population_obs = args.annotation_col
embedding_method = args.embedding_method
diff_method = args.diff_method

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

def read_oor_design_output(simdir, ref_design, embedding_method, diff_method, population_obs):
    perturb_pop = simdir.split(population_obs)[1].split('_queryBatch')[0]
    print(f'Reading {perturb_pop}\n')
    if ref_design == 'AR' and embedding_method == 'scArches' and diff_method == 'milo': 
        h5ad_file = simdir + f'/ar_design.h5ad'
    else:
        h5ad_file = simdir + f'/{ref_design}_design.{embedding_method}_{diff_method}.h5ad'
    if diff_method == 'milo':
        try:
            adata = milopy.utils.read_milo_adata(
                h5ad_file, backed=True)
        except:
            print(f'skipping {ref_design} design')
            return(None)
    else:
        try:
            adata = sc.read_h5ad(h5ad_file, backed=True)
        except:
            print(f'skipping {ref_design} design')
            return(None)
    adata.obs['OOR_state_name'] = perturb_pop
    return(adata)

def parse_design(adata, ref_design):
    if 'sample_adata' not in adata.uns:
        harmonize_output(adata)
    perturb_pop = adata.obs['OOR_state_name'].unique()[0]
    tpr_df = FDR_TPR_FPR.FDR_TPR_FPR(adata)
    auprc_df = auprc.auprc(adata, return_curve=True)
    nhoods_df = adata.uns['sample_adata'].var.copy()
    tpr_df['design'] = ref_design
    tpr_df['OOR_state_name'] = perturb_pop
    auprc_df['design'] = ref_design
    auprc_df['OOR_state_name'] = perturb_pop
    nhoods_df['design'] = ref_design
    nhoods_df['OOR_state_name'] = perturb_pop
    return(nhoods_df, tpr_df, auprc_df)


nhoods_df_all = pd.DataFrame()
tpr_df_all = pd.DataFrame()
auprc_df_all = pd.DataFrame()

for d in ['ACR', "CR", "AR"]:
    adata = read_oor_design_output(simdir, d, embedding_method, diff_method, population_obs)
    if adata is not None:
        nhoods_df, tpr_df, auprc_df = parse_design(adata, d)
        nhoods_df_all = pd.concat([nhoods_df_all, nhoods_df])
        tpr_df_all = pd.concat([tpr_df_all, tpr_df])
        auprc_df_all = pd.concat([auprc_df_all, auprc_df])

# print(auprc_df_all.head())
nhoods_df_all.to_csv(simdir + f'/nhoods_obs.{embedding_method}_{diff_method}.csv')
tpr_df_all.to_csv(simdir + f'/TPR_res.{embedding_method}_{diff_method}.csv')
auprc_df_all.to_csv(simdir + f'/AUPRC_res.{embedding_method}_{diff_method}.csv')
