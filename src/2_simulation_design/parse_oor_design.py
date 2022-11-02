import pandas as pd
import milopy
from oor_benchmark.metrics import FDR_TPR_FPR
from oor_benchmark.metrics import auprc


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("sim_dir",
                    type=str,
                    help="path to simulation directory")
parser.add_argument("--annotation_col",
                    default="cell_type",
                    type=str,
                    help="obs column for annotation")
args = parser.parse_args()

# Parse args
simdir = args.sim_dir
population_obs = args.annotation_col


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


def parse_design(simdir, ref_design, population_obs):
    perturb_pop = simdir.split(population_obs)[1].split('_queryBatch')[0]
    print(f'Reading {perturb_pop}\n')
    try:
        adata = milopy.utils.read_milo_adata(
            simdir + f'/{ref_design.lower()}_design.h5ad', backed=True)
    except:
        print(f'skipping {ref_design} design')
        return(None)

    harmonize_output(adata)
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
    nhoods_df, tpr_df, auprc_df = parse_design(simdir, d, population_obs)
    nhoods_df_all = pd.concat([nhoods_df_all, nhoods_df])
    tpr_df_all = pd.concat([tpr_df_all, tpr_df])
    auprc_df_all = pd.concat([auprc_df_all, auprc_df])

# print(auprc_df_all.head())
nhoods_df_all.to_csv(simdir + '/nhoods_obs.csv')
tpr_df_all.to_csv(simdir + '/TPR_res.csv')
auprc_df_all.to_csv(simdir + '/AUPRC_res.csv')
