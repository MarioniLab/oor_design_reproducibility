
import scanpy as sc
import pandas as pd
import numpy as np
import milopy.utils
import sklearn.metrics
# pip install leidenalg

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("simdir",
                    type=str,
                    help="path to simulation directory")
parser.add_argument("--filter_dataset", default=None,
                    help="Which type of dataset in adata.obs['dataset_group'] to filter")
args = parser.parse_args()

simdir = args.simdir
filter_dataset = args.filter_dataset

if filter_dataset is not None:
    b = False
else:
    b = True

acr_adata = milopy.utils.read_milo_adata(
    simdir + '/acr_design.h5ad', backed=b)
ar_adata = milopy.utils.read_milo_adata(
    simdir + '/ar_design.h5ad', backed=b)
cr_adata = milopy.utils.read_milo_adata(
    simdir + '/cr_design.h5ad', backed=b)
ct = acr_adata.obs[acr_adata.obs['OOR_state'] == 1]['cell_type'].unique()[0]
design_dict = {'CR': cr_adata, 'ACR': acr_adata, 'AR': ar_adata}
if filter_dataset is not None:
    for d in design_dict:
        if any(design_dict[d].obs['dataset_group'] == filter_dataset):
            design_dict[d] = design_dict[d][design_dict[d].obs['dataset_group']
                                            == filter_dataset].copy()
            sc.pp.neighbors(design_dict[d], n_neighbors=50, use_rep='X_scVI')
        else:
            design_dict.pop(d)

int_quality_df = pd.DataFrame()
for d, adata in design_dict.items():
    for res in np.arange(0.5, 2.25, 0.25):
        print(f'design: {d} / resolution: {res}')
        sc.tl.leiden(adata, resolution=res, key_added=f'leiden_res{res}')
        nmi_ct = sklearn.metrics.normalized_mutual_info_score(
            adata.obs[f'leiden_res{res}'], adata.obs['cell_type'])
        nmi_batch = sklearn.metrics.normalized_mutual_info_score(
            adata.obs[f'leiden_res{res}'], adata.obs['sample_id'])
        df = pd.DataFrame([d, res, nmi_ct, nmi_batch, ct]).T
        df.columns = ['design', 'resolution',
                      'NMI_celltype', 'NMI_batch', "OOR_state_name"]
        int_quality_df = pd.concat([int_quality_df, df])

if filter_dataset is not None:
    int_quality_df.to_csv(
        simdir + f'/integration_quality_filter{filter_dataset}.csv')
else:
    int_quality_df.to_csv(simdir + '/integration_quality_NMI.csv')
