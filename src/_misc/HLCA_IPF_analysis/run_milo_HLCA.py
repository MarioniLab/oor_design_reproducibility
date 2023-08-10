from milopy.utils import write_milo_adata
import scanpy as sc
import milopy
import os

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("design",
                    type=str,
                    help="ID of reference design (CR, ACR, AR)")
parser.add_argument("--write_tmp", action='store_true',
                    help="Save intermediate results in tmp directory")
args = parser.parse_args()

design = args.design
write_tmp = args.write_tmp
    
# Define constants
sample_col = 'sample'
annotation_col = 'Celltype_HLCA'
milo_design = '~disease'
data_dir = '/lustre/scratch117/cellgen/team205/ed6/HLCA/'


tmp_dir = data_dir + "/tmp/"
if not os.path.exists(tmp_dir):
    os.mkdir(tmp_dir)

print('Reading data\n')
try:
    adata = sc.read_h5ad(f'{tmp_dir}/Kaminski_2020_oor_design.{design}.h5ad')
except FileNotFoundError:
    adata = sc.read_h5ad(f'{data_dir}/Kaminski_2020_oor_design.{design}.h5ad')

    if design == 'AR':
        ## Remove the control cells 
        remove_obs = adata.obs_names[(adata.obs['study'] == 'Kaminski_2020') & (adata.obs['disease'] == 'Control')]
        adata = adata[~adata.obs_names.isin(remove_obs)].copy()

    print('Running KNN graph\n')
    # Make KNN graph for Milo neigbourhoods
    n_samples = adata.obs[sample_col].unique().shape[0]
    #  Set max to 200 or memory explodes for large datasets
    k = min([n_samples * 5, 200])
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=k)

    if write_tmp:
        adata.write_h5ad(f'{tmp_dir}/Kaminski_2020_oor_design.{design}.h5ad')

# Run Milo
print('Running Milo\n')
adata.obs['disease'] = adata.obs['disease'].str.replace('-| ', '_') # fix names or R complains
milopy.core.make_nhoods(adata, prop=0.01)
milopy.core.count_nhoods(adata, sample_col=sample_col)
# milopy.utils.annotate_nhoods(adata[adata.obs["disease"] == 'Control'], annotation_col)
milopy.core.DA_nhoods(adata, design=milo_design, model_contrasts='diseaseIPF-diseaseControl')

# Write output
print('Saving output\n')
write_milo_adata(adata, f'{data_dir}/Kaminski_2020_oor_design.{design}.h5ad')