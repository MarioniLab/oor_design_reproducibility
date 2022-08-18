
import scanpy as sc
import milopy

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("sim_dir",
                    type=str,
                    help="path to simulation directory")
parser.add_argument("--prop",
                    type=str, default='fixed',
                    help="either fixed (10%) or target")
args = parser.parse_args()


def run_milo(pac_design_adata,
             query_group,
             reference_group,
             sample_col='sample_id',
             annotation_col='cell_type',
             target_n_nhoods = None
             ):
    if target_n_nhoods:
        ## Define prop based on no of cells to have ~ no of nhoods across designs
        prop = target_n_nhoods/pac_design_adata.n_obs
    else:
        prop=0.1
    milopy.core.make_nhoods(pac_design_adata, prop=prop)
    milopy.core.count_nhoods(pac_design_adata, sample_col=sample_col)
    milopy.utils.annotate_nhoods(
        pac_design_adata[pac_design_adata.obs['dataset_group'] == reference_group], annotation_col)
    pac_design_adata.obs['is_query'] = pac_design_adata.obs['dataset_group'] == query_group
    milopy.core.DA_nhoods(pac_design_adata, design='~ is_query')

def main(sim_dir, prop):
    if prop == 'target':
        target_n_nhoods = 5000
    else:
        target_n_nhoods = None
    pac_design_adata = sc.read_h5ad(sim_dir + '/pac_design.h5ad')
    run_milo(pac_design_adata, query_group='query', reference_group='ctrl', target_n_nhoods=target_n_nhoods)
    milopy.utils.write_milo_adata(
        pac_design_adata, sim_dir + '/pac_design.h5ad')

    pa_design_adata = sc.read_h5ad(sim_dir + '/pa_design.h5ad')
    run_milo(pa_design_adata, query_group='query', reference_group='atlas', target_n_nhoods=target_n_nhoods)
    milopy.utils.write_milo_adata(pa_design_adata, sim_dir + '/pa_design.h5ad')

    pc_design_adata = sc.read_h5ad(sim_dir + '/pc_design.h5ad')
    run_milo(pc_design_adata, query_group='query', reference_group='ctrl', target_n_nhoods=target_n_nhoods)
    milopy.utils.write_milo_adata(pc_design_adata, sim_dir + '/pc_design.h5ad')

main(args.sim_dir, args.prop)