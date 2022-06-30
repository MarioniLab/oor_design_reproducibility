
import scanpy as sc
import milopy

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("sim_dir",
                    type=str,
                    help="path to simulation directory")
args = parser.parse_args()


def run_milo(pac_design_adata,
             query_group,
             reference_group,
             sample_col='sample_id',
             annotation_col='cell_type'
             ):
    milopy.core.make_nhoods(pac_design_adata, prop=0.1)
    milopy.core.count_nhoods(pac_design_adata, sample_col=sample_col)
    milopy.utils.annotate_nhoods(
        pac_design_adata[pac_design_adata.obs['dataset_group'] == reference_group], annotation_col)
    pac_design_adata.obs['is_query'] = pac_design_adata.obs['dataset_group'] == query_group
    milopy.core.DA_nhoods(pac_design_adata, design='~ is_query')
    # milopy.utils.build_nhood_graph(pac_design_adata)


def main(sim_dir):
    pac_design_adata = sc.read_h5ad(sim_dir + '/pac_design.h5ad')
    run_milo(pac_design_adata, query_group='query', reference_group='ctrl')
    milopy.utils.write_milo_adata(
        pac_design_adata, sim_dir + '/pac_design.h5ad')

    pa_design_adata = sc.read_h5ad(sim_dir + '/pa_design.h5ad')
    run_milo(pa_design_adata, query_group='query', reference_group='atlas')
    milopy.utils.write_milo_adata(pa_design_adata, sim_dir + '/pa_design.h5ad')

    pc_design_adata = sc.read_h5ad(sim_dir + '/pc_design.h5ad')
    run_milo(pc_design_adata, query_group='query', reference_group='ctrl')
    milopy.utils.write_milo_adata(pc_design_adata, sim_dir + '/pc_design.h5ad')

main(args.sim_dir)