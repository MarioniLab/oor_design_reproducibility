import os
import numpy as np
import scanpy as sc
import anndata
import milopy

import diff2atlas

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
args = parser.parse_args()


def clean_pop_name(string):
    return(''.join(e if e.isalnum() else '_' for e in string))


def filter_genes_scvi(adata_PC_train):
    # Filter genes not expressed anywhere
    sc.pp.filter_genes(adata_PC_train, min_cells=1)

    # Select HVGs
    if 'log1p' not in adata_PC_train.uns.keys():
        sc.pp.normalize_per_cell(adata_PC_train)
        sc.pp.log1p(adata_PC_train)

    sc.pp.highly_variable_genes(
        adata_PC_train,
        n_top_genes=5000,
        subset=True
    )


def run_milo(adata_design, query_group, reference_group,
             sample_col='sample_id',
             annotation_col='cell_type',
             design='~ is_query'
             ):
    milopy.core.make_nhoods(adata_design, prop=0.1)
    milopy.core.count_nhoods(adata_design, sample_col=sample_col)
    milopy.utils.annotate_nhoods(
        adata_design[adata_design.obs['dataset_group'] == reference_group], annotation_col)
    adata_design.obs['is_query'] = adata_design.obs['dataset_group'] == query_group
    milopy.core.DA_nhoods(adata_design, design=design)


def main(sim_dir, n_controls, n_querys, random_seed):
    # Load query and control
    adata_query = sc.read_h5ad(sim_dir + '/query.h5ad')
    adata_ctrl = sc.read_h5ad(sim_dir + '/ctrl.h5ad')

    # Make new dir for outputs
    outdir = os.path.join(sim_dir, f'ctrl_size_analysis_nquery{n_querys}/')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    sim_id = 'ctrl_size{n_controls}_seed{random_seed}'.format(**locals())
    outdir = outdir + sim_id + '/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Subset query donors with most cells in query_specific population (to always have less query than ctrl donors)
    # np.random.seed(random_seed)
    # query_donors = adata_query.obs['donor_id'].unique()
    # sample_querys = np.random.choice(query_donors, n_querys, replace=False)
    perturb_pop = sim_dir.split(
        '/')[-1].split('perturb_cell_type')[-1].split('_queryBatch')[0]
    sample_querys = adata_query.obs[adata_query.obs['cell_type'] == perturb_pop].value_counts(
        'donor_id')[0:n_querys].index.values
    adata_query = adata_query[adata_query.obs['donor_id'].isin(sample_querys)]

    # Subset control donors
    np.random.seed(random_seed)
    ctrl_donors = adata_ctrl.obs['donor_id'].unique()
    sample_ctrls = np.random.choice(ctrl_donors, n_controls, replace=False)
    adata_ctrl = adata_ctrl[adata_ctrl.obs['donor_id'].isin(sample_ctrls)]

    adata_merge = anndata.concat([adata_query, adata_ctrl])
    adata_merge.layers['counts'] = adata_merge.X.copy()
    adata_query.layers['counts'] = adata_query.X.copy()
    adata_ctrl.layers['counts'] = adata_ctrl.X.copy()

    ## --- PC design joint --- ##
    adata_jointPC = adata_merge.copy()
    adata_jointPC_train = adata_merge.copy()
    filter_genes_scvi(adata_jointPC_train)

    # Train model
    vae_ref = diff2atlas.model_wrappers.train_scVI(
        adata_jointPC_train, outdir + f"/model_ctrlquery_joint/", batch_key='sample_id')

    adata_jointPC.obsm['X_scVI'] = vae_ref.get_latent_representation()

    ## --- PC design with mapping --- ##
    adata_PC = adata_merge.copy()
    adata_ctrl_train = adata_ctrl.copy()
    filter_genes_scvi(adata_ctrl_train)

    # Train model
    vae_ctrl = diff2atlas.model_wrappers.train_scVI(
        adata_ctrl_train, outdir + f"/model_ctrl/", batch_key='sample_id')

    # Fit query
    adata_query_fit = adata_query.copy()
    vae_fit_query2ctrl = diff2atlas.model_wrappers.fit_scVI(
        vae_ctrl,
        adata_query_fit,
        outfile=outdir + "/model_fit_query2ctrl/")

    adata_PC.obsm['X_scVI'] = np.vstack([
        vae_fit_query2ctrl.get_latent_representation(),
        vae_ctrl.get_latent_representation()
    ])

    ## -- PAC design -- ##
    adata_PAC = adata_merge.copy()

    # map ctrl to atlas and load query model
    adata_ctrl_fit = adata_ctrl.copy()
    adata_query_fit = adata_query.copy()
    vae_fit_ctrl2atlas = diff2atlas.model_wrappers.fit_scVI(
        sim_dir + '/model_reference',
        adata_ctrl_fit,
        outfile=outdir + "/model_fit_ctrl2atlas/")
    vae_fit_query2atlas = diff2atlas.model_wrappers.fit_scVI(
        sim_dir + '/model_reference',
        adata_query_fit,
        outfile=outdir + "/model_fit_query2atlas/")
    # vae_fit_query2atlas = scvi.model.SCVI.load(
    #     sim_dir + "/model_fit_query2atlas/")

    adata_PAC.obsm['X_scVI'] = np.vstack([
        vae_fit_query2atlas.get_latent_representation(),
        vae_fit_ctrl2atlas.get_latent_representation()
    ])

    # Check that all the cells have the same order
    assert (adata_jointPC.obs_names == adata_merge.obs_names).all()
    assert (adata_PC.obs_names == adata_merge.obs_names).all()
    assert (adata_PAC.obs_names == adata_merge.obs_names).all()

    sc.pp.neighbors(adata_jointPC, use_rep='X_scVI',
                    n_neighbors=(n_controls+n_querys)*5)
    sc.pp.neighbors(adata_PC, use_rep='X_scVI',
                    n_neighbors=(n_controls+n_querys)*5)
    sc.pp.neighbors(adata_PAC, use_rep='X_scVI',
                    n_neighbors=(n_controls+n_querys)*5)

    # Run milo analysis
    adata_PAC.obs['dataset_group'] = np.where(
        adata_PAC.obs['is_test'] == 1, 'query', 'ctrl')
    adata_PC.obs['dataset_group'] = np.where(
        adata_PC.obs['is_test'] == 1, 'query', 'ctrl')
    adata_jointPC.obs['dataset_group'] = np.where(
        adata_jointPC.obs['is_test'] == 1, 'query', 'ctrl')

    run_milo(adata_PAC, query_group='query', reference_group='ctrl')
    run_milo(adata_jointPC, query_group='query', reference_group='ctrl')
    run_milo(adata_PC, query_group='query', reference_group='ctrl')

    milopy.utils.write_milo_adata(adata_PAC, outdir + '/PAC_design.h5ad')
    milopy.utils.write_milo_adata(adata_PC, outdir + '/PC_design.h5ad')
    milopy.utils.write_milo_adata(
        adata_jointPC, outdir + '/jointPC_design.h5ad')


# resdir = '/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/'
# sim_dir = resdir  + [x for x in os.listdir(resdir) if x.startswith('qPBMC')][2]
# n_controls = 3
# random_seed = 12345 ## seed for sampling of control cells
main(args.sim_dir, args.n_controls, args.n_querys, args.random_seed)
