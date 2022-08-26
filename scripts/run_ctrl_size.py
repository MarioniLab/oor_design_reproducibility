import os
import numpy as np
import scanpy as sc
import anndata
import scvi

import diff2atlas

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("sim_dir",
                    type=str,
                    help="path to simulation directory")
parser.add_argument("n_controls",
                    type=int,
                    help="number of control donors to use")
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


def main(sim_dir, n_controls, random_seed):
    # Load query and control
    adata_query = sc.read_h5ad(sim_dir + '/query.h5ad')
    adata_ctrl = sc.read_h5ad(sim_dir + '/ctrl.h5ad')

    # Make new dir for outputs
    outdir = os.path.join(sim_dir, 'ctrl_size_analysis/')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    sim_id = 'ctrl_size{n_controls}_seed{random_seed}'.format(**locals())
    outdir = outdir + sim_id + '/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Subset control donors
    np.random.seed(random_seed)
    ctrl_donors = adata_ctrl.obs['donor_id'].unique()
    sample_ctrls = np.random.choice(ctrl_donors, n_controls)
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
    vae_fit_ctrl2atlas = diff2atlas.model_wrappers.fit_scVI(
        sim_dir + '/model_reference',
        adata_ctrl_fit,
        outfile=outdir + "/model_fit_ctrl2atlas/")
    vae_fit_query2atlas = scvi.model.SCVI.load(
        sim_dir + "/model_fit_query2atlas/")

    adata_PAC.obsm['X_scVI'] = np.vstack([
        vae_fit_query2atlas.get_latent_representation(),
        vae_fit_ctrl2atlas.get_latent_representation()
    ])

    # Check that all the cells have the same order
    assert (adata_jointPC.obs_names == adata_merge.obs_names).all()
    assert (adata_PC.obs_names == adata_merge.obs_names).all()
    assert (adata_PAC.obs_names == adata_merge.obs_names).all()

    sc.pp.neighbors(adata_jointPC, use_rep='X_scVI', n_neighbors=100)
    sc.pp.neighbors(adata_PC, use_rep='X_scVI', n_neighbors=100)
    sc.pp.neighbors(adata_PAC, use_rep='X_scVI', n_neighbors=100)

    adata_jointPC.write_h5ad(outdir + "jointPC_design.h5ad")
    adata_PC.write_h5ad(outdir + "PC_design.h5ad")
    adata_PAC.write_h5ad(outdir + "PAC_design.h5ad")


# resdir = '/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/'
# sim_dir = resdir  + [x for x in os.listdir(resdir) if x.startswith('qPBMC')][2]
# n_controls = 3
# random_seed = 12345 ## seed for sampling of control cells
main(args.sim_dir, args.n_controls, args.random_seed)
