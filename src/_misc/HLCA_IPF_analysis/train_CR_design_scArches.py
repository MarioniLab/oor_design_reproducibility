import pandas as pd
import scanpy as sc
import scvi
from anndata import AnnData

data_dir = '/lustre/scratch117/cellgen/team205/ed6/HLCA/'

## Read CR design object
adata = sc.read_h5ad(data_dir + "Kaminski_2020_oor_design.CR.h5ad")

## Filter genes
adata_train = adata[:,adata.var['mapping_gene']].copy()

## Set up model
outdir = data_dir + "/model_ctrl/"

# split in query and control
adata_query = adata_train[adata_train.obs['disease'] != 'Control'].copy()
adata_train = adata_train[adata_train.obs['disease'] == 'Control'].copy()

## Params from https://github.com/LungCellAtlas/HLCA_reproducibility/blob/main/notebooks/1_building_and_annotating_the_atlas_core/05_scANVI_integration.ipynb    
condition_key = 'dataset'
cell_type_key = 'scanvi_labels'
unlabeled_category = "unlabeled"

vae_epochs = 500
scanvi_epochs = 200

early_stopping_kwargs = {
    "early_stopping_patience": 10
}

adata_train.obs['dataset'] = 'Kaminski_2020'
scvi.model.SCANVI.setup_anndata(
    adata_train, batch_key=condition_key, labels_key=cell_type_key, unlabeled_category=unlabeled_category)

vae = scvi.model.SCANVI(
    adata_train,
    n_layers=2,
    n_latent = 30, # to allow for capturing more heterogeneity
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
    gene_likelihood="nb", # because we have UMI data
)

vae.train(
    max_epochs = vae_epochs, 
    **early_stopping_kwargs
)

vae.save(outdir, save_anndata=True, overwrite=True)
adata_train.obsm["X_scVI"] = vae.get_latent_representation()

# Map query dataset 
vae_q = scvi.model.SCANVI.load_query_data(adata_query, vae, inplace_subset_query_vars=True)
vae_q.train(plan_kwargs={"weight_decay": 0.0})

# Get latent embeddings
X_scVI_ref = pd.DataFrame(vae.get_latent_representation(), index=vae.adata.obs_names)
X_scVI_q = pd.DataFrame(vae_q.get_latent_representation(), index=vae_q.adata.obs_names)
X_scVI = pd.concat([X_scVI_q, X_scVI_ref], axis=0)
adata.obsm["X_scVI"] = X_scVI.loc[adata.obs_names].values

adata.write_h5ad(data_dir + "Kaminski_2020_oor_design.CR.scArches.h5ad")