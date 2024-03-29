# Disease-state identification with healthy single-cell references 
[![DOI](https://zenodo.org/badge/503372743.svg)](https://zenodo.org/badge/latestdoi/503372743)

This repository contains notebooks and scripts to reproduce analyses benchmarking the use of control and atlas datasets as references for identification of disease-associated cell states (see [manuscript](https://doi.org/10.1101/2022.11.10.515939)).

The workflow for disease-state identification and evaluation of out-of-reference detection is available as a [python package](https://github.com/emdann/oor_benchmark). 

## Repository structure

- `diff2atlas` - utility module 
- `metadata` - metadata tables used for analysis
- `src` - analysis notebooks and scripts 
  - `1_PBMC_data_preprocessing/` - preprocessing and harmonization of PBMC dataset
  - `2_simulation_design/` - out-of-reference detection benchmark on simulations
  - `3_simulation_ctrl_atlas_size` - out-of-reference detection robustness to atlas and control dataset size
  - `3b_crosstissue_atlas` - out-of-reference detection robustness with tissue-matched or cross-tissue atlas
  - `4_COVID_design` - reference design comparison on COVID-19 dataset 
  - `5_IPF_HLCA_design` - reference design comparison on IPF lung dataset

## Data

Processed datasets and scVI models used in this analysis are available via [figshare](https://doi.org/10.6084/m9.figshare.21456645). For references of the original datasets collected see [study metadata](https://github.com/MarioniLab/oor_design_reproducibility/blob/master/metadata/suppl_table_studies.csv). 

For simulation analysis
- `PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad` - harmonized object of healthy PBMC profiles from 13 studies, used for OOR identification benchmark with simulations
- `model_PBMC_merged.normal.subsample500cells.zip` - scVI model trained on healthy PBMC profiles (used for joint annotation) (trained with scvi-tools v0.16.2, see [notebooks](https://github.com/MarioniLab/oor_design_reproducibility/blob/master/src/1_PBMC_data_preprocessing/20220601_PBMC_scVI.ipynb) for training parameters)
- Results from simulation analysis are shared in .csv files (`OOR_simulations_*.csv`)
    - `*.nhood_results_all.csv` - neighbourhood level Milo results (with fraction of OOR state)
    - `*.TPRFPRFDR_results_all.csv` - TPR/FDR/FPR for each simulation
    - `*.AUPRC_results_all.csv` - AUPRC for each simulation

For COVID-19 analysis
- `PBMC_COVID.subsample500cells.atlas.h5ad` - atlas dataset (PBMCs from healthy individuals from 12 studies)
- `PBMC_COVID.subsample500cells.covid.h5ad`- disease dataset (PBMCs from COVID-19 patients from [Stephenson et al. 2021](https://www.nature.com/articles/s41591-021-01329-2))
- `PBMC_COVID.subsample500cells.ctrl.h5ad` - control dataset (PBMCs from healthy individuals from [Stephenson et al. 2021](https://www.nature.com/articles/s41591-021-01329-2))
- `PBMC_COVID.subsample500cells.design.query_PC_refA.post_milo.h5ad` - ACR design processed object with Milo results
- `PBMC_COVID.subsample500cells.design.query_PC_refA.post_milo.nhood_adata.h5ad` - ACR design processed object with Milo results (nhood AnnData)
- `PBMC_COVID.subsample500cells.design.query_P_refC.post_milo.h5ad` - CR design processed object with Milo results
- `PBMC_COVID.subsample500cells.design.query_P_refC.post_milo.nhood_adata.h5ad` - CR design processed object with Milo results (nhood AnnData)
- `model_COVID19_reference_atlas_scvi0.16.2.zip` - scVI model trained on atlas dataset (used for ACR design) (trained with scvi-tools v0.16.2, see [script](https://github.com/MarioniLab/oor_design_reproducibility/blob/master/src/4_COVID_design/COVID_train_references.py) for training parameters)

For IPF analysis
- `IPF_HLCA.ACR_design.post_milo.h5ad` - ACR design processed object with Milo results. Includes annotation of aberrant basal-like states (`adata.obs['basal_like_annotation']`)
- `IPF_HLCA.ACR_design.post_milo.nhood_adata.h5ad` - ACR design processed object with Milo results (nhood AnnData)
- `IPF_HLCA.CR_design.post_milo.h5ad` - CR design processed object with Milo results.
- `IPF_HLCA.CR_design.post_milo.nhood_adata.h5ad` - CR design processed object with Milo results (nhood AnnData)
- `IPF_HLCA.AR_design.post_milo.h5ad` - AR design processed object with Milo results.
- `IPF_HLCA.AR_design.post_milo.nhood_adata.h5ad` - AR design processed object with Milo results (nhood AnnData)

For cross-tissue atlas analysis
- `model_TabulaSapiens_scvi0.20.0.zip` - scVI model trained on Tabula Sapiens dataset (trained with scvi-tools v0.20.0, see [script](https://github.com/MarioniLab/oor_design_reproducibility/blob/revision-1.0/src/3b_crosstissue_atlas/train_atlas.py) for training parameters)

## Citation

> Dann E., Teichmann S.A. and Marioni J.C. Precise identification of cell states altered in disease with healthy single-cell references. biorXiv https://doi.org/10.1101/2022.11.10.515939

## Contact

For any questions, please post an [issue](https://github.com/MarioniLab/oor_design_reproducibility/issues?q=is%3Aissue+is%3Aopen+sort%3Aupdated-desc).


