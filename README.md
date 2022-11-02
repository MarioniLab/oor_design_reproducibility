# Disease-state identification with healthy single-cell references  

This repository contains notebooks and scripts to reproduce analyses benchmarking the use of control and atlas datasets as references for identification of disease-associated cell states (see [manuscript]()).

The workflow for disease-state identification and evaluation of out-of-reference detection is available as a [python package](https://github.com/emdann/oor_benchmark). 

## Repository structure

- `diff2atlas` - utility module 
- `metadata` - study and sample level metadata
- `src` - analysis notebooks and scripts 
  - `1_PBMC_data_preprocessing/` - preprocessing and harmonization of PBMC dataset
  - `2_simulation_design/` - out-of-reference detection benchmark on simulations
  - `3_simulation_ctrl_atlas_size` - out-of-reference detection robustness to atlas and control dataset size
  - `4_COVID_design` - reference design comparison on COVID-19 dataset 

## Data

Processed datasets and scVI models used in this analysis are available via [figshare](). For references of the original datasets collected see [study metadata](). 

For simulation analysis
- `PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad` - harmonized object of healthy PBMC profiles from 13 studies, used for OOR identification benchmark with simulations
- `model_PBMC_merged.normal.subsample500cells.zip` - scVI model trained on healthy PBMC profiles (used for joint annotation) 

For COVID-19 analysis
- `PBMC_COVID.subsample500cells.atlas.h5ad` - atlas dataset (PBMCs from healthy individuals from 12 studies)
- `PBMC_COVID.subsample500cells.covid.h5ad`- disease dataset (PBMCs from COVID-19 patients from [Stephenson et al. 2021](https://www.nature.com/articles/s41591-021-01329-2))
- `PBMC_COVID.subsample500cells.ctrl.h5ad` - control dataset (PBMCs from healthy individuals from [Stephenson et al. 2021](https://www.nature.com/articles/s41591-021-01329-2))
- `PBMC_COVID.subsample500cells.design.query_PC_refA.post_milo.h5ad` - ACR design processed object with Milo results (load with [`milopy.utils.read_milo_adata`](https://milopy.readthedocs.io/en/latest/autoapi/milopy/utils/index.html#milopy.utils.read_milo_adata)).
- `PBMC_COVID.subsample500cells.design.query_PC_refA.post_milo.nhood_adata.h5ad` - ACR design processed object with Milo results (nhood AnnData)
- `PBMC_COVID.subsample500cells.design.query_P_refC.post_milo.h5ad` - CR design processed object with Milo results (load with [`milopy.utils.read_milo_adata`](https://milopy.readthedocs.io/en/latest/autoapi/milopy/utils/index.html#milopy.utils.read_milo_adata)).
- `PBMC_COVID.subsample500cells.design.query_P_refC.post_milo.nhood_adata.h5ad` - CR design processed object with Milo results (nhood AnnData)
- `model_COVID19_reference_atlas_scvi0.16.2` - scVI model trained on atlas dataset (used for ACR design)

## Citation

> Coming soon

## Contact

For any questions, please post an [issue](https://github.com/emdann/diff2atlas/issues?q=is%3Aissue+is%3Aopen+sort%3Aupdated-desc).


