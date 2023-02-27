# Disease-state identification with healthy single-cell references  

This repository contains notebooks and scripts to reproduce analyses benchmarking the use of control and atlas datasets as references for identification of disease-associated cell states (see [manuscript](https://doi.org/10.1101/2022.11.10.515939)).

The workflow for disease-state identification and evaluation of out-of-reference detection is available as a [python package](https://github.com/emdann/oor_benchmark). 

## Repository structure

- `diff2atlas` - utility module 
- `metadata` - metadata tables used for analysis
  - `PBMC_study_metadata.csv` - metadata for PBMC studies used in simulation and COVID-19 analysis
  - `PBMC_sample_metadata.csv` - sample level metadata for PBMC studies used in simulation and COVID-19 analysis
  - `efotraits_EFO_0004314-studies-2023-01-24.csv ` - table of GWAS study IDs and metadata for lung function (FEV, EFO_0004314) (see [analysis of aberrant basal-like cells in IPF]()) 
  -  `opentargets_drugs.EFO_0003818.tsv` - table of validated target genes for drugs approved or in trial for lung diseasee GWAS loci and L2G predicted associated genes for lung function (FEV, EFO_0004314) (see [analysis of aberrant basal-like cells in IPF]())
  - `opentargets_drugs.EFO_0003818.tsv` - table of validated target genes for drugs approved or in trial for lung disease (EFO_OOO3818) (see [analysis of aberrant basal-like cells in IPF]())
- `src` - analysis notebooks and scripts 
  - `1_PBMC_data_preprocessing/` - preprocessing and harmonization of PBMC dataset
  - `2_simulation_design/` - out-of-reference detection benchmark on simulations
  - `3_simulation_ctrl_atlas_size` - out-of-reference detection robustness to atlas and control dataset size
  - `3b_crosstissue_atlas` - out-of-reference detection robustness with tissue-matched or cross-tissue atlas
  - `4_COVID_design` - reference design comparison on COVID-19 dataset 
  - `5_IPF_HLCA_design` - reference design comparison on IPF lung dataset

## Data

Processed datasets and scVI models used in this analysis are available via [figshare](https://doi.org/10.6084/m9.figshare.21456645.v1). For references of the original datasets collected see [study metadata](https://github.com/MarioniLab/oor_design_reproducibility/blob/master/metadata/PBMC_study_metadata.csv). 

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

> Dann E., Teichmann S.A. and Marioni J.C. Precise identification of cell states altered in disease with healthy single-cell references. biorXiv https://doi.org/10.1101/2022.11.10.515939

## Contact

For any questions, please post an [issue](https://github.com/MarioniLab/oor_design_reproducibility/issues?q=is%3Aissue+is%3Aopen+sort%3Aupdated-desc).


