# diff2atlas

## Motivation

Large single-cell datasets capturing cells from several healthy individuals (atlases) promise to serve as reference datasets for discovery of cellular phenotypes associated with disease, genetic or environmental perturbations. Several computational models have been proposed to contextualize new data with existing atlases, by learning common latent dimensions or cell type labels.

However, it remains unclear how such models should be used for accurate identification of disease-specific changes: contrasting cells from disease conditions to an all-encompassing atlas (perturbation-atlas design) fails to disentangle whether deviation from the reference is driven by the condition of interest or technical and biological variation within the atlas. Conversely, comparing only to matching controls (perturbation-control design) fails to distinguish where disease effects differ from natural variation. 

**When does this matter:** in the presence of biological variation in the atlas - aggregated gene expression profiles or cell type/state proportions are not representative of true gene expression profiles and proportions found in tissues. This could lead to false positives (differences driven by context rather than condition of interest), or false negatives (differences missed because specific to the context)

