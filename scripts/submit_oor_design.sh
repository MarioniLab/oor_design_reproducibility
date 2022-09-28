#!/bin/bash

conda activate oor-benchmark
cd /nfs/team205/ed6/bin/diff2atlas/scripts

outdir=/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/tmp/

for p in $(cat /lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/PBMC_merged.normal.subsample500cells.clean_celltypes.txt); do
    for d in ACR AR CR; do
        echo "python run_oor_design.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-scvi-%J.out -e logfile-scvi-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    done
done

## Test AR design requirements
p=central_memory_CD4_positive_alpha_beta_T_cell
d=AR
echo "python run_oor_design.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-scvi-%J.out -e logfile-scvi-%J.err -q gpu-normal -M200000 -R "select[mem>200000] rusage[mem=200000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
