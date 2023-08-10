#!/bin/bash

conda activate scvi-env
cd /nfs/team205/ed6/bin/diff2atlas/scripts
outdir=/lustre/scratch126/cellgen/team205/ed6/PBMC_COVID/

## Train reference models
# A = Atlas (for ACR design - scArches)
# C = Control (for CR design - scArches)
# PC = Perturbation+Control (for CR design - scVI)
# PAC = Perturbation+Atlas+Control (for ACR design - scVI)
for d in 'A' 'C' 'PC' 'PAC'; do     
    echo "python oor_design_reproducibility/src/4_COVID_design/COVID_train_references.py ${d} --datadir ${outdir}" | \
    bsub -G teichlab -o logfile-scvi-%J.out -e logfile-scvi-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1"
done

## Map query data

# ACR design - scArches
echo "python oor_design_reproducibility/src/4_COVID_design/COVID_map_query.py model_reference_A PC --datadir ${outdir}" | \
bsub -G teichlab -o logfile-scarches-%J.out -e logfile-scarches-%J.err -q gpu-normal -M25000 -R "select[mem>25000] rusage[mem=25000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1"

# CR design - scArches
echo "python COVID_map_query.py model_reference_C P --datadir ${outdir}" | \
    bsub -G teichlab -o logfile-scarches-%J.out -e logfile-scarches-%J.err -q gpu-normal -M25000 -R "select[mem>25000] rusage[mem=25000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1"
