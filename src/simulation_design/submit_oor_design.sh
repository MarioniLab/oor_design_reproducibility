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

simdirs=$(ls -d $outdir/qPBMC_500cells_demo_perturb_cell_type*)
for s in $simdirs; do
    # echo "
    python parse_oor_design.py ${s}
    # " | \
    #         bsub -G teichlab -o logfile-scvi-%J.out -e logfile-scvi-%J.err -q normal -M25000 -R "select[mem>25000] rusage[mem=25000]" 
done

