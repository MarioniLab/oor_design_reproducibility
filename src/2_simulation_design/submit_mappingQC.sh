#!/bin/bash

conda activate oor-benchmark
cd /nfs/team205/ed6/bin/diff2atlas/scripts

outdir=/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/tmp/
h5ad_files=$(ls $outdir/qPBMC_500cells_demo_perturb_cell_type*/ar_design.h5ad)

for s in $h5ad_files; do
    python run_mappingQC.py $s; 
    done

for s in $h5ad_files; do
    echo "python run_mappingQC_reconstruction.py $s" | \
        bsub -G teichlab -o logfile-mappingQC-%J.out -e logfile-mappingQC-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
done

