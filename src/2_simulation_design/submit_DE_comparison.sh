#!/bin/bash

conda activate oor-benchmark
cd /nfs/team205/ed6/bin/diff2atlas/src/2_simulation_design/
outdir=/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/tmp/

dirs=$(ls -d $outdir/qPBMC_500cells_demo_perturb_cell_type*seed2022/)

for d in $dirs; do \
    echo "python pseudobulk_for_DE.py ${d} ACR scArches" | \
        bsub -G teichlab -o logfile-pbulkde-%J.out -e logfile-pbulkde-%J.err -M50000 -R "select[mem>50000] rusage[mem=50000]" 
    echo "python pseudobulk_for_DE.py ${d} CR scArches" | \
        bsub -G teichlab -o logfile-pbulkde-%J.out -e logfile-pbulkde-%J.err -M50000 -R "select[mem>50000] rusage[mem=50000]" 
    echo "python pseudobulk_for_DE.py ${d} CR scVI" | \
        bsub -G teichlab -o logfile-pbulkde-%J.out -e logfile-pbulkde-%J.err -M50000 -R "select[mem>50000] rusage[mem=50000]" 
    echo "python pseudobulk_for_DE.py ${d} AR scArches" | \
        bsub -G teichlab -o logfile-pbulkde-%J.out -e logfile-pbulkde-%J.err -M50000 -R "select[mem>50000] rusage[mem=50000]" 
done 

for d in $dirs; do \
    python run_DE_comparison.py ${d} AR scArches
    # " | \
    #     bsub -G teichlab -o logfile-de-%J.out -e logfile-de-%J.err -M50000 -R "select[mem>50000] rusage[mem=50000]" 
    # python run_DE_comparison.py ${d} CR scArches
        # bsub -G teichlab -o logfile-de-%J.out -e logfile-de-%J.err -M50000 -R "select[mem>50000] rusage[mem=50000]" 
    # python run_DE_comparison.py ${d} CR scVI 
        # bsub -G teichlab -o logfile-de-%J.out -e logfile-de-%J.err -M50000 -R "select[mem>50000] rusage[mem=50000]" 
    # echo "python run_DE_comparison.py ${d} AR scArches" | \
    #     bsub -G teichlab -o logfile-de-%J.out -e logfile-de-%J.err -M50000 -R "select[mem>50000] rusage[mem=50000]" 
done 