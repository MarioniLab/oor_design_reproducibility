#!/bin/bash

conda activate oor-benchmark
cd /nfs/team205/ed6/bin/diff2atlas/src/2_simulation_design/

outdir=/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/tmp/

for p in $(cat /lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/PBMC_merged.normal.subsample500cells.clean_celltypes.txt); do
    for d in ACR CR; do
        # echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --embedding_method scVI --outpath ${outdir}" | \
        #     bsub -G teichlab -o logfile-scvi-%J.out -e logfile-scvi-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
        # echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --diff_method meld --outpath ${outdir}" | \
        #     bsub -G teichlab -o logfile-meld-%J.out -e logfile-meld-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
        echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --diff_method cna --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-cna-%J.out -e logfile-cna-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    done
done

d=AR
for p in $(cat /lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/PBMC_merged.normal.subsample500cells.clean_celltypes.txt); do
    echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --diff_method meld --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-meld-%J.out -e logfile-meld-%J.err -q gpu-normal -M100000 -R "select[mem>100000] rusage[mem=100000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
done

simdirs=$(ls -d $outdir/qPBMC_500cells_demo_perturb_cell_type*)
for s in $simdirs; do
    python parse_oor_design_v2.py $s --embedding_method scArches --diff_method meld
    # python parse_oor_design_v2.py $s --embedding_method scArches --diff_method cna
    # echo "
    # python parse_oor_design.py ${s}
    # " | \
    #         bsub -G teichlab -o logfile-scvi-%J.out -e logfile-scvi-%J.err -q normal -M25000 -R "select[mem>25000] rusage[mem=25000]" 
done