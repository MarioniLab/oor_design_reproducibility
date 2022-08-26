#!/bin/bash

conda activate scvi-env
cd /nfs/team205/ed6/bin/diff2atlas/scripts

outdir=/lustre/scratch117/cellgen/team205/ed6/PBMC_COVID/

## Train reference models
for d in 'A' 'C' 'AC' 'PA' 'PC' 'PAC'; do     
    echo "python COVID_train_reference.py ${d} --datadir ${outdir}" | \
    bsub -G teichlab -o logfile-scvi-%J.out -e logfile-scvi-%J.err -q gpu-normal -M25000 -R "select[mem>25000] rusage[mem=25000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1"
done

## Map query
for q in 'P' 'PC'; do
    echo "python COVID_map_query.py model_reference_A ${q} --datadir ${outdir}" | \
    bsub -G teichlab -o logfile-scarches-%J.out -e logfile-scarches-%J.err -q gpu-normal -M25000 -R "select[mem>25000] rusage[mem=25000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1"
done

echo "python COVID_map_query.py model_reference_AC P --datadir ${outdir}" | \
    bsub -G teichlab -o logfile-scarches-%J.out -e logfile-scarches-%J.err -q gpu-normal -M25000 -R "select[mem>25000] rusage[mem=25000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1"


echo "python COVID_map_query.py model_reference_C P --datadir ${outdir}" | \
    bsub -G teichlab -o logfile-scarches-%J.out -e logfile-scarches-%J.err -q gpu-normal -M25000 -R "select[mem>25000] rusage[mem=25000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1"




# for p in $(cat /lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/PBMC_merged.normal.subsample500cells.clean_celltypes.txt); do
#     echo "python run_design_comparison.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} --outpath ${outdir}" | \
#         bsub -G teichlab -o logfile-scvi-%J.out -e logfile-scvi-%J.err -q gpu-normal -M25000 -R "select[mem>25000] rusage[mem=25000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
# done

# for d in $(ls -d $outdir*/); do
#     echo "python run_design_merge.py ${d}/" | \
#         bsub -G teichlab -o logfile-merge-%J.out -e logfile-merge-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
# done

# for d in $(ls -d $outdir*/); do
#     echo "python run_milo.py ${d}/" | \
#         bsub -G teichlab -o logfile-milo-%J.out -e logfile-milo-%J.err -q gpu-normal -M30000 -R "select[mem>30000] rusage[mem=30000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
# done

