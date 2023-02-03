#!/bin/bash

conda activate oor-benchmark
cd /nfs/team205/ed6/bin/diff2atlas/src/2_simulation_design/

outdir=/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/tmp/

for p in $(cat /lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/PBMC_merged.normal.subsample500cells.clean_celltypes.txt); do
    for d in ACR CR; do
        echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --diff_method meld --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-meld-%J.out -e logfile-meld-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
        echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --diff_method cna --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-cna-%J.out -e logfile-cna-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
        echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --embedding_method scArches --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-milo-%J.out -e logfile-milo-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    done
done

for p in $(cat /lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/PBMC_merged.normal.subsample500cells.clean_celltypes.txt); do
    for d in ACR CR; do
        echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --embedding_method scVI --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-scvi-%J.out -e logfile-scvi-%J.err -q gpu-normal -M100000 -R "select[mem>100000] rusage[mem=100000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    done
done

for p in $(cat /lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/PBMC_merged.normal.subsample500cells.clean_celltypes.txt); do
    d=AR
    # echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --diff_method meld --outpath ${outdir}" | \
    #     bsub -G teichlab -o logfile-cna-%J.out -e logfile-cna-%J.err -q gpu-normal -M150000 -R "select[mem>150000] rusage[mem=150000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --diff_method cna --outpath ${outdir}" | \
        bsub -G teichlab -o logfile-cna-ar-%J.out -e logfile-cna-ar-%J.err -q gpu-normal -M100000 -R "select[mem>100000] rusage[mem=100000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
done

# d=AR
# for p in $(cat /lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/PBMC_merged.normal.subsample500cells.clean_celltypes.txt); do
#     echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --diff_method meld --outpath ${outdir}" | \
#             bsub -G teichlab -o logfile-meld-%J.out -e logfile-meld-%J.err -q gpu-normal -M100000 -R "select[mem>100000] rusage[mem=100000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
# done

d=AR
for p in $(cat /lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/PBMC_merged.normal.subsample500cells.clean_celltypes.txt); do  
    if ! [ -f "$outdir/qPBMC_500cells_demo_perturb_cell_type${p}_queryBatchdataset_id10_1038_s41591_021_01329_2_seed2022/AR_design.scArches_meld.h5ad" ]; then 
        echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --diff_method meld --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-meldAR-%J.out -e logfile-meldAR-%J.err -q gpu-normal -M100000 -R "select[mem>100000] rusage[mem=100000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    fi
    if ! [ -f "$outdir/qPBMC_500cells_demo_perturb_cell_type${p}_queryBatchdataset_id10_1038_s41591_021_01329_2_seed2022/AR_design.scArches_cna.h5ad" ]; then 
        echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --diff_method cna --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-cnaAR-%J.out -e logfile-cnaAR-%J.err -q gpu-normal -M100000 -R "select[mem>100000] rusage[mem=100000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    fi
    done

## Test with more cells

for p in $(cat /lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/tmp/PBMC_merged.normal.subsample1000cells.clean_celltypes.txt); do
    for d in ACR CR; do
        # echo "python run_oor_design_v2.py PBMC_merged.normal.subsample1000cells.clean_celltypes.h5ad ${p} ${d} --diff_method meld --outpath ${outdir}" | \
        #     bsub -G teichlab -o logfile-meld-%J.out -e logfile-meld-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
        # echo "python run_oor_design_v2.py PBMC_merged.normal.subsample1000cells.clean_celltypes.h5ad ${p} ${d} --embedding_method scArches --n_cells_subsample 1000 --outpath ${outdir}" | \
        #     bsub -G teichlab -o logfile-milo-%J.out -e logfile-milo-%J.err -q gpu-normal -M75000 -R "select[mem>75000] rusage[mem=75000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
        echo "python run_oor_design_v2.py PBMC_merged.normal.subsample1000cells.clean_celltypes.h5ad ${p} ${d} --embedding_method scVI --n_cells_subsample 1000 --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-milo-scvi-1k-%J.out -e logfile-milo-scvi-1k-%J.err -q gpu-normal -M75000 -R "select[mem>75000] rusage[mem=75000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
        # echo "python run_oor_design_v2.py PBMC_merged.normal.subsample1000cells.clean_celltypes.h5ad ${p} ${d} --diff_method cna --n_cells_subsample 1000 --outpath ${outdir}" | \
        #     bsub -G teichlab -o logfile-cna-%J.out -e logfile-cna-%J.err -q gpu-normal -M75000 -R "select[mem>75000] rusage[mem=75000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    done
done


## Parse results ##
simdirs=$(ls -d $outdir/qPBMC_500cells_demo_perturb_cell_type*_seed2022)
for s in $simdirs; do
    python parse_oor_design_v2.py $s --embedding_method scArches --diff_method meld
    python parse_oor_design_v2.py $s --embedding_method scArches --diff_method cna
    python parse_oor_design_v2.py $s --embedding_method scArches --diff_method milo
done

simdirs=$(ls -d $outdir/qPBMC_500cells_demo_perturb_cell_type*_seed2022)
for s in $simdirs; do
    python parse_oor_design_v2.py $s --embedding_method scVI --diff_method milo
done