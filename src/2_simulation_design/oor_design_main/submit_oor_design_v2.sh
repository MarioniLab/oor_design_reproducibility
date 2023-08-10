#!/bin/bash

conda activate oor-benchmark
cd /nfs/team205/ed6/bin/diff2atlas/src/2_simulation_design/

outdir=/lustre/scratch126/cellgen/team205/ed6/PBMC_CZI_integration_filtered/tmp2/

# --- Run design simulations (removed) --- #

for p in $(cat /lustre/scratch126/cellgen/team205/ed6/PBMC_CZI_integration_filtered/PBMC_merged.normal.subsample500cells.clean_celltypes.txt); do
    ## ACR/CR scArches design ##
    for d in ACR CR; do
        echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --diff_method meld --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-meld-%J.out -e logfile-meld-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
        echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --diff_method cna --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-cna-%J.out -e logfile-cna-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
        echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} ${d} --embedding_method scArches --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-milo-%J.out -e logfile-milo-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    done
    ## CR scVI design ##
    echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} CR --embedding_method scVI --outpath ${outdir}" | \
        bsub -G teichlab -o logfile-scvi-%J.out -e logfile-scvi-%J.err -q gpu-normal -M100000 -R "select[mem>100000] rusage[mem=100000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} CR --diff_method meld --embedding_method scVI --outpath ${outdir}" | \
        bsub -G teichlab -o logfile-scvi-%J.out -e logfile-scvi-%J.err -q gpu-normal -M100000 -R "select[mem>100000] rusage[mem=100000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1"
    echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} CR --diff_method cna --embedding_method scVI --outpath ${outdir}" | \
        bsub -G teichlab -o logfile-scvi-%J.out -e logfile-scvi-%J.err -q gpu-normal -M100000 -R "select[mem>100000] rusage[mem=100000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    ## AR design ##
    echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} AR --diff_method meld --outpath ${outdir}" | \
        bsub -G teichlab -o logfile-meld-ar-%J.out -e logfile-meld-ar-%J.err -q gpu-normal -M150000 -R "select[mem>150000] rusage[mem=150000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} AR --diff_method cna --outpath ${outdir}" | \
        bsub -G teichlab -o logfile-cna-ar-%J.out -e logfile-cna-ar-%J.err -q gpu-normal -M100000 -R "select[mem>100000] rusage[mem=100000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1"
    echo "python run_oor_design_v2.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad ${p} AR --outpath ${outdir}" | \
        bsub -G teichlab -o logfile-ar-%J.out -e logfile-ar-%J.err -q gpu-normal -M100000 -R "select[mem>100000] rusage[mem=100000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1"
done

# --- Run design simulations (mixed effect) --- #

cts="natural_killer_cell naive_B_cell CD14_low_CD16_positive_monocyte central_memory_CD4_positive_alpha_beta_T_cell effector_memory_CD8_positive_alpha_beta_T_cell naive_thymus_derived_CD4_positive_alpha_beta_T_cell memory_B_cell naive_thymus_derived_CD8_positive_alpha_beta_T_cell"
for p in $cts; do
    ## ACR/CR scArches design ##
    for d in ACR CR; do
       echo "python run_oor_mixed_effect.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad classical_monocyte ${p} ${d} --embedding_method scArches --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-mixed-%J.out -e logfile-mixed-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    done
    ## AR design ##
    echo "python run_oor_mixed_effect.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad classical_monocyte ${p} AR --embedding_method scArches --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-mixed-%J.out -e logfile-mixed-%J.err -q gpu-normal -M100000 -R "select[mem>100000] rusage[mem=100000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    ## CR scVI design ##
    echo "python run_oor_mixed_effect.py PBMC_merged.normal.subsample500cells.clean_celltypes.h5ad classical_monocyte ${p} CR --embedding_method scVI --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-mixed-%J.out -e logfile-mixed-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
done


# --- Parse results --- #

simdirs=$(ls -d $outdir/qPBMC_500cells_demo_perturb_cell_type*_seed2022)
for s in $simdirs; do
    echo "python parse_oor_design_v2.py $s --diff_method milo" | \
        bsub -G teichlab -o logfile-parsemilo-%J.out -e logfile-parsemilo-%J.err -q normal -M25000 -R "select[mem>25000] rusage[mem=25000]"
    echo "python parse_oor_design_v2.py $s --diff_method cna" | \
        bsub -G teichlab -o logfile-parsecna-%J.out -e logfile-parsecna-%J.err -q normal -M25000 -R "select[mem>25000] rusage[mem=25000]" 
    echo "python parse_oor_design_v2.py $s --diff_method meld" | \
        bsub -G teichlab -o logfile-parsemeld-%J.out -e logfile-parsemeld-%J.err -q normal -M25000 -R "select[mem>25000] rusage[mem=25000]"
done

simdirs=$(ls -d $outdir/qPBMC_500cells_demo_mixed_cell_type*_seed2022)
for s in $simdirs; do
    python parse_oor_design_v2.py $s --diff_method milo --mixed_effect
done