#!/bin/bash

conda activate oor-benchmark
cd /nfs/team205/ed6/bin/diff2atlas/scripts

outdir=/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/suboptimal_atlas_analysis/

for p in conventional_dendtritic_cell natural_killer_cell platelet; do
    for d in ACR AR; do
        for n in $(seq 3 11); do
        echo "python run_suboptimal_atlas.py ${p} ${d} ${n} --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-suboptimal-%J.out -e logfile-suboptimal-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
        done
    done
done

echo "python run_suboptimal_atlas.py conventional_dendtritic_cell CR 3 --outpath ${outdir}" | \
            bsub -G teichlab -o logfile-suboptimal-%J.out -e logfile-suboptimal-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 

## With OOR population
for p in natural_killer_cell classical_monocyte CD14_low_CD16_positive_monocyte central_memory_CD4_positive_alpha_beta_T_cell effector_memory_CD8_positive_alpha_beta_T_cell naive_thymus_derived_CD4_positive_alpha_beta_T_cell naive_B_cell memory_B_cell; do
    for d in ACR AR; do
        # for n in $(seq 1 12); do
        for n in $(seq 1 2); do
        echo "python run_suboptimal_atlas.py ${p} ${d} ${n} --outpath ${outdir} --oor_ct_oi"  | \
            bsub -G teichlab -o logfile-suboptimal-%J.out -e logfile-suboptimal-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
        done
    done
done

# With CR design
for p in natural_killer_cell classical_monocyte CD14_low_CD16_positive_monocyte central_memory_CD4_positive_alpha_beta_T_cell effector_memory_CD8_positive_alpha_beta_T_cell naive_thymus_derived_CD4_positive_alpha_beta_T_cell naive_B_cell memory_B_cell; do
    echo "python run_suboptimal_atlas.py ${p} CR 3 --outpath ${outdir} --oor_ct_oi" | \
            bsub -G teichlab -o logfile-suboptimal-%J.out -e logfile-suboptimal-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
        done

echo "python run_suboptimal_atlas.py natural_killer_cell CR 3 --outpath ${outdir} --oor_ct_oi" | \
            bsub -G teichlab -o logfile-suboptimal-%J.out -e logfile-suboptimal-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
echo "python run_suboptimal_atlas.py platelet CR 3 --outpath ${outdir} --oor_ct_oi" | \
            bsub -G teichlab -o logfile-suboptimal-%J.out -e logfile-suboptimal-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 

conda activate oor-benchmark
cd /nfs/team205/ed6/bin/diff2atlas/scripts
outdir=/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/suboptimal_atlas_analysis/
simdirs=$(ls -d $outdir/suboptimal_atlas_queryBatchdataset_id10_1038_s41591_021_01329_2_seed2022_keep*)
for s in $simdirs; do
    # echo "
    python parse_oor_design.py ${s}
    # " | \
    #         bsub -G teichlab -o logfile-scvi-%J.out -e logfile-scvi-%J.err -q normal -M25000 -R "select[mem>25000] rusage[mem=25000]" 
done

