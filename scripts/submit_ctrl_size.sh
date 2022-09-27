#!/bin/bash

conda activate scvi-env
cd /nfs/team205/ed6/bin/diff2atlas/scripts

outdir=/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/

## Selecting cell types with high TPR

for ct in natural_killer_cell classical_monocyte central_memory_CD4_positive_alpha_beta_T_cell effector_memory_CD8_positive_alpha_beta_T_cell; do
    # for d in $(ls -d $outdir/qPBMC*/); do
    d=${outdir}/qPBMC_500cells_demo_perturb_cell_type${ct}_queryBatchdataset_id10_1038_s41591_021_01329_2_seed2022
    # for nq in $(seq 3 15); do
    for nc in $(seq 5 12); do
        for s in $(seq 12345 12347); do 
            echo "python run_ctrl_size.py $d $nc 5 $s"  | \
            bsub -G teichlab -o logfile-merge-%J.out -e logfile-merge-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
        done
    done
done

# for d in $(ls -d $outdir*/); do
#     echo "python run_milo.py ${d} --prop target" | \
#         bsub -G teichlab -o logfile-milo-%J.out -e logfile-milo-%J.err -q gpu-normal -M30000 -R "select[mem>30000] rusage[mem=30000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
# done

