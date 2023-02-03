#!/bin/bash

conda activate oor-benchmark
cd /nfs/team205/ed6/bin/diff2atlas/src/3_simulation_ctrl_atlas_size/

outdir=/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/tmp/
cts="natural_killer_cell classical_monocyte CD14_low_CD16_positive_monocyte central_memory_CD4_positive_alpha_beta_T_cell effector_memory_CD8_positive_alpha_beta_T_cell naive_thymus_derived_CD4_positive_alpha_beta_T_cell"
dirs=$(for c in $cts; do dirname $outdir/*${c}*/ACR_design.scArches_milo.h5ad; done)

for d in $dirs; do
    for nc in $(seq 3 12); do
        for nq in $(seq 5 2 9); do
            for s in $(seq 12345 12349); do 
                echo "python run_oor_ctrl_size.py ${d}/ ${nc} ${nq} ${s}" | \
                    bsub -G teichlab -o logfile-ctrlsize-%J.out -e logfile-ctrlsize-%J.err -q gpu-normal -M75000 -R "select[mem>75000] rusage[mem=75000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
            done
        done    
    done
done
 
## To merge results
awk '(NR == 1) || (FNR > 1)' ${outdir}/*/ctrl_size_analysis_nquery7/*/*_results.csv > ${outdir}/ctrl_size_nquery7_results.csv