#!/bin/bash

conda activate oor-benchmark
cd /nfs/team205/ed6/bin/diff2atlas/scripts

outdir=/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/tmp/
dirs=$(dirname $outdir/*/acr_design.h5ad)
# dirs=$(cat missing_cts_ctrl_size.txt)

## Selecting cell types with high TPR
for d in $dirs; do
    for nc in $(seq 3 12); do
        for s in $(seq 12348 12349); do 
            echo "python run_oor_ctrl_size.py ${d}/ ${nc} 7 ${s}" | \
                bsub -G teichlab -o logfile-ctrlsize-%J.out -e logfile-ctrlsize-%J.err -q gpu-normal -M50000 -R "select[mem>50000] rusage[mem=50000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
        done
    done    
done

# ## On GPU cloud ##
# for nc in $(seq 3 12); do
#         for s in $(seq 12345 12347); do 
#             python run_oor_ctrl_size.py ${d}/ ${nc} 7 ${s}                
#         done
#     done

## To merge results
awk '(NR == 1) || (FNR > 1)' ${outdir}/*/ctrl_size_analysis_nquery7/*/*_results.csv > ${outdir}/ctrl_size_nquery7_results.csv
