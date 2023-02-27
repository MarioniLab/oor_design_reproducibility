#!/bin/bash

conda activate oor-benchmark
cd /nfs/team205/ed6/bin/diff2atlas/src/4b_COVID_full_design

for m in scArches; do
    # for d in ACR CR; do
    # d=CR
    # echo "python run_COVID_design.py ${d} ${m}" | \
    #         bsub -G teichlab -o logfile-COVID-%J.out -e logfile-COVID-%J.err -q gpu-normal -M75000 -R "select[mem>75000] rusage[mem=75000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    d=ACR
    echo "python run_COVID_design.py ${d} ${m}" | \
            bsub -G team283 -o logfile-COVID-%J.out -e logfile-COVID-%J.err -q gpu-cellgeni -M150000 -R "select[mem>150000] rusage[mem=150000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    # done
done