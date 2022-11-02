#!/bin/bash

conda activate scvi-env
cd /nfs/team205/ed6/bin/query2reference_uncertainty/notebooks/PBMC_metanalysis

dataset_ids=$(cat /nfs/team205/ed6/data/PBMC_CZI_integration_filtered/PBMC_sample_metadata.normal.csv | tail -n +2 | cut -f 9 -d ','| sort | uniq)
for d in $dataset_ids; do
    echo "python split_PBMC_dataset.py ${d} normal" | \
        bsub -G teichlab -o logfile-${d}-normal.out -e logfile-${d}-normal.err -q normal -M50000 -R "select[mem>50000] rusage[mem=50000]" 
    done

dataset_ids=$(cat /nfs/team205/ed6/data/PBMC_CZI_integration_filtered/PBMC_sample_metadata.COVID.csv | tail -n +2 | cut -f 9 -d ','| sort | uniq)
for d in $dataset_ids; do
    echo "python split_PBMC_dataset.py ${d} COVID" | \
        bsub -G teichlab -o logfile-${d}-covid.out -e logfile-${d}-covid.err -q normal -M50000 -R "select[mem>50000] rusage[mem=50000]" 
    done

dataset_ids=$(cat /nfs/team205/ed6/data/PBMC_CZI_integration_filtered/PBMC_sample_metadata.lupus.csv | tail -n +2 | cut -f 9 -d ','| sort | uniq)
for d in $dataset_ids; do
    echo "python split_PBMC_dataset.py ${d} lupus" | \
        bsub -G teichlab -o logfile-${d}-lupus.out -e logfile-${d}-lupus.err -q normal -M50000 -R "select[mem>50000] rusage[mem=50000]" 
    done
    
