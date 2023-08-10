#!/bin/bash

conda activate scvi-env
cd /nfs/team205/ed6/bin/diff2atlas/src/1_PBMC_data_preprocessing

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

# Prep variations of query dataset
echo "python split_PBMC_dataset.py 10_1038_s41591_021_01329_2 normal --n_cells_sample 1000" | \
    bsub -G teichlab -o logfile-query1k.out -e logfile-query1k.err -q normal -M50000 -R "select[mem>50000] rusage[mem=50000]" 
echo 'python split_PBMC_dataset.py 10_1038_s41591_021_01329_2 normal --n_cells_sample 1500' |
    bsub -G teichlab -o logfile-query15k.out -e logfile-query15k.err -q normal -M50000 -R "select[mem>50000] rusage[mem=50000]" 