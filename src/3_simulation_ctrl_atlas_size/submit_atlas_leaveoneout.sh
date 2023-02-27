#!/bin/bash

conda activate oor-benchmark
cd /nfs/team205/ed6/bin/diff2atlas/src/3_simulation_ctrl_atlas_size

outdir=/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/suboptimal_atlas_analysis/
dataset_ids=$(cat /lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/tmp/atlas_dataset_IDs.txt)

for p in natural_killer_cell classical_monocyte naive_thymus_derived_CD8_positive_alpha_beta_T_cell central_memory_CD4_positive_alpha_beta_T_cell effector_memory_CD8_positive_alpha_beta_T_cell naive_thymus_derived_CD4_positive_alpha_beta_T_cell naive_B_cell memory_B_cell; do
    for d in $dataset_ids; do
        echo "python run_atlas_leaveoneout.py ${p} ${d} --outpath ${outdir} --oor_ct_oi"  | \
            bsub -G teichlab -o logfile-leaveoneout-%J.out -e logfile-leaveoneout-%J.err -q gpu-normal -M75000 -R "select[mem>75000] rusage[mem=75000]" -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" 
    done
done

outdir=/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/suboptimal_atlas_analysis/
simdirs=$(ls -d $outdir/leaveoneout* | xargs -n 1 basename)
for s in $simdirs; do
    if ! [ -f "${outdir}/$s/auprc_df.csv" ]; then
    # echo "
      python parse_suboptimal_atlas.py ${s}
    fi
    # " | \
    #         bsub -G teichlab -o logfile-scvi-%J.out -e logfile-scvi-%J.err -q normal -M25000 -R "select[mem>25000] rusage[mem=25000]" 
done

