#!/bin/bash
#SBATCH --job-name=filterVcf
#SBATCH --mem 20G
#SBATCH --time 12:00:00
#SBATCH -o logs/filtervcf_10mb.%A_%a.log
#SBATCH --cpus-per-task=1
#SBATCH --array=1-161
#sbatch scripts/filtervcf_chr_intervals_allkokanee.sh 

module load bcftools
module load vcftools

cd $ANNA

REGION=`cat reference/onerka_chr_10mb_intervals.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
OUTPUTPATH="variant_calling/allkokanee/bcftools"

vcftools --vcf "$OUTPUTPATH"/allkokanee_"$SLURM_ARRAY_TASK_ID".tmp.vcf \
    --minQ 30 \
    --minGQ 20 \
    --minDP 3 \
    --maxDP 50 \
    --max-alleles 2 \
    --min-alleles 2 \
    --max-missing 0.7 \
    --remove-indels \
    --recode \
    --stdout > "$OUTPUTPATH"/allkokanee_"$SLURM_ARRAY_TASK_ID".vcf
