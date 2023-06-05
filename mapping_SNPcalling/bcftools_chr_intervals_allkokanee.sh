#!/bin/bash
#SBATCH --job-name=bcftools_chr
#SBATCH --output=logs/bcftools_chr.%A_%a.out
#SBATCH --error=logs/bcftools_chr.%A_%a.err
#SBATCH --time=1-00:00
#SBATCH --nodes=1
#SBATCH --mem=110GB
#SBATCH --cpus-per-task=1
#SBATCH --array=1-161
#sbatch scripts/bcftools_chr_intervals_allkokanee.sh 
module load bcftools
module load vcftools

cd $ANNA

REGION=`cat reference/onerka_chr_10mb_intervals.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
GENOME="reference/onerka_chr.fa"
OUTPUTPATH="variant_calling/allkokanee/bcftools"
BAM="sample_lists/kokanee_allbams.txt"

bcftools mpileup -Ou -f $GENOME -b $BAM -q 30 -Q 30 -r $REGION -I -a AD,DP,SP,ADF,ADR -d 50 | bcftools call -G - -mv -Ov > "$OUTPUTPATH"/allkokanee_"$SLURM_ARRAY_TASK_ID".tmp.vcf
