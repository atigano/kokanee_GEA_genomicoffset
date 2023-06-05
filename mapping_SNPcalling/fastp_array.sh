#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --output=fastp.%A_%a.out
#SBATCH --error=fastp.%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --mem=110GB
#SBATCH --cpus-per-task=1
#SBATCH --array=1-xxx
#sbatch scripts/fastp_array.sh sample_lists/kokanee_files.txt
module load fastp

cd $ANNA
SAMPLELIST=$1
SAMPLE=`cat $SAMPLELIST | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

fastp -i renamed/${SAMPLE}_1.fastq.gz -I renamed/${SAMPLE}_2.fastq.gz -o trimmed/${SAMPLE}_trimmed_1.fastq.gz -O trimmed/${SAMPLE}_trimmed_2.fastq.gz --adapter_sequence AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --adapter_sequence_r2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG -G -h ${SAMPLE}.html
