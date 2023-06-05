#!/bin/bash
#SBATCH --job-name=sort
#SBATCH --output=logs/sort.%A_%a.out
#SBATCH --error=logs/sort.%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --mem=110GB
#SBATCH --cpus-per-task=6
#SBATCH --array=1-xx
#sbatch sortindex_array.sh sample_lists/kokanee_files.txt

module load samtools

cd $ANNA
SAMPLELIST=$1
SAMPLE=`cat $SAMPLELIST | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

samtools view -F 0x4 -b bamfiles/${SAMPLE}_2sort.bam | samtools sort -@6 - -o goodbam/$SAMPLE.bam

samtools index goodbam/$SAMPLE.bam
