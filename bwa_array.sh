#!/bin/bash
#SBATCH --job-name=BWA_MEM
#SBATCH --output=logs/BWA_MEM.%A_%a.out
#SBATCH --error=logs/BWA_MEM.%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --mem=110GB
#SBATCH --cpus-per-task=6
#SBATCH --array=1-xxx
#sbatch bwa_array.sh sample_lists/kokanee_files.txt

module load bwa
module load samtools
module load samblaster

cd $ANNA
SAMPLELIST=$1
SAMPLE=`cat $SAMPLELIST | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

bwa mem -t 6 reference/onerka.fa \
trimmed/${SAMPLE}_trimmed_1.fastq.gz \
trimmed/${SAMPLE}_trimmed_2.fastq.gz \
| samblaster --removeDups | samtools view -S -h -b -@6 -o bamfiles/${SAMPLE}_2sort.bam
