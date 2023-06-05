#!/bin/bash
#SBATCH --job-name=plink_ldprune
#SBATCH --mem 50G
#SBATCH --time 12:00:00
#SBATCH -o logs/%x-%j.log
#SBATCH --cpus-per-task=1
#sbatch scripts/plink_ldprune_allkokanee_norep_0.5.sh

module load plink/1.9b_6.21-x86_64
module load vcftools
cd $ANNA/variant_calling/allkokanee/bcftools
awk 'BEGIN{OFS="\t"} !/#/ {sub(/\./, $1"_"$2, $3)}1' allkokanee_nomiss30_norep.vcf > allkokanee_nomiss30_norep_annot.vcf
plink --vcf allkokanee_nomiss30_norep_annot.vcf --make-bed --allow-extra-chr --double-id --no-sex --no-pheno --out allkokanee_nomiss30_norep_annot
 
WINDOW=200kb
SNP=100
R2=0.5
plink --bed allkokanee_nomiss30_norep_annot.bed \
--bim allkokanee_nomiss30_norep_annot.bim \
--fam allkokanee_nomiss30_norep_annot.fam \
--indep-pairwise $WINDOW $SNP $R2 --allow-extra-chr \
--out allkokanee_nomiss30_norep_annot_ld

vcftools --vcf allkokanee_nomiss30_norep_annot.vcf --exclude allkokanee_nomiss30_norep_annot_ld.prune.out --recode --out allkokanee_nomiss30_norep_annot_ld
