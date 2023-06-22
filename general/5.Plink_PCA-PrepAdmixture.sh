#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=30000mb
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

# Merge VCFs found in Autosomes.list
bcftools concat --file-list Autosomes.list --threads 10 -Oz -o autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.vcf.gz
bcftools index --threads 10 autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.vcf.gz

GROUP=Samples_C-CM-NoBlue.inds

# CHR_W, LD-pruning doesn't make sense 
~/modules/plink --keep ${GROUP}.inds --pca --mac 2 --threads 20 --vcf snp/chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK.vcf.gz --allow-extra-chr --set-missing-var-ids @:# --make-bed --double-id --out autosomes/chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK.${GROUP}

# Now for autosomes
VCF=autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.vcf.gz

# AUTOSOMES: LD-prune, remove SNPs missing genotypes 20%
~/modules/plink --keep ${GROUP} --indep-pairwise 50 10 0.1 --mac 2 --threads 10 --vcf ${VCF} \
        --geno 0.2 --keep-allele-order --allow-extra-chr --set-missing-var-ids @:# --make-bed --double-id --out autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.SETC-P50.10.1
~/modules/plink --keep ${GROUP} --extract autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.SETC-P50.10.1.prune.in --recode vcf bgz \
        --pca --geno 0.2 --mac 2 --threads 10 --vcf ${VCF} --allow-extra-chr --set-missing-var-ids @:# \
        --double-id --keep-allele-order --make-bed --out autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.SETC-P50.10.1-LD
bcftools index autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.SETC-P50.10.1-LD.vcf.gz
