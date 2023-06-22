#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

VCF=autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.vcf.gz
GROUP=Samples_C-CM-NoBlue.inds

#LD-prune, remove SNPs missing genotypes 20%
~/modules/plink --keep ${GROUP} --indep-pairwise 50 10 0.1 --mac 2 --threads 10 --vcf ${VCF} \
	--geno 0.2 --keep-allele-order --allow-extra-chr --set-missing-var-ids @:# --make-bed --double-id --out autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.SETC-P50.10.1
~/modules/plink --keep ${GROUP} --extract autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.SETC-P50.10.1.prune.in --recode vcf bgz \
	--pca --geno 0.2 --mac 2 --threads 10 --vcf ${VCF} --allow-extra-chr --set-missing-var-ids @:# \
	--double-id --keep-allele-order --make-bed --out autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.SETC-P50.10.1-LD
bcftools index autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.SETC-P50.10.1-LD.vcf.gz
