#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

autosomes=autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.vcf.gz
chrW=WLD/W_SNPS_FST1.vcf.gz
chrZ=snp/chr_Z.IF-GF-MM2-BP-ANN-AC2-HAPMASK.vcf.gz
GROUP=populations/HUNGARY.pop.ind

#We will simply randomly subset 10 high quality W SNPs which show the phenotype 
# Convert VCF to PLINK format for all SNPs
plink --vcf ${autosomes} --keep ${GROUP} --keep-allele-order --allow-extra-chr --set-missing-var-ids @:# --double-id --threads 10 --geno 0.2 --mac 2 --make-bed --out WLD/autosomes
plink --vcf ${chrW} --keep ${GROUP} --keep-allele-order --allow-extra-chr --set-missing-var-ids @:# --double-id --threads 10 --make-bed --out WLD/chr_W
plink --vcf ${chrZ} --keep ${GROUP} --keep-allele-order --allow-extra-chr --set-missing-var-ids @:# --double-id --threads 10 --geno 0.2 --mac 2 --make-bed --out WLD/chr_Z
plink --bfile WLD/autosomes --bmerge WLD/chr_Z --keep-allele-order --allow-extra-chr --set-missing-var-ids @:# --double-id --threads 10  --make-bed --out WLD/allchromosomes
