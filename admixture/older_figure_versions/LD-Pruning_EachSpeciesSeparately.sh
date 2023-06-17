#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

mkdir ld_pruned
VCF=autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.vcf.gz

#LD-prune, remove SNPs missing genotypes 20%
for GROUP in 'CO' 'CC'; do 
~/modules/plink --keep ${GROUP}.inds --indep-pairwise 50 10 0.1 --mac 2 --threads 10 --vcf ${VCF} \
	--geno 0.2 --keep-allele-order --allow-extra-chr --set-missing-var-ids @:# --make-bed --double-id --out ld_pruned/Autosomes.IF-GF-MM2-BP-ANN-AC2.${GROUP}
~/modules/plink --keep ${GROUP}.inds --extract ld_pruned/Autosomes.IF-GF-MM2-BP-ANN-AC2.${GROUP}.prune.in --recode vcf bgz \
	--pca --geno 0.2 --mac 2 --threads 10 --vcf ${VCF} --allow-extra-chr --set-missing-var-ids @:# \
	--double-id --keep-allele-order --make-bed --out ld_pruned/Autosomes.IF-GF-MM2-BP-ANN-AC2.${GROUP}-LD
bcftools index ld_pruned/Autosomes.IF-GF-MM2-BP-ANN-AC2.${GROUP}-LD.vcf.gz
done 

#merge the ld-pruned files and then create a full bed again from the LD pruned sites 
bcftools merge --threads 10 -Oz -o ld_pruned/Autosomes.IF-GF-MM2-BP-ANN-AC2.All-LD.vcf.gz ld_pruned/Autosomes.IF-GF-MM2-BP-ANN-AC2.CC-LD.vcf.gz ld_pruned/Autosomes.IF-GF-MM2-BP-ANN-AC2.CO-LD.vcf.gz
~/modules/plink --threads 10 --vcf ld_pruned/Autosomes.IF-GF-MM2-BP-ANN-AC2.All-LD.vcf.gz \
	--allow-extra-chr --set-missing-var-ids @:# --make-bed --double-id --out ld_pruned/Autosomes.IF-GF-MM2-BP-ANN-AC2.All-LD
