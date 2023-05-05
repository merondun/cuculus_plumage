#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=30000mb
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

GROUP=$1
#CHR_W, LD-prune, but do pca on both
~/modules/plink --keep ${GROUP}.inds --pca --indep-pairwise 50 5 0.1 --mac 2 --threads 20 --vcf snp/chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK.vcf.gz --allow-extra-chr --set-missing-var-ids @:# --make-bed --double-id --out autosomes/chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK.${GROUP}
~/modules/plink --keep ${GROUP}.inds --extract autosomes/chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK.${GROUP}.prune.in --recode vcf bgz --pca --mac 2 --threads 20 --vcf snp/chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK.vcf.gz --allow-extra-chr --set-missing-var-ids @:# --double-id --make-bed --out autosomes/chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK-LD.${GROUP}

#AUTOSOMES
~/modules/plink --keep ${GROUP}.inds --pca --indep-pairwise 50 5 0.1 --mac 2 --threads 20 --vcf snp/chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK.vcf.gz --allow-extra-chr --set-missing-var-ids @:# --make-bed --double-id --out autosomes/chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK.${GROUP}
~/modules/plink --keep ${GROUP}.inds --extract autosomes/chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK.${GROUP}.prune.in --recode vcf bgz --pca --mac 2 --threads 20 --vcf snp/chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK.vcf.gz --allow-extra-chr --set-missing-var-ids @:# --double-id --make-bed --out autosomes/chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK-LD.${GROUP}
