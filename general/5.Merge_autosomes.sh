#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=50000mb
#SBATCH --cpus-per-task=30
#SBATCH --time=72:00:00

#LD pruned set
bcftools concat --file-list Autosomes.list --threads 10 -Oz -o autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.vcf.gz
bcftools index --threads 10 autosomes/Autosomes.IF-GF-MM2-BP-ANN-AC2.vcf.gz
