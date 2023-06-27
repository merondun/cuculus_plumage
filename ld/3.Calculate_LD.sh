#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

# Loop through targets 
CHR=$1

echo "WORKING ON CHROMOSOME: ${CHR}"	
# Extract SNPs on current autosome
plink --bfile WLD/allchromosomes --chr ${CHR} --make-bed --keep-allele-order --allow-extra-chr --set-missing-var-ids @:# --double-id --threads 10 --out WLD/${CHR}.tmp

# Calculate LD between SNPs on chr_W and current autosome
plink --bfile WLD/chr_W --bmerge WLD/${CHR}.tmp --keep-allele-order --allow-extra-chr --set-missing-var-ids @:# --double-id --threads 10 --r2 inter-chr gz --out WLD/chr_W__${CHR}
zcat WLD/chr_W__${CHR}.ld.gz | awk '$1 != $4' > WLD/chr_W__${CHR}.ldout
#rm WLD/chr_W__${CHR}.ld.gz

