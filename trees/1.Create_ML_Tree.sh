#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --time=200:00:00

# The script is taking an argument which would be the chromosome number
CHR=$1
echo "${CHR} PLUMAGE"

# This block of code filters the VCF files based on the provided samples file
# The output is a VCF file for the chromosome with genotypes having less than 5% missing data
bcftools view --samples-file Samples_B-CM-NoBlue.list --force-samples snp/${CHR}.IF-GF-MM2-BP-ANN-HAPMASK.vcf.gz  | \
        bcftools view -e 'F_MISSING>0.05' --min-alleles 2 --max-alleles 2 -Oz -o snp/${CHR}.IF-GF-MM05-BP-ANN-HAPMASK-SetB-CM-NoBlue.vcf.gz
bcftools index snp/${CHR}.IF-GF-MM05-BP-ANN-HAPMASK-SetB-CM-NoBlue.vcf.gz

# The VCF file is then converted into a phylip file for phylogenetic analysis using a custom script (vcf2phylip.py)
python ~/modules/vcf2phylip.py -i snp/${CHR}.IF-GF-MM05-BP-ANN-HAPMASK-SetB-CM-NoBlue.vcf.gz --output-folder mltree

# iqtree is then used to generate a phylogenetic tree based on the converted VCF data
# The script runs iqtree twice - once for the full phylip file and once for a 'variant sites only' version
# The -o flag denotes the outgroup taxon for the tree. 
iqtree --redo -T 20 -s mltree/${CHR}.IF-GF-MM05-BP-ANN-HAPMASK-SetB-CM-NoBlue.min4.phy -o 387_CP_MBW_RUS_F --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000
iqtree --redo -T 20 -s mltree/${CHR}.IF-GF-MM05-BP-ANN-HAPMASK-SetB-CM-NoBlue.min4.phy.varsites.phy -o 387_CP_MBW_RUS_F --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000
