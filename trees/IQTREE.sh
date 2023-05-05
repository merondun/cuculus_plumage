#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=3
#SBATCH --time=48:00:00

#ML TREE
CHR=$1
echo "${CHR} TREE"

#ensure that there's no SNPs which fall in our coverage mask which aren't likely haploid W based on male coverage, also only retain samples in Retain.list 
bedtools subtract -header -a snp/${CHR}.IF-GF-MM2-BP-ANN-HAP.vcf.gz -b /dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/${CHR}.consmask | \
  bcftools view --force-samples --samples-file Retain.list -Ou | \
  bcftools view --min-alleles 2 --max-alleles 2 -Oz -o mltree/${CHR}_WithChinaMask.vcf.gz
#Tree
python ~/modules/vcf2phylip.py -i mltree/${CHR}_WithChinaMask.vcf.gz --output-folder mltree
#do it 
iqtree --redo -T 3 -s mltree/${CHR}_WithChinaMask.min4.phy -o 387_CP_MBW_RUS_F --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000

