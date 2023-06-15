#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=20
#SBATCH --time=200:00:00

#make sure you're in the non-polarized directory
cd /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/withmic

CHR=$1
echo "${CHR} PLUMAGE"

#Either subset samples, or filter additionally (here keep only Set B and genotypes with less than 5% missing data). Make sure you use the haploid and masked data for sex chromosomes
bcftools view --samples-file Samples_B-CM-NoBlue.list --force-samples snp/${CHR}.IF-GF-MM2-BP-ANN-HAPMASK.vcf.gz  | \
        bcftools view -e 'F_MISSING>0.05' --min-alleles 2 --max-alleles 2 -Oz -o snp/${CHR}.IF-GF-MM05-BP-ANN-HAPMASK-SetB-CM-NoBlue.vcf.gz
bcftools index snp/${CHR}.IF-GF-MM05-BP-ANN-HAPMASK-SetB-CM-NoBlue.vcf.gz

#Tree, convert vcf to phylip and then iqtree. Based on haploidization there could be invariant sites so the second iqtree command will automatically run on the generated variants-only phylip
python ~/modules/vcf2phylip.py -i snp/${CHR}.IF-GF-MM05-BP-ANN-HAPMASK-SetB-CM-NoBlue.vcf.gz --output-folder mltree
iqtree --redo -T 20 -s mltree/${CHR}.IF-GF-MM05-BP-ANN-HAPMASK-SetB-CM-NoBlue.min4.phy -o 387_CP_MBW_RUS_F --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000
iqtree --redo -T 20 -s mltree/${CHR}.IF-GF-MM05-BP-ANN-HAPMASK-SetB-CM-NoBlue.min4.phy.varsites.phy -o 387_CP_MBW_RUS_F --seqtype DNA -m MFP+ASC -alrt 1000 -B 1000
