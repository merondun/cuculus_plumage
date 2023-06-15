#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=2
#SBATCH --time=72:00:00

CHR=$1

mkdir divergence divergence/allbp
popdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/withmic/populations

#annotate, take the INFO filtered SNPs, not the max-missing set yet 
java -Xmx25g -jar ~/modules/snpEff/snpEff.jar -v -stats vcfqc/${CHR}_RAW_annotation.html bCucCan1 snp/${CHR}.IF-GF.vcf.gz | bgzip -c > snp/${CHR}.IF-GF-ANN.vcf.gz

#Change ploidy
bcftools +fixploidy snp/${CHR}.IF-GF-ANN.vcf.gz -- -f 1 | bcftools view -Oz -o snp/${CHR}.IF-GF-ANN-HAP.vcf.gz

#remove mask
bedtools subtract -header -a snp/${CHR}.IF-GF-ANN-HAP.vcf.gz -b /dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/${CHR}.consmask | \
  bcftools view --min-alleles 2 --max-alleles 2 -Oz -o snp/${CHR}.IF-GF-ANN-HAPMASK.vcf.gz
bcftools index snp/${CHR}.IF-GF-ANN-HAPMASK.vcf.gz

for GROUP in $(cat Pairwise_Plumage.list); do
echo "WORKING ON CHR: ${CHR} and ${GROUP}"

p1=$(echo ${GROUP} | sed 's/_.*//g')
p2=$(echo ${GROUP} | sed 's/.*_//g')

cat populations/${p1}.pop.list populations/${p2}.pop.list > divergence/allbp/${CHR}_${GROUP}.inds
bcftools view --force-samples --samples-file divergence/allbp/${CHR}_${GROUP}.inds snp/${CHR}.IF-GF-ANN-HAPMASK.vcf.gz -Ou | \
  bcftools view --min-alleles 2 --max-alleles 2 -i 'F_MISSING<0.1' -Oz -o divergence/allbp/${CHR}.${GROUP}.vcf.gz

#calculate fst
~/modules/vcftools/bin/vcftools --haploid --gzvcf divergence/allbp/${CHR}.${GROUP}.vcf.gz --weir-fst-pop populations/${p1}.pop.list --weir-fst-pop populations/${p2}.pop.list --out divergence/allbp/${CHR}_${GROUP}
awk '{OFS="\t"}{print $1, $2-1, $2, $3}' divergence/allbp/${CHR}_${GROUP}.weir.fst | grep -v 'nan' > divergence/allbp/${CHR}_${GROUP}.bed

if [[ $CHR = 'chr_MT' ]]
then
  #intersect with VCF file
  bedtools intersect -a divergence/allbp/${CHR}_${GROUP}.bed -b snp/${CHR}.IF-GF-ANN.vcf.gz -wao | \
    awk '{OFS="\t"}{if ($3 == $6) print $1, $2, $3, $4, $12}' | sed 's/DP.*ANN=.|//g' | tr '|' '\t' | \
    awk -v g=${GROUP} '{OFS="\t"}{print $1, $2, $3, $4, $5, "FILLER", g}' > divergence/allbp/${CHR}_${GROUP}.SNPs.bed
else
  #intersect with VCF file
  bedtools intersect -a divergence/allbp/${CHR}_${GROUP}.bed -b snp/${CHR}.IF-GF-ANN.vcf.gz -wao | \
    awk '{OFS="\t"}{if ($3 == $6) print $1, $2, $3, $4, $12}' | sed 's/DP.*ANN=.|//g' | tr '|' '\t' | \
    awk -v g=${GROUP} '{OFS="\t"}{print $1, $2, $3, $4, $5, $10, g}' > divergence/allbp/${CHR}_${GROUP}.SNPs.bed
fi

bedtools intersect -c -a divergence/allbp/${CHR}_${GROUP}.SNPs.bed -b /dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/${CHR}.consmask > divergence/allbp/${CHR}_${GROUP}.FIXED.bed

done

