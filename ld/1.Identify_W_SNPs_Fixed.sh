# Set populations
p1=HH
p2=HG
GROUP=HH_HG
CHR=chr_W

# Create pop file with all samples for subsetting 
cat populations/${p1}.pop.list populations/${p2}.pop.list > divergence/allbp/${CHR}_${GROUP}.inds

# Subset vcf, retain high quality SNPs
bcftools view --force-samples --samples-file divergence/allbp/${CHR}_${GROUP}.inds snp/${CHR}.IF.ANN-HAP.vcf.gz -Ou | bcftools view --min-alleles 2 --max-alleles 2 -i 'F_MISSING<0.01' --min-af 0.2 -Oz -o divergence/allbp/${CHR}.${GROUP}.vcf.gz
bcftools index divergence/allbp/${CHR}.${GROUP}.vcf.gz
bcftools index -n divergence/allbp/${CHR}.${GROUP}.vcf.gz #17116 SNPs

# Calculate FST
~/modules/vcftools/bin/vcftools --haploid --gzvcf divergence/allbp/${CHR}.${GROUP}.vcf.gz --weir-fst-pop populations/${p1}.pop.list --weir-fst-pop populations/${p2}.pop.list --out divergence/allbp/${CHR}_${GROUP}
#Weir and Cockerham weighted Fst estimate: 0.91282

# Grab 5 random W SNPs which show FST 1 with very low missing data, these will all show the same pattern
awk '$3 == 1' divergence/allbp/chr_W_HH_HG.weir.fst | shuf -n 5 | awk '{OFS="\t"}{print $1, $2-1, $2}' | bedtools sort -i - > WLD/W_SNPS_FST1.bed
bcftools view -R WLD/W_SNPS_FST1.bed snp/${CHR}.IF.ANN-HAP.vcf.gz -Oz -o WLD/W_SNPS_FST1.vcf.gz
