#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=20000mb
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00

#This will apply some hard INFO filters, require 2x for sex chrom genotype and 4x for autosomes. Set Z/W/mtDNA genotypes to most frequent allele UNLESS it appers to be truly heterozygous, in which case make it missing (binomial test AD<1e-5 
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
cd /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes

mkdir vcfqc snp allsites phased lostruct lostruct/input

CHR=$1
echo "${CHR} PLUMAGE"

# Perform initial variant calling filtering

# Display status message
echo "FILTER VCF"

# Use bcftools to perform sample-wise filtering on the vcf file. 
# The resulting vcf will only contain samples listed in the AllSamples.list file.
bcftools view --force-samples --samples-file AllSamples.list -Oz ../../merged/${CHR}.vcf.gz -o snp/${CHR}.vcf.gz

# Index the vcf file for faster access
bcftools index --threads 10  snp/${CHR}.vcf.gz

# Calculate some statistics for the next round of filtering
# Count the number of samples in the vcf file
SAMPN=$(bcftools query -l snp/${CHR}.vcf.gz | wc -l)

# Compute the median coverage depth across all sites in the vcf file
AVGDP=$(bcftools query -f '%DP\n' snp/${CHR}.vcf.gz | datamash median 1 | datamash round 1)

# Calculate thresholds for filtering based on depth of coverage
DPLOW=$(($AVGDP/3))
DPHI=$(($AVGDP*2))

# Perform a second round of filtering based on various quality metrics
bcftools view --types snps --threads 10 -e "QUAL < 20 || INFO/DP > $DPHI || INFO/DP < $SAMPN || MQ < 30 || RPBZ < -3 || RPBZ > 3" -Oz -o snp/${CHR}.IF.vcf.gz snp/${CHR}.vcf.gz
bcftools index --threads 10 snp/${CHR}.IF.vcf.gz

# Determine the minimum depth of coverage based on the type of chromosome
if [[ $CHR = 'chr_W' || $CHR == 'chr_Z' ]]
then
    MINDP=2 # Lower depth for sex chromosomes
else
    MINDP=4  # Higher depth for autosomes and mitochondria
fi

# Different processing steps for sex chromosomes and mitochondria versus autosomes
if [[ $CHR = 'chr_W' || $CHR == 'chr_MT' || $CHR == 'chr_Z' ]]
then
    # For haploid chromosomes, convert low-coverage and heterozygous sites to missing
    echo "HAPLOID, SETTING HETEROZYGOUS SITES TO MISSING"
    bcftools +setGT -Ou -o - snp/${CHR}.IF.vcf.gz -- -t q -i "FMT/DP < ${MINDP}" -n "./." | \
    bcftools +setGT -Ou -- -t "b:AD<1e-5" -n "./." | \
    bcftools +setGT -Oz -o snp/${CHR}.IF-GF.vcf.gz -- --target-gt q --new-gt M -i 'GT=="het"'
    bcftools index --threads 10 snp/${CHR}.IF-GF.vcf.gz
else
    # For diploid chromosomes, convert low-coverage sites to missing
    echo "DIPLOID, SETTING LOW COVERAGE TO MISSING"
    bcftools +setGT -Oz -o snp/${CHR}.IF-GF.vcf.gz snp/${CHR}.IF.vcf.gz -- -t q -i "FMT/DP < ${MINDP}" -n "./."
    bcftools index --threads 10 snp/${CHR}.IF-GF.vcf.gz
fi

#FORMAT SNP summary stats
echo -e "CHROM\tID\tnREF\tnALT\tnHET\tnTs\tnTv\tavgDP\tSingletons\tMissing_Sites" > vcfqc/${CHR}_sample.stats
bcftools stats -S- snp/${CHR}.IF-GF.vcf.gz | grep 'PSC' | tr ' ' '_' | awk -v c=${CHR} '{OFS="\t"}{print c, $3,$4,$5,$6,$7,$8,$10,$11,$14}' | sed '1,2d' >> vcfqc/${CHR}_sample.stats
