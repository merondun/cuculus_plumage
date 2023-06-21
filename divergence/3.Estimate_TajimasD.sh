#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --time=38:00:00

# Create necessary directories
mkdir divergence divergence/tajima divergence/out

# This script accepts two parameters: Group name and iteration number
GROUP=$1
IT=$2

# Define the directory for population data
popdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/populations

# Extract distinct populations from popfile
n1=$(awk '{print $2}' $popdir/${GROUP}.popfile | sort | uniq | head -n 1)
n2=$(awk '{print $2}' $popdir/${GROUP}.popfile | sort | uniq | tail -n 1)

# Depending on the group, subsample individuals from each population differently
if [[ ${GROUP} == 'HUNGARY' ]]; then
        echo "HUNGARIAN POPULATION, MAKING RUFOUS 60% GREY 40%"
        awk -v p=${n1} '$2 == p' $popdir/${GROUP}.popfile | shuf | head -n 4 > tmp/tajima_${GROUP}_${IT}.popfile
        awk -v t=${n2} '$2 == t' $popdir/${GROUP}.popfile | shuf | head -n 6 >> tmp/tajima_${GROUP}_${IT}.popfile
else
        echo "NON-HUNGARIAN POPULATION, MAKING RUFOUS 10% GREY 90%"
        awk -v p=${n1} '$2 == p' $popdir/${GROUP}.popfile | shuf | head -n 9 | awk '{print $1}' > tmp/tajima_${GROUP}_${IT}.popfile
        awk -v t=${n2} '$2 == t' $popdir/${GROUP}.popfile | shuf | head -n 1 | awk '{print $1}' >> tmp/tajima_${GROUP}_${IT}.popfile
fi

# Also subset the full population
awk '{print $1}' $popdir/${GROUP}.popfile > tmp/tajima_${GROUP}_${IT}_FULL.popfile

# Cycle through all chromosomes
for CHR in $(cat Chromosomes.list); do

        # Cycle through different window sizes
        for win in 5000 50000 500000; do
        echo "WORKING ON ${CHR} FOR WIN SIZE ${win}"

        # If the chromosome is a sex chromosome, calculate haploid metrics
        if [[ ${CHR} == 'chr_W' || ${CHR} == 'chr_Z' ]]; then
                echo "SEX CHROMOSOME, CALCULATING HAPLOID METRICS"
                input=snp/${CHR}.IF-GF-MM2-BP-ANN-HAPMASK.vcf.gz
                # Tajima's D calculation for subset population
                ~/modules/vcftools/bin/vcftools --haploid --gzvcf ${input} --keep tmp/tajima_${GROUP}_${IT}.popfile --out divergence/work/${CHR}_${GROUP}_${IT}_${win}_SUB --TajimaD ${win}
                awk -v g=${GROUP} -v w=${win} -v i=${IT} '{OFS="\t"}{print $0, g, i, w, "Subset"}' divergence/work/${CHR}_${GROUP}_${IT}_${win}_SUB*.Tajima.D > divergence/out/${CHR}_${GROUP}_${IT}_${win}_SUB.tajima

                # Tajima's D calculation for full population
                ~/modules/vcftools/bin/vcftools --haploid --gzvcf ${input} --keep tmp/tajima_${GROUP}_${IT}_FULL.popfile --out divergence/work/${CHR}_${GROUP}_${IT}_${win}_ALL --TajimaD ${win}
                awk -v g=${GROUP} -v w=${win} -v i=${IT} '{OFS="\t"}{print $0, g, i, w, "All"}' divergence/work/${CHR}_${GROUP}_${IT}_${win}_ALL*.Tajima.D > divergence/out/${CHR}_${GROUP}_${IT}_${win}_ALL.tajima

        # If the chromosome is an autosome, calculate diploid metrics
        else
                echo "AUTOSOME, CALCULATING DIPLOID METRICS"
                input=snp/${CHR}.IF-GF-MM2-BP-ANN.vcf.gz
                # Tajima's D calculation for subset population
                ~/modules/vcftools/bin/vcftools --gzvcf ${input} --keep tmp/tajima_${GROUP}_${IT}.popfile --out divergence/work/${CHR}_${GROUP}_${IT}_${win}_SUB --TajimaD ${win}
                awk -v g=${GROUP} -v w=${win} -v i=${IT} '{OFS="\t"}{print $0, g, i, w, "Subset"}' divergence/work/${CHR}_${GROUP}_${IT}_${win}_SUB*.Tajima.D > divergence/out/${CHR}_${GROUP}_${IT}_${win}_SUB.tajima

                # Tajima's D calculation for full population
                ~/modules/vcftools/bin/vcftools --gzvcf ${input} --keep tmp/tajima_${GROUP}_${IT}_FULL.popfile --out divergence/work/${CHR}_${GROUP}_${IT}_${win}_ALL --TajimaD ${win}
                awk -v g=${GROUP} -v w=${win} -v i=${IT} '{OFS="\t"}{print $0, g, i, w, "All"}' divergence/work/${CHR}_${GROUP}_${IT}_${win}_ALL*.Tajima.D > divergence/out/${CHR}_${GROUP}_${IT}_${win}_ALL.tajima

        fi

        done
done


