#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=15
#SBATCH --time=200:00:00

# This script accepts two parameters: Group name and iteration number
GROUP=$1
IT=$2

# Define the directories for population data, genotype data, and output
popdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/populations
genodir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/allsites
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/divergence

# Create necessary directories
mkdir tmp divergence divergence/work divergence/out

# Popfile contains $ID $POPULATION
popfile=$popdir/${GROUP}.popfile
# Extract distinct populations from popfile
n1=$(awk '{print $2}' $popdir/${GROUP}.popfile | sort | uniq | head -n 1)
n2=$(awk '{print $2}' $popdir/${GROUP}.popfile | sort | uniq | tail -n 1)

# Randomly subsample 5 individuals from each population for the specified iteration
awk -v p=${n1} '$2 == p' $popdir/${GROUP}.popfile | shuf | head -n 5 > tmp/${GROUP}_${IT}.popfile
awk -v t=${n2} '$2 == t' $popdir/${GROUP}.popfile | shuf | head -n 5 >> tmp/${GROUP}_${IT}.popfile

# Loop over all chromosomes
for CHR in $(cat Chromosomes.list); do
echo "WORKING ON CHR: ${GROUP} and ${CHR}"

# If the chromosome is a sex chromosome, calculate haploid metrics
if [[ ${CHR} == 'chr_W' || ${CHR} == 'chr_Z' ]]; then
        echo "SEX CHROMOSOME, CALCULATING HAPLOID METRICS"
        input=$genodir/${CHR}.HAPMASKgeno.gz
        ploidy=1
# If the chromosome is an autosome, calculate diploid metrics
else
        echo "AUTOSOME, CALCULATING DIPLOID METRICS"
        input=$genodir/${CHR}.geno.gz
        ploidy=2
fi

        # Loop through windows (here only one window size is defined as 50000 but more could be added)
        for win in 50000; do
        winsize=${win}
        # Maximum number of missing sites per window is 1/5 of the window size
        missize=$(echo ${winsize} | awk -v w=${winsize} '{print w/5}')

        # Calculate dxy and fst for each window using popgenWindows.py and output to a csv file
        popgenWindows.py -w $winsize -m $missize -g ${input} \
                -o ${outdir}/work/${CHR}__${GROUP}_${IT}_${win}_SUBM.csv.gz \
                -f phased --windType coordinate --ploidy ${ploidy} -T 15 -p ${n1} -p ${n2} --popsFile tmp/${GROUP}_${IT}.popfile

        # Convert the output csv to a tab-delimited file and add some extra columns
        zcat ${outdir}/work/${CHR}__${GROUP}_${IT}_${win}_SUBM.csv.gz | \
                awk -v g=${GROUP} -v i=${IT} -v w=${winsize} '{OFS="\t"}{print $0, g, i, w, "SubsetMask"}' | tr ',' '\t' > ${outdir}/out/${CHR}__${GROUP}_${IT}_${win}_SUBM.dxy

        # Repeat the process for the full sample set for comparison
        popgenWindows.py -w $winsize -m $missize -g ${input} \
                -o ${outdir}/work/${CHR}__${GROUP}_${IT}_${win}_ALLM.csv.gz \
                -f phased --windType coordinate --ploidy ${ploidy} -T 15 -p ${n1} -p ${n2} --popsFile $popdir/${GROUP}.popfile
        zcat ${outdir}/work/${CHR}__${GROUP}_${IT}_${win}_ALLM.csv.gz | \
                awk -v g=${GROUP} -v i=${IT} -v w=${winsize} '{OFS="\t"}{print $0, g, i, w, "AllMask"}' | tr ',' '\t' > ${outdir}/out/${CHR}__${GROUP}_${IT}_${win}_ALLM.dxy

        done
done

