#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=5
#SBATCH --time=200:00:00

# Define variables with the paths to different directories and files
vcfdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/2022_11/raw_vcf
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
ints=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/2022_11/snpcalling/intervals

# Assign command line arguments to variables
CHR=$1  # This is the genomic interval ID
GROUP=$2  # This is a text file containing a list of samples

# Create a directory to store the variant call format (vcf) files
mkdir $vcfdir

# Set the max depth parameter based on the chromosome value
if [[ $CHR == 0102 ]]
then
        maxdp=20000  # Set a higher max depth for chromosome 102
else
        maxdp=100  # Use a default max depth for other chromosomes
fi

# Inform the user about the interval and max depth being used
echo "WORKING ON INTERVAL ${CHR} WITH MAX DEPTH ${maxdp}"

# Run bcftools mpileup to create a pileup of bases at each position in the interval,
# and then call variants from the pileup using bcftools call.
# Save the results in a gzipped vcf file in the vcf directory.
bcftools mpileup --max-depth ${maxdp} -C 50 -threads 5 -f ${genome} -b ${GROUP}.bam -R $ints/$CHR.bed -a "AD,ADF,ADR,DP,SP" | \
        bcftools call --ploidy 2 --threads 5 -m -Oz -o ${vcfdir}/${CHR}_bcft.vcf.gz


