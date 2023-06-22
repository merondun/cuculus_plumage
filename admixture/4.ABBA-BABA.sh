#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=10000mb
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

# Calculate ABBA-BABA for each chromosome and with canorus grey or canorus rufous as the target
CHR=$1
GROUP=$2

# Using Simon martins' scripts https://github.com/simonhmartin/genomics_general and the tutorial found here : https://github.com/simonhmartin/tutorials/tree/master/ABBA_BABA_whole_genome
~/modules/genomics_general/freq.py -t 5 -g allsites/${CHR}.geno.gz -p OG -p OH -p ${GROUP} -p CP --popsFile populations/${GROUP}_ABBA.pop --target derived -o ABBA/${CHR}-${GROUP}.derFreq.tsv.gz

