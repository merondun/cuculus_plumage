#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=12:00:00

# Calculate ABBA-BABA for each chromosome and with canorus grey or canorus rufous as the target
CHR=$1
GROUP=$2

# Using Simon martins' scripts https://github.com/simonhmartin/genomics_general and the tutorial found here : https://github.com/simonhmartin/tutorials/tree/master/ABBA_BABA_whole_genome

~/modules/genomics_general/freq.py -t 5 -g ../input/${CHR}.geno.gz -p OG -p OH -p CG -p CH -p CP --popsFile ${GROUP}.pop --target derived -o ABBA/${CHR}-${GROUP}.derFreq.tsv.gz
