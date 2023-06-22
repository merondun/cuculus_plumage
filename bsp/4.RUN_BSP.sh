#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=12:00:00

#Create XML in beauti:
#
#* Gamma site model, 4 categories, GTR w/ estimated frequencies
#* Strict clock with 0.00000000454 as a rate 
#* Coalescent bayesian skyline model 30000000 chains, log every 3k 

GROUP=$1
mkdir ${GROUP}
cd ${GROUP}
beast -threads 3 -overwrite -java ../${GROUP}.xml
