#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

RUN=$1
K=$2
echo "${RUN} at ${K}"
mkdir admixture

cd admixture
admixture -j10 --cv=5 ../autosomes/${RUN}.bed ${K} > ${RUN}.log${K}.out
