#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00

TREE=Retained_Genes_Between10k-75k_BSP.fa

#tree 
iqtree --redo -T 20 -s ${TREE} --seqtype DNA -m MFP -alrt 1000 -B 1000
