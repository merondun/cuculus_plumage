#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00

#tree 
iqtree --redo -T 20 -s chr_W.BSP.AllSites.50bp.fa --seqtype DNA -m MFP -alrt 1000 -B 1000
iqtree --redo -T 20 -s chr_W.BSP.AllSites.10bp.fa --seqtype DNA -m MFP -alrt 1000 -B 1000
