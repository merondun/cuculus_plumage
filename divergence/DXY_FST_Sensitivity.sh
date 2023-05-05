#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --time=200:00:00

GROUP=$1
IT=$2

popdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/populations
genodir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/allsites
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/divergence

mkdir tmp divergence divergence/work divergence/out

#subsample to parity n = 5
popfile=$popdir/${GROUP}.popfile
n1=$(awk '{print $2}' $popdir/${GROUP}.popfile | sort | uniq | head -n 1)
n2=$(awk '{print $2}' $popdir/${GROUP}.popfile | sort | uniq | tail -n 1)
awk -v p=${n1} '$2 == p' $popdir/${GROUP}.popfile | shuf | head -n 5 > tmp/${GROUP}_${IT}.popfile
awk -v t=${n2} '$2 == t' $popdir/${GROUP}.popfile | shuf | head -n 5 >> tmp/${GROUP}_${IT}.popfile

for CHR in $(cat Chromosomes.list); do
echo "WORKING ON CHR: ${GROUP} and ${CHR}"

if [[ ${CHR} == 'chr_W' || ${CHR} == 'chr_Z' ]]; then
        echo "SEX CHROMOSOME, CALCULATING HAPLOID METRICS"
        input=$genodir/${CHR}.HAPgeno.gz
        ploidy=1
else
        echo "AUTOSOME, CALCULATING DIPLOID METRICS"
        input=$genodir/${CHR}.geno.gz
        ploidy=2
fi 

        #cycle through windows
        for win in 5000 50000 500000; do
        winsize=${win}
        missize=$(echo ${winsize} | awk -v w=${winsize} '{print w/5}')

        echo "Populations: ${n1} and ${n2}, for ${CHR} and iteration ${IT} working on ${winsize} window with ${missize} missing sites max"

        #then calculate dxy/fst
        popgenWindows.py -w $winsize -m $missize -g ${input} \
                -o ${outdir}/work/${CHR}__${GROUP}_${IT}_${win}_SUB.csv.gz \
                -f phased --windType coordinate --ploidy ${ploidy} -T 5 -p ${n1} -p ${n2} --popsFile tmp/${GROUP}_${IT}.popfile
        zcat ${outdir}/work/${CHR}__${GROUP}_${IT}_${win}_SUB.csv.gz | \
                awk -v g=${GROUP} -v i=${IT} -v w=${winsize} '{OFS="\t"}{print $0, g, i, w, "Subset"}' | tr ',' '\t' > ${outdir}/out/${CHR}__${GROUP}_${IT}_${win}_SUB.dxy

        #also perform on the full sample set for comparison 
        popgenWindows.py -w $winsize -m $missize -g ${input} \
                -o ${outdir}/work/${CHR}__${GROUP}_${IT}_${win}_ALL.csv.gz \
                -f phased --windType coordinate --ploidy ${ploidy} -T 5 -p ${n1} -p ${n2} --popsFile $popdir/${GROUP}.popfile
        zcat ${outdir}/work/${CHR}__${GROUP}_${IT}_${win}_ALL.csv.gz | \
                awk -v g=${GROUP} -v i=${IT} -v w=${winsize} '{OFS="\t"}{print $0, g, i, w, "All"}' | tr ',' '\t' > ${outdir}/out/${CHR}__${GROUP}_${IT}_${win}_ALL.dxy

        done 
done 
