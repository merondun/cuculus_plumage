#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --time=200:00:00

GROUP=$1

popdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/randopops
genodir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/allsites
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/divergence


#randomize this first separately. the popfile simply shows $ID $POP, it will create 10 new files with with the first column shuffled with random pop ids 
awk '{print $2}' HH_HG.popfile > temp2

for i in {1..10}
do
    #shuffle lines in the file and pick first 10
    awk '{print $1}' HH_HG.popfile | shuf > temp.list

    #add the pop names to the shuffled files 
    paste temp.list temp2 > "rf${i}.popfile"
done

mkdir tmp divergence divergence/work divergence/out

#subsample to parity n = 5
popfile=$popdir/${GROUP}.popfile
n1=$(awk '{print $2}' $popdir/${GROUP}.popfile | sort | uniq | head -n 1)
n2=$(awk '{print $2}' $popdir/${GROUP}.popfile | sort | uniq | tail -n 1)

for CHR in $(cat RANDCHR.list); do
echo "WORKING ON CHR: ${GROUP} and ${CHR}"

if [[ ${CHR} == 'chr_W' || ${CHR} == 'chr_Z' ]]; then
        echo "SEX CHROMOSOME, CALCULATING HAPLOID METRICS"
        input=$genodir/${CHR}.HAPMASKgeno.gz
        ploidy=1
else
        echo "AUTOSOME, CALCULATING DIPLOID METRICS"
        input=$genodir/${CHR}.geno.gz
        ploidy=2
fi

        #cycle through windows
        for win in 50000; do
        winsize=${win}
        missize=$(echo ${winsize} | awk -v w=${winsize} '{print w/5}')

        echo "Populations: ${n1} and ${n2}, for ${CHR} and working on ${winsize} window with ${missize} missing sites max"

        #then calculate dxy/fst
        popgenWindows.py -w $winsize -m $missize -g ${input} \
                -o ${outdir}/work/${CHR}__${GROUP}_${win}_RAND.csv.gz \
                -f phased --windType coordinate --ploidy ${ploidy} -T 5 -p ${n1} -p ${n2} --popsFile $popfile
        zcat ${outdir}/work/${CHR}__${GROUP}_${win}_RAND.csv.gz | \
                awk -v g=${GROUP} '{OFS="\t"}{print $0, g}' | tr ',' '\t' > ${outdir}/out/${CHR}__${GROUP}_${win}_RAND.dxy

        done
done
