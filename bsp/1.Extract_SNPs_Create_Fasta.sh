#Grab a bed file with gene coordinates from the gff annotation 
awk '$1 == "chr_W" && $3 == "gene"' GCA_017976375.1_bCucCan1.pri_genomic.CHR.gff | grep 'protein_coding' | awk '{OFS="\t"}{print $1, $4, $5, $7, $9}' | sed 's/ID=gene-//g' | sed 's/;.*//g' | awk '{OFS="\t"}{print $1, $2, $3, $5, ".", $4}'  > ~/symlinks/toes/bsp/2023JUNE_GENES/Gene_Coordinates.bed


### Extract fastas for each gene for each sample 

#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00

#Run this in an empty directory. Specify the paths to the standard files and it will output a file for each gene (indicated by name in field 4 of the bed annotation), with the sample names from the vcf as the sequence header. This is primarily useful for haploidized VCFs e.g. chrW or mtDNA 

fasta=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/chr_W/GCA_017976375.1_bCucCan1.pri_chr_W.fa
annotation=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/bsp/2023JUNE_GENES/Gene_Coordinates.bed
samples=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/bsp/2023JUNE_GENES/BSP_Individuals.list
vcf=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/allsites/chr_W.AllSitesHAPMASK.vcf.gz

rm *fa
#loop through genes
for gene in $(awk '{print $4}' ${annotation}); do

	echo ${gene}
	coords=$(grep -w ${gene} ${annotation} | awk '{print $1":",$2,"-",$3"}')

        #loop through samples
        for sample in $(cat ${samples}); do
        echo ${sample}
        samtools faidx ${fasta} ${coords} | bcftools consensus --haplotype A ${vcf} --sample ${sample} | sed "1 s/^>.*/>${sample}/" >> ${gene}.fa

        done
done
