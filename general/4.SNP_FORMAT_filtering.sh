#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --mem=20000mb
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00


#PLUMAGE
genome=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa
#filter snp according to missingness genotypes, and prune any SNPs that occur within 5bp of one another. Annotate with snpeff, output some stats, and also create lostruct input. Also output invariant sites for dxy
cd /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes

mkdir vcfqc snp allsites phased lostruct lostruct/input

CHR=$1
echo "${CHR} PLUMAGE"


#filter for missingness
bcftools view --samples-file Retain.list --force-samples -Ou snp/${CHR}.IF-GF.vcf.gz | bcftools view --min-ac 1 --threads 4 -i 'F_MISSING<0.2' -Oz -o snp/${CHR}.IF-GF-MM2.vcf.gz
bcftools index --threads 4 snp/${CHR}.IF-GF-MM2.vcf.gz
#filter for SNPs within 5bp for MT, 50bp for autosomes [KEEPS THE ALLELE WITH THE HIGHEST MAF, DEFAULT]
if [[ $CHR = 'chr_MT' ]]
then
        echo "PRUNING SNPS WITHIN 5BP"
        bcftools +prune -Oz -o snp/${CHR}.IF-GF-MM2-BP.vcf.gz --window 5bp --nsites-per-win 1 snp/${CHR}.IF-GF-MM2.vcf.gz
        bcftools index --threads 4 snp/${CHR}.IF-GF-MM2-BP.vcf.gz
else
        echo "PRUNING SNPS WITHIN 5BP"
        bcftools +prune -Oz -o snp/${CHR}.IF-GF-MM2-BP.vcf.gz --window 5bp --nsites-per-win 1 snp/${CHR}.IF-GF-MM2.vcf.gz
        bcftools index --threads 4 snp/${CHR}.IF-GF-MM2-BP.vcf.gz
fi

#annotate
java -Xmx25g -jar ~/modules/snpEff/snpEff.jar -v -stats vcfqc/${CHR}_annotation.html bCucCan1 snp/${CHR}.IF-GF-MM2-BP.vcf.gz | bgzip -c > snp/${CHR}.IF-GF-MM2-BP-ANN.vcf.gz

#remove singletons
bcftools view --threads 4 -U -i 'TYPE=="snp" & MAC >= 2' -Oz -o snp/${CHR}.IF-GF-MM2-BP-ANN-AC2.vcf.gz snp/${CHR}.IF-GF-MM2-BP-ANN.vcf.gz
bcftools index --threads 4 snp/${CHR}.IF-GF-MM2-BP-ANN-AC2.vcf.gz

#summarize
sites=$(bcftools index -n snp/${CHR}.vcf.gz)
snps=$(bcftools view --min-alleles 2 --max-alleles 2 snp/${CHR}.vcf.gz | grep -v '#' | wc -l)
infofilt=$(zcat snp/${CHR}.IF.vcf.gz | grep -v '#' | wc -l)
genofilt=$(zcat snp/${CHR}.IF-GF-MM2.vcf.gz | grep -v '#' | wc -l)
bpfilt=$(zcat snp/${CHR}.IF-GF-MM2-BP.vcf.gz | grep -v '#' | wc -l)
acfilt=$(zcat snp/${CHR}.IF-GF-MM2-BP-ANN-AC2.vcf.gz | grep -v '#' | wc -l)
echo -e "CHR\tsites\tsnps\tinfofilter\tgenofilter\tbpfilter\tacfilter" > vcfqc/${CHR}.nsites.txt
echo "${CHR} ${sites} ${snps} ${infofilt} ${genofilt} ${bpfilt} ${acfilt}" >> vcfqc/${CHR}.nsites.txt

#separate invariant sites
echo "FILTERING AND MERGING INVARIANT SITES FOR ${CHR}"
bcftools view --samples-file Retain.list --force-samples -Ou snp/${CHR}.vcf.gz | bcftools view --max-ac 0 --threads 4 -Oz -o allsites/${CHR}.1N.vcf.gz
bcftools index --threads 4 allsites/${CHR}.1N.vcf.gz

#re-merge the invariant and filtered SNPs
bcftools concat --threads 4 -Ou allsites/${CHR}.1N.vcf.gz snp/${CHR}.IF-GF-MM2-BP-ANN-AC2.vcf.gz | bcftools sort -Oz -o allsites/${CHR}.AllSites.vcf.gz
bcftools index --threads 4 allsites/${CHR}.AllSites.vcf.gz

#create simon's divergence input file from the filtered vcf and the raw vcf
echo "CREATE SIMON MARTIN DIVERGENCE OUTPUT "
python ~/modules/genomics_general/VCF_processing/parseVCF.py --ploidyMismatchToMissing --ploidy 2 --skipIndels -i allsites/${CHR}.AllSites.vcf.gz | bgzip > allsites/${CHR}.geno.gz

# subset samples
bcftools view --threads 10 --force-samples --samples-file INGROUP.list snp/${CHR}.IF-GF-MM2-BP-ANN-AC2.vcf.gz -Ou | \
        bcftools view --min-alleles 2 --max-alleles 2 --min-ac 2 -Oz -o phased/${CHR}.IF-GF-MM2-BP-ANN-AC2-IG.vcf.gz
bcftools index phased/${CHR}.IF-GF-MM2-BP-ANN-AC2-IG.vcf.gz

# phase vcf with beagle
java -jar ~/modules/beagle.28Jun21.220.jar gt=phased/${CHR}.IF-GF-MM2-BP-ANN-AC2-IG.vcf.gz out=phased/${CHR}.IF-GF-MM2-BP-ANN-AC2-IG.phased nthreads=5 window=40 overlap=2 impute=true
bcftools index --threads 5 phased/${CHR}.IF-GF-MM2-BP-ANN-AC2-IG.phased.vcf.gz

# prep lostruct
bcftools view -Ou -o - phased/${CHR}.IF-GF-MM2-BP-ANN-AC2-IG.phased.vcf.gz | \
    bcftools convert -Ob  > lostruct/input/${CHR}.bcf
bcftools index lostruct/input/${CHR}.bcf
