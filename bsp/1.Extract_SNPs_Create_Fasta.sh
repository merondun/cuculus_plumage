#extract SNPs with less than 5% missing genotypes for the individuals we will use for BSP 
bcftools view --samples-file bsp/BSP_Individuals.list --force-samples -Ou vcfs/chr_W.IF-GF-MM2-BP-ANN-HAPMASK__SetB.vcf.gz | 
        bcftools view --min-alleles 2 --max-alleles 2 -i 'F_MISSING<0.05' -Oz -o bsp/chr_W.IF-GF-MM05-BP-ANN-HAP_BSP-All.vcf.gz
#extend each SNP by 25bp in each direction (50bp fragments), and merge the resultant beds 
bedtools slop -i bsp/chr_W.IF-GF-MM05-BP-ANN-HAP_BSP-All.positions -b 25 -g general/Cuckoo_Genome_Chromosome_Lengths.txt | bedtools merge -i - > bsp/chr_W.IF-GF-MM05-BP-ANN-HAP_BSP-All.50bp.bed

#for sensitivity also check 10 bp
bedtools slop -i bsp/chr_W.IF-GF-MM05-BP-ANN-HAP_BSP-All.positions -b 5 -g general/Cuckoo_Genome_Chromosome_Lengths.txt | bedtools merge -i - > bsp/chr_W.IF-GF-MM05-BP-ANN-HAP_BSP-All.10bp.bed

#intersect with the full SNP file, before we do any clustered BP-pruning or missingness 
bedtools intersect -header -a snp/chr_W.IF.ANN-HAP.vcf.gz -b ~/merondun/cuculus_plumage/bsp/chr_W.IF-GF-MM05-BP-ANN-HAP_BSP-All.50bp.bed |
        bcftools view --samples-file bsp/BSP.keep --force-samples -Oz -o bsp/chr_W.BSP.2N.50bp.vcf.gz
#and 10bp
bedtools intersect -header -a snp/chr_W.IF.ANN-HAP.vcf.gz -b ~/merondun/cuculus_plumage/bsp/chr_W.IF-GF-MM05-BP-ANN-HAP_BSP-All.10bp.bed |
        bcftools view --samples-file bsp/BSP.keep --force-samples -Oz -o bsp/chr_W.BSP.2N.10bp.vcf.gz

#extract invariant sites
bedtools intersect -header -a allsites/chr_W.1N-HAP.vcf.gz -b ~/merondun/cuculus_plumage/bsp/chr_W.IF-GF-MM05-BP-ANN-HAP_BSP-All.50bp.bed |
        bcftools view --samples-file bsp/BSP.keep --force-samples -Oz -o bsp/chr_W.BSP.1N.50bp.vcf.gz
#and 10bp
bedtools intersect -header -a allsites/chr_W.1N-HAP.vcf.gz -b ~/merondun/cuculus_plumage/bsp/chr_W.IF-GF-MM05-BP-ANN-HAP_BSP-All.10bp.bed |
    bcftools view --samples-file bsp/BSP.keep --force-samples -Oz -o bsp/chr_W.BSP.1N.10bp.vcf.gz

#merge VCFs 
bcftools concat bsp/chr_W.BSP.2N.50bp.vcf.gz bsp/chr_W.BSP.1N.50bp.vcf.gz -Ou | 
        bcftools sort -Oz -o bsp/chr_W.BSP.AllSites.50bp.vcf.gz
bcftools concat bsp/chr_W.BSP.2N.10bp.vcf.gz bsp/chr_W.BSP.1N.10bp.vcf.gz -Ou | 
        bcftools sort -Oz -o bsp/chr_W.BSP.AllSites.10bp.vcf.gz
#index
for i in $(ls bsp/*vcf.gz); do bcftools index -f ${i}; done
for i in $(ls bsp/*vcf.gz); do bcftools index -n ${i}; done #this will count SNP, ensure that the files merged correctly

#Tree, first convert vcf to phylip
python ~/modules/vcf2phylip.py -i bsp/chr_W.BSP.AllSites.50bp.vcf.gz --output-folder bsp/
python ~/modules/vcf2phylip.py -i bsp/chr_W.BSP.AllSites.10bp.vcf.gz --output-folder bsp/

#then use awk to convert phylip to fasta
awk 'NR>1 {print ">"$1"\n"$2}' bsp/chr_W.BSP.AllSites.50bp.min4.phy > bsp/chr_W.BSP.AllSites.50bp.fa
awk 'NR>1 {print ">"$1"\n"$2}' bsp/chr_W.BSP.AllSites.10bp.min4.phy > bsp/chr_W.BSP.AllSites.10bp.fa
