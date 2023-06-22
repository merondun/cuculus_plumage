#extract SNPs with less than 5% missing genotypes for the individuals we will use for BSP 
bcftools view --samples-file bsp/BSP_Individuals.list -Ou vcfs/chr_W.IF-GF-MM2-BP-ANN-HAPMASK__SetB.vcf.gz | \
        bcftools view -e 'F_MISSING>0.2' -v snps -m2 -M2 --min-ac 1 -Oz -o bsp/chr_W.IF-GF-MM2-BP-SNPs.vcf.gz

bcftools query -f "%CHROM\t%POS\t%POS\n" bsp/chr_W.IF-GF-MM2-BP-SNPs.vcf.gz| awk '{OFS="\t"}{print $1, $2-1, $3}' > bsp/chr_W.IF-GF-MM2-BP-ANN-HAPMASK__SetB.positions.bed

#extend each SNP by 50bp in each direction (100bp fragments), and merge the resultant beds 
bedtools slop -i bsp/chr_W.IF-GF-MM2-BP-ANN-HAPMASK__SetB.positions.bed -b 50 -g general/Cuckoo_Genome_Chromosome_Lengths.txt | bedtools merge -i - > bsp/chr_W.IF-GF-MM2-BP-ANN-HAPMASK.100bp.bed

#extract invariant sites
bedtools intersect -header -a /dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/allsites/chr_W.AllSitesHAPMASK.vcf.gz -b ~/merondun/cuculus_plumage/bsp/chr_W.IF-GF-MM2-BP-ANN-HAPMASK.100bp.bed |
        bcftools view --samples-file ~/merondun/cuculus_plumage/bsp/BSP_Individuals.list --force-samples -Oz -o bsp/chr_W.BSP.100bp.vcf.gz

#index
for i in $(ls bsp/*100bp.vcf.gz); do bcftools index -f ${i}; done
for i in $(ls bsp/*100bp.vcf.gz); do bcftools index -n ${i}; done #this will count SNP, ensure that the files merged correctly

#Tree, first convert vcf to phylip
python ~/modules/vcf2phylip.py -i bsp/chr_W.BSP.100bp.vcf.gz --output-folder bsp/

#then use awk to convert phylip to fasta
awk 'NR>1 {print ">"$1"\n"$2}' bsp/chr_W.BSP.100bp.min4.phy > bsp/chr_W.BSP.100bp.fa

#Randomly subsample n = 5 samples
for j in $(cat ../Groups.list); do 
cat ${j}.pop | shuf -n 5 > sub.tmp
seqtk subseq ../chr_W.BSP.100bp.fa sub.tmp | seqret -sequence stdin -outseq 100_mac1_sub/${j}.nex -osformat nexus
done

#Grab all samples
for j in $(cat ../Groups.list); do 
seqtk subseq ../chr_W.BSP.100bp.fa ${j}.pop | seqret -sequence stdin -outseq 100_mac1_full/${j}.nex -osformat nexus
done
