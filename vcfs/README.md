# chrW Data: VCF Files and Related Documentation

This repository holds VCF (Variant Call Format) files for the chrW chromosome data. The corresponding Maximum Likelihood (ML) tree was generated from these VCF files utilizing scripts located in the `/trees/` directory, which run the `vcf2phylip` and `iqtree` programs in sequence.

## File Structure

### 1. `chr_W.IF-GF-MM2-BP-ANN-HAPMASK__SetB.vcf.gz`

This VCF file holds chrW chromosome data for **Set B**, and its naming structure denotes specific data processing stages:

- **IF**: Info-Filtered
- **GF**: Genotype-Filtered
- **MM2**: Max-Missing rate of 20%
- **BP**: No Single Nucleotide Polymorphisms (SNPs) within 5 base pairs (BP) of one another
- **ANN**: Annotated via snpeff
- **HAPMASK**: Haploidized using `bcftools` and masked to remove sites overlapping the sex-chromosome mask, which identifies male to female coverage discrepancies
- **SetB**: Denotes Sample Set B, comprising higher quality individuals. 

### 2. `chr_W.IF-GF-MM2-BP-ANN-HAPMASK__SetC.vcf.gz`

Following the same structure as above, this file contains data for **Sample Set C**, which includes the complete dataset with toepad samples.

### 3. `chr_W.IF-GF-MM2-BP-ANN-HAPMASK__SetB.positions.bed`

This bed file signifies the positions covered in the VCF file. The 1-based VCF coordinates were converted to 0-based bed coordinates using:

`bcftools query -f "%CHROM\t%POS\t%POS\n" chr_W.IF-GF-MM2-BP-ANN-HAPMASK__SetB.vcf.gz | awk '{OFS="\t"}{print $1, $2-1, $3}' > chr_W.IF-GF-MM2-BP-ANN-HAPMASK__SetB.positions.bed`

