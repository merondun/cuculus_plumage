# Bioinformatics Pipeline

This repository contains scripts and files used for general processing of the data before analysis. The scripts should be run in the order of their numbered prefixes.

## Scripts

- `0.a.Sample_summaries_spatial_plotting.R`: This script is used for creating spatial distribution plots for samples.
- `0.Trimming_Mapping.sh`: This script performs initial trimming and mapping of sequence reads.
- `1.SNP_Calling.sh`: This script calls SNPs from mapped reads.
- `2.Chromosome_Merging.sh`: This script merges chromosome-wide SNPs.
- `3.SNP_INFO_filtering.sh`: This script applies filtering to the SNPs based on information from INFO field in the VCF.
- `4.SNP_FORMAT_filtering.sh`: This script applies further filtering to the SNPs based on genotype format fields in the VCF.
- `5.Merge_Autosomes.sh`: This script is used to merge the SNP calls from individual autosomes.
- `6.Plink_PCA-PrepAdmixture.sh`: This script prepares data for PCA and Admixture analysis using plink.
- `7.SNP_Summaries.R`: This script generates summary statistics for SNPs.

## Other Files

- `Cuckoo_Genome_Chromosome_Lengths.txt`: This file contains lengths of chromosomes in the cuckoo genome.
- `SNP_Numbers.txt`: This text file contains the numbers of SNPs per chromosome.
- `SNP_Stats.txt`: This text file contains various statistics for the SNPs.

## Spatial Distribution Plots

- `Spatial-Distribution_chr_MT_SETC_1DEGREEJITTER.pdf`
- `Spatial-Distribution_chr_W_SETC_1DEGREEJITTER.pdf`
- `Spatial-Distribution_chr_W_SETC-MM05_1DEGREEJITTER.pdf`

These PDFs contain plots showing the spatial distribution of samples, with a degree of jitter added to prevent overplotting.


