# PCA and Wilcoxon Analysis Results

This directory contains results from Principal Component Analysis (PCA) and Wilcoxon Analysis on chromosome data. Here's a brief explanation of the files:

## File Descriptions

### 1. `PCA_Wilcoxon_Analysis.R`

This is the R script used to perform PCA and Wilcoxon Analysis on the chromosome data.

### 2. `Autosomes.IF-GF-MM2-BP-ANN-AC2.AnalysisBIngroup.eigenvec` and `chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK.AnalysisBIngroup.eigenvec`

These are the eigenvectors resulting from PCA on autosomes and chrW chromosome data respectively. The naming convention indicates the filters and operations applied to the data:
- **IF**: Info-Filtered
- **GF**: Genotype-Filtered
- **MM2**: Max-Missing rate of 20%
- **BP**: No SNPS within 5 base pairs of one another
- **ANN**: Annotated via snpeff
- **AC2**: Minor allele cont >= 2
- **HAPMASK**: Haploidized using `bcftools` and masked to remove sites overlapping the sex-chromosome mask, which identifies male to female coverage discrepancies
- **AnalysisBIngroup**: Denotes analysis was performed on sample Set B.

### 3. `Autosomes.IF-GF-MM2-BP-ANN-AC2.AnalysisBIngroup.eigenval` and `chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK.AnalysisBIngroup.eigenval`

These files contain eigenvalues corresponding to the eigenvectors from the PCA. The naming convention is similar to the eigenvectors files.


