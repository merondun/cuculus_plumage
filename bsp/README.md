# Bayesian Skyline Plot Analyses

This dir involves generating Bayesian Skyline Plots (BSPs) for population genetic analysis, along with sensitivity assessments via subsampling random individuals.

## File Descriptions

1. `Extract_SNPs_Create_Fasta.sh`: Extracts SNPs, flanks them so we include invariant sites, and converts them into FASTA format for subsequent analyses.
2. `ML_Tree_To_Check.sh`: Constructs a Maximum Likelihood (ML) tree for initial data checks.
3. `Plot_ML_Tree.R`: Provides a visual representation of the ML tree.
4. `Random_Subsample_BSP.sh`: Implements the process of randomly subsampling individuals to assess sensitivity in the Bayesian Skyline Plot analyses.
5. `Full_Subsample_BSP.sh`: Implements subsampling across the entire dataset for a comprehensive assessment.
6. `RUN_BSP.sh`: Executes the main Bayesian Skyline Plot analyses.
7. `Plot_BSP.R`: Produces visual representations of the Bayesian Skyline Plots.

The project also includes several data files:

- `BSP_Individuals.list`: A list of individuals included in the BSP analysis.
- `chr_W.IF-GF-MM05-BP-ANN-HAP_BSP-All.10bp.bed`: Contains BED format data for 10bp flank regions around SNPs.
- `chr_W.IF-GF-MM05-BP-ANN-HAP_BSP-All.50bp.bed`: Contains BED format data for 50bp flanking regions around SNPs.
- `chr_W.IF-GF-MM05-BP-ANN-HAP_BSP-All.positions`: Contains position information for all SNPs included in the analysis.
- `chr_W.IF-GF-MM05-BP-ANN-HAP_BSP-All.vcf.gz`: The compressed VCF file that contains all SNP data.
- `chr_W.IF-GF-MM05-BP-ANN-HAP_BSP-All.vcf.gz.csi`: An index file for the VCF file.

