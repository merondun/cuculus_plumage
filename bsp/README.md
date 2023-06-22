# Bayesian Skyline Plot Analyses

This dir involves generating Bayesian Skyline Plots (BSPs) for population genetic analysis, along with sensitivity assessments via subsampling random individuals.

## File Descriptions

1. `Extract_SNPs_Create_Fasta.sh`: Extracts SNPs, flanks them so we include invariant sites, and converts them into FASTA format for subsequent analyses.
2. `ML_Tree_To_Check.sh`: Constructs a Maximum Likelihood (ML) tree for initial data checks.
3. `Plot_ML_Tree.R`: Provides a visual representation of the ML tree.
6. `RUN_BSP.sh`: Executes the main Bayesian Skyline Plot analyses.
7. `Plot_BSP.R`: Produces visual representations of the Bayesian Skyline Plots.

The project also includes several data files:

- `BSP_Individuals.list`: A list of individuals included in the BSP analysis.
- `chr_W.IF-GF-MM2-BP-SNPs.vcf.gz`: VCF file without invariant sites filtered for the individuals, this was then flanked by 50bp each side to create the fasta.
- `chr_W.BSP.100bp.fa.gz`: MSA Fasta file created from the above scripts and used for BSP analysis.

