# ML-Trees

This project conducts phylogenetic analyses using [IQ-TREE](http://www.iqtree.org/) to create maximum likelihood (ML) trees with SH-aLRT and bootstrap support values for nodes. Scripts are included for visualizing the ML trees and for co-plotting mitochondrial DNA and W chromosome data.

## Directory Structure

The directories are organized by different sets of genetic data, each allowing a different percentage of missing genotypes:

- `max-missing-05_SetB`: Contains data for Set B, with a maximum allowed missing genotype rate of 5%.
- `max-missing-05_SetC`: Contains data for Set C, with a maximum allowed missing genotype rate of 5%.
- `max-missing-20_SetC`: Contains data for Set C, allowing for a higher missing genotype rate of up to 20%.
- `max-missing-20_SetB`: Contains data for Set B, also allowing a maximum missing genotype rate of 20%.

## File Descriptions

- `Create_ML_Tree.sh`: A shell script that utilizes IQ-TREE to generate Maximum Likelihood (ML) trees from VCFs.

- `Plot_MLtrees_and_Cophyloplot.R`: An R script that uses the resulting ML trees from the previous script and visualizes them. It also creates a cophylogenetic plot of mitochondrial DNA and W chromosome data.

