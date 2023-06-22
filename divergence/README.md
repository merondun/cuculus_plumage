# Estimating Population Genetic Statistics and Visualization

This dir involves the estimation of population genetic statistics such as DXY, Fst, and Tajima's D. The statistics are estimated for specific genomic windows and then visualized. 

## File Descriptions

1. `Estimate_DXY_Fst_Windows.sh`: This script estimates DXY and Fst statistics for defined genomic windows.
2. `Estimate_DXY_Fst_Windows-Randomized.sh`: Similar to the above script but includes an additional step of randomizing the data to see if chrW shows the pattern for random samples. 
3. `Estimate_TajimasD.sh`: Estimates Tajima's D statistic for genomic windows. 
4. `Plot_DXY_Fst_Heatmaps_Manhattan.R`: This R script visualizes the computed DXY and Fst values in the form of heatmaps and Manhattan plots.
5. `Plot_Ranomdized_FST.R`: Plots the Fst values obtained from the randomized data to examine the effect of potential bias.
6. `Plot_Genotypes.R`: An additional script that generates visual representations of genotype data from a VCF.

Data files include:

- `DXY_FST_Input_50KB-MaxMissing20.txt`: Input file for DXY and Fst estimations, representing 50KB genomic windows and allowing a maximum of 20% missing data.
- `Randomized_DXY_FST_Input_50KB-MaxMissing20.txt`: The randomized version of the above file, used to assess potential biases in the data.


