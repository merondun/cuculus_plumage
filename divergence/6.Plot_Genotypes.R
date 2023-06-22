library(vcfR)
library(tidyverse)
library(viridis)


# Specify the VCF file location
vcf_file <- "~/merondun/cuculus_plumage/divergence/Fixed_Grey-Rufous.vcf"

# Load VCF data into an object
vcf <- read.vcfR(vcf_file)

# Extract the genotypes from the VCF file
gt_data <- vcfR2tidy(vcf, format_fields = 'GT' )
gt = gt_data$gt

# Rename cols
names(gt) = c('Key','site','ID','Genotype','Allele')

# Load metadata, as we want to group by morph
md = read.table('~/merondun/cuculus_plumage/SimpleMetadata.txt',header=TRUE,comment.char = '',sep='\t')

# Join the metadata and genotype data
gt = left_join(md,gt) #ensure that the 'ID' from the vcf matches the 'ID' from your metadata

# Convert the 'Group' column to a factor ans specicyf the levels the you want
gt$Group = factor(gt$Group,levels=c('CG','OG','CH','OH'))

# Create a plot of genotypes across individuals and morphs
gtp = gt %>% filter(Allele != '.' &
                      site < 2000000) %>%  # to reduce the number of SNPs for visualization..
  ggplot(aes(x = IDtree, y = as.factor(site), fill = Allele)) +
  geom_tile() +
  xlab("Individuals") +
  facet_grid(Group~.,space='free',scales='free')+
  scale_fill_viridis("Genotype", discrete=TRUE) +
  theme_minimal(base_size=10) + ylab('chrW: SNP')+xlab('') +
  theme(axis.text.x = element_text(size=4,angle = 90, vjust = 0.5, hjust = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        strip.text.y.right = element_text(angle = 0))+
  coord_flip()
gtp
