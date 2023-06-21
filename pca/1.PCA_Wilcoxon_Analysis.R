# Change the current working directory to where the files are located
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/autosomes')

# Add a new library path
.libPaths('~/mambaforge/envs/r/lib/R/library')

# Load required libraries
library(tidyverse)
library(viridis)
library(GGally)
library(broom)
library(car)
library(ggpubr)

# Initialize variables
counter = 0
pdat=NULL
pcdat=NULL
RUNS= c('Autosomes.IF-GF-MM2-BP-ANN-AC2.AnalysisBIngroup','chr_W.IF-GF-MM2-BP-ANN-AC2-HAPMASK.AnalysisBIngroup')

# Loop through each run
for (RUN in RUNS) {
	counter = counter + 1

	# Load metadata
	mdk = read.table('../EntireMetadata.txt',header=TRUE,comment.char='',sep='\t')
	mdk$Species_Latin = factor(mdk$Species_Latin,levels=c('C. canorus','C. optatus','C. poliocephalus'))

	# Import eigenvectors
	eigs = read.table(paste0(RUN,'.eigenvec'),header=FALSE) %>% select(-V2)
	names(eigs) = c('ID',paste0('PC',seq(1,ncol(eigs)-1)))
	eigs = eigs[1:7]

	# Merge metadata and eigenvectors
	pcd = left_join(eigs,mdk)

	# Import eigenvalues
	relval = read.table(paste0(RUN,'.eigenval'),header=FALSE)
	relval = relval %>% mutate(VE = V1/sum(V1),
	                         label = paste0('PC',row_number(),' ',signif(VE,3)*100,'%'))

    # Plot the first 2 Principal Component (PC) axes, with different colors representing different plumage, and different shapes for different species
  	pcd %>% ggplot(aes(x=PC1,y=PC2,fill=Plumage,shape=Species_Latin))+
    geom_point(size=2)+
    scale_shape_manual(values=pcd$Shape,breaks=pcd$Species_Latin)+
    scale_fill_manual(values=pcd$PlumageColor,breaks=pcd$Plumage)+
    xlab(paste0('PC1: ',round(relval$VE[[1]],3)*100,'%'))+ylab(paste0('PC2: ',round(relval$VE[[2]],3)*100,'%'))+
    theme_classic()+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(legend.text = element_text(size = 8),legend.title = element_text(size = 8),legend.key.size = unit(0.3, 'cm'),legend.position="right")

	# Transform the data so that it can be easily visualized, with values normalized to fall between 0 and 1
	pcl = pcd %>% select(ID,contains(c('PC','Plumage','Species_Latin','Shape'))) %>% pivot_longer(!c(ID,Plumage,PlumageColor,Species_Latin,Shape)) %>%
    filter(!grepl('327',ID)) %>%
    group_by(name) %>%
    mutate(centered_value = value - min(value),
           scaled_value = centered_value / (max(value) - min(value)),
           axis = gsub('PC','',name)) %>% mutate_at('axis',as.numeric)

	# Plot the scaled PC scores by axis for each individual, with different colors for different plumage, and different shapes for different species
	lines = pcl %>% ggplot(aes(x=as.factor(axis),y=scaled_value,fill=Plumage,col=Plumage,group=ID,shape=Species_Latin))+
	geom_point(col='black',size=2)+
	geom_line(show.legend = FALSE)+
	scale_color_manual(values=pcl$PlumageColor,breaks=pcl$Plumage)+
	scale_shape_manual(values=pcl$Shape,breaks=pcl$Species_Latin)+
	scale_fill_manual(values=pcl$PlumageColor,breaks=pcl$Plumage)+
	scale_x_discrete(labels = relval[1:10,3]) +
	guides(fill=guide_legend(override.aes=list(shape=21)))+
	theme_classic()+xlab('')+ylab('Scaled PC Scores')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
	lines

	# Perform a Wilcoxon rank sum test to see if there are significant differences in the scaled PC scores across species and plumage
	p_values <- pcl %>%
	group_by(axis) %>%
	summarise(species = wilcox.test(scaled_value ~ Species_Latin)$p.value,
	          plumage = wilcox.test(scaled_value ~ Plumage)$p.value) %>%
	pivot_longer(!(axis))

	# Adjust the p-values for multiple testing using the Bonferroni correction
	p_values_adjusted <- p_values %>%
	mutate(padj = p.adjust(value, method = "bonferroni",n = 24),
	       significance = ifelse(padj < 0.05,'*','n.s.'),
	       chromosome = RUN)
  
	# Append the adjusted p-values to an existing data frame, `pdat`
	pdat = rbind(pdat,p_values_adjusted)

	# Plot the scaled Principal Component (PC) scores with significance annotations
	pp = pcl %>% ggplot(aes(x=as.factor(axis),y=scaled_value,fill=Plumage,shape=Species_Latin))+
	geom_point(size=2,position=position_dodge(width=0.5))+
	geom_text(data = p_values_adjusted,
	          aes(x = axis, y = 1, col = name, size=significance,label = significance),
	          inherit.aes = FALSE,vjust=0,position=position_dodge(width=0.5)) +
	scale_fill_manual(values=pcl$PlumageColor,breaks=pcl$Plumage)+xlab('')+ylab('Scaled PC Scores')+
	scale_shape_manual(values=pcl$Shape,breaks=pcl$Species_Latin)+
	scale_color_manual(values=c('grey20','grey60'))+
	scale_size_manual(values=c(6,4))+
	scale_x_discrete(labels = relval[1:6,3]) +
	theme_classic()

	# Assign the ggplot object to a variable dynamically named (e.g., p1, p2, etc.), based on the counter variable
	assign(paste0('p',counter),pp)

	# Append the reshaped data to an existing data frame, `pcdat`, adding an "Analysis" column
	pcdat = rbind(pcdat,pcl %>% mutate(Analysis = RUN))
}

# Save the arranged plots to a pdf file
pdf('Comparison_A-W_NOLD.pdf',height=4,width=6)
ggarrange(p1,p2,nrow=2,ncol=1,common.legend = TRUE)
dev.off()

# Save the prepared data for future use
pcdat = pcdat %>% mutate(Analysis = gsub('\\..*','',Analysis))
write.table(pcdat,file='PC_Scores_NO-LD.txt',quote=F,sep='\t',row.names=F)
pwdat = pdat %>% mutate(Analysis = gsub('\\..*','',chromosome)) %>% select(-c(chromosome,value,significance)) %>% ungroup %>% pivot_wider(names_from=name,values_from=padj)
write.table(pwdat,file='Pvalues_NO-LD.txt',quote=F,sep='\t',row.names=F)

# Function to check ANOVA assumptions: normality of residuals and homogeneity of variances
check_anova_assumptions <- function(data, axis) {
  # Filter data for the axis of interest
  data_axis <- data %>%
    filter(axis == axis)

  # Fit a two-way ANOVA model
  model <- aov(scaled_value ~ Plumage * Species_Latin, data = data_axis)

  # Check normality assumption using a QQ plot
  qq_plot <- ggplot(data_axis, aes(sample = scaled_value)) +
    geom_qq() +
    geom_qq_line() +
    ggtitle(paste("QQ plot for", axis))
  print(qq_plot)

  # Check homoscedasticity assumption using Levene's test
  levene_test <- car::leveneTest(scaled_value ~ Plumage * Species_Latin, data = data_axis)
  print(levene_test)
}

# Check ANOVA assumptions for each axis
unique_axes <- unique(data$axis)
for (axis in unique_axes) {
  cat("Checking assumptions for", axis, "\n")
  check_anova_assumptions(data, axis)
}

