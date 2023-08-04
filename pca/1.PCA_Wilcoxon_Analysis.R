#### PCA from PLINK, analyze significance with a Wilcoxon test since it grossly violates normality 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/autosomes')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(GGally)
library(broom)
library(car)
library(ggpubr)

### W and Autosome PC1 side-by-side

#Read metadata
mdk = read.table('../EntireMetadata.txt',header=TRUE,comment.char='',sep='\t')
mdk$Species_Latin = factor(mdk$Species_Latin,levels=c('C. canorus','C. optatus','C. poliocephalus'))

#import eigenvectors AUTOSOMES
eigs = read.table('Autosomes.IF-GF-MM2-BP-ANN-AC2.TransSpecies.eigenvec',header=FALSE) %>% select(-V2)
names(eigs) = c('ID',paste0('PC',seq(1,ncol(eigs)-1)))
eigs = eigs[1:2]

#import eigenvalues AUTOSOMES
relval = read.table('Autosomes.IF-GF-MM2-BP-ANN-AC2.TransSpecies.eigenval',header=FALSE)
relval = relval %>% mutate(VE = V1/sum(V1),label = paste0('Autosome: PC',row_number(),' ',signif(VE,3)*100,'%'))

#only grab the first axis, scale PC1 between 0 - 1
pcd = left_join(eigs,mdk)
pcl1 = pcd %>% select(ID,contains(c('PC','Plumage','Species_Latin','Shape'))) %>% pivot_longer(!c(ID,Plumage,PlumageColor,Species_Latin,Shape)) %>% 
  filter(!grepl('327',ID)) %>% 
  group_by(name) %>% 
  mutate(centered_value = value - min(value),
         scaled_value = centered_value / (max(value) - min(value)),
         axis = gsub('PC','',name),
         compartment = relval$label[1]) %>% mutate_at('axis',as.numeric)

#also W chr 
eigs = read.table('chr_W.IF-GF-MM2-BP-ANN-AC2.TransSpecies.eigenvec',header=FALSE) %>% select(-V2)
names(eigs) = c('ID',paste0('PC',seq(1,ncol(eigs)-1)))
eigs = eigs[1:2]

#import eigenvalues
relval = read.table('chr_W.IF-GF-MM2-BP-ANN-AC2.TransSpecies.eigenval',header=FALSE)
relval = relval %>% mutate(VE = V1/sum(V1),label = paste0('W: PC',row_number(),' ',signif(VE,3)*100,'%'))

#only grab the first axis
pcd = left_join(eigs,mdk)
pcl2 = pcd %>% select(ID,contains(c('PC','Plumage','Species_Latin','Shape'))) %>% pivot_longer(!c(ID,Plumage,PlumageColor,Species_Latin,Shape)) %>% 
  filter(!grepl('327',ID)) %>% 
  group_by(name) %>% 
  mutate(centered_value = value - min(value),
         scaled_value = centered_value / (max(value) - min(value)),
         axis = gsub('PC','',name),
         compartment = relval$label[1]) %>% mutate_at('axis',as.numeric)

#bind the autosome and W data together
pcl = rbind(pcl1,pcl2)

#perform wilcoxon rank sum test on the PC axes for both species and plumage 
p_values <- pcl %>%
  group_by(axis,compartment) %>%
  summarise(species = wilcox.test(scaled_value ~ Species_Latin)$p.value,
            plumage = wilcox.test(scaled_value ~ Plumage)$p.value) %>%
  pivot_longer(!c(axis,compartment))

# Apply Benjamini-Hochberg correction to the p-values
p_values_adjusted <- p_values %>%
  mutate(padj = p.adjust(value, method = "bonferroni",n = 4),
         significance = ifelse(padj < 0.05,'*','n.s.'))
pdat = p_values_adjusted

#make the plot !
pp = pcl %>% ggplot(aes(x=compartment,y=scaled_value,fill=Plumage,shape=Species_Latin))+
  geom_point(stroke=0.5,size=3,position=position_dodge(width=0.5))+
  geom_text(data = p_values_adjusted,
            aes(x = compartment, y = 1, col = name, size=significance,label = significance),
            inherit.aes = FALSE,vjust=-0.5,position=position_dodge(width=0.5)) +
  scale_fill_manual(values=pcl$PlumageColor,breaks=pcl$Plumage)+xlab('')+ylab('Scaled PC Scores')+
  scale_shape_manual(values=pcl$Shape,breaks=pcl$Species_Latin)+
  scale_color_manual(values=c('grey20','grey60'))+
  scale_size_manual(values=c(6,4))+
  theme_bw() +
  theme(legend.position='none')
pp

pdf('Comparison_A-W_SidebySide_2023AUG03-WIDE.pdf',height=2,width=5.5)
pp
dev.off()


counter = 0 
pdat=NULL
pcdat=NULL
RUNS= c('Autosomes.IF-GF-MM2-BP-ANN-AC2.TransSpecies','chr_W.IF-GF-MM2-BP-ANN-AC2.TransSpecies')
for (RUN in RUNS) {
  
  counter = counter + 1   
  #Read metadata
  mdk = read.table('../EntireMetadata.txt',header=TRUE,comment.char='',sep='\t')
  mdk$Species_Latin = factor(mdk$Species_Latin,levels=c('C. canorus','C. optatus','C. poliocephalus'))
  #import eigenvectors
  eigs = read.table(paste0(RUN,'.eigenvec'),header=FALSE) %>% select(-V2)
  names(eigs) = c('ID',paste0('PC',seq(1,ncol(eigs)-1)))
  eigs = eigs[1:7]
  #only grab the first 6 axes
  pcd = left_join(eigs,mdk)
  #import eigenvalues
  relval = read.table(paste0(RUN,'.eigenval'),header=FALSE)
  relval = relval %>% mutate(VE = V1/sum(V1),
                             label = paste0('PC',row_number(),' ',signif(VE,3)*100,'%'))
  
  #take a peek at first 2 axes
  pcd %>% ggplot(aes(x=PC1,y=PC2,fill=Plumage,shape=Species_Latin))+
    geom_point(size=3)+
    scale_shape_manual(values=pcd$Shape,breaks=pcd$Species_Latin)+
    scale_fill_manual(values=pcd$PlumageColor,breaks=pcd$Plumage)+
    xlab(paste0('PC1: ',round(relval$VE[[1]],3)*100,'%'))+ylab(paste0('PC2: ',round(relval$VE[[2]],3)*100,'%'))+
    theme_classic()+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(legend.text = element_text(size = 8),legend.title = element_text(size = 8),legend.key.size = unit(0.3, 'cm'),legend.position="right")
  
  #eigenvalues
  pcl = pcd %>% select(ID,contains(c('PC','Plumage','Species_Latin','Shape'))) %>% pivot_longer(!c(ID,Plumage,PlumageColor,Species_Latin,Shape)) %>% 
    filter(!grepl('327',ID)) %>% 
    group_by(name) %>% 
    mutate(centered_value = value - min(value),
           scaled_value = centered_value / (max(value) - min(value)),
           axis = gsub('PC','',name)) %>% mutate_at('axis',as.numeric)
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
  
  # Perform Wilcoxon rank sum test for each axis, separately for Species_Latin and Plumage
  p_values <- pcl %>%
    group_by(axis) %>%
    summarise(species = wilcox.test(scaled_value ~ Species_Latin)$p.value,
              plumage = wilcox.test(scaled_value ~ Plumage)$p.value) %>%
    pivot_longer(!(axis))
  
  # Apply Benjamini-Hochberg correction to the p-values
  p_values_adjusted <- p_values %>%
    mutate(padj = p.adjust(value, method = "bonferroni",n = 24),
           significance = ifelse(padj < 0.05,'*','n.s.'),
           chromosome = RUN)
  pdat = rbind(pdat,p_values_adjusted)
  
  pp = pcl %>% ggplot(aes(x=as.factor(axis),y=scaled_value,fill=Plumage,shape=Species_Latin))+
    geom_point(stroke=0.5,size=3,position=position_dodge(width=0.5))+
    geom_text(data = p_values_adjusted,
              aes(x = axis, y = 1, col = name, size=significance,label = significance),
              inherit.aes = FALSE,vjust=0,position=position_dodge(width=0.5)) +
    #annotate("text", x = -0.01, y = 1.1, label = "Plumage", color = "grey50",hjust=-.5,vjust=0) +
    scale_fill_manual(values=pcl$PlumageColor,breaks=pcl$Plumage)+xlab('')+ylab('Scaled PC Scores')+
    scale_shape_manual(values=pcl$Shape,breaks=pcl$Species_Latin)+
    scale_color_manual(values=c('grey20','grey60'))+
    scale_size_manual(values=c(6,4))+
    scale_x_discrete(labels = relval[1:6,3]) +
    theme_classic() 
  assign(paste0('p',counter),pp)
  pcdat = rbind(pcdat,pcl %>% mutate(Analysis = RUN))
}

pdf('Comparison_A-W_2023AUG03.pdf',height=4,width=7)
ggarrange(p2,p1,nrow=2,ncol=1,common.legend = TRUE)
dev.off()

#save data
pcdat = pcdat %>% mutate(Analysis = gsub('\\..*','',Analysis))
write.table(pcdat,file='PC_Scores_2023AUG03.txt',quote=F,sep='\t',row.names=F)
pwdat = pdat %>% mutate(Analysis = gsub('\\..*','',chromosome)) %>% select(-c(chromosome,value,significance)) %>% ungroup %>% pivot_wider(names_from=name,values_from=padj)
write.table(pwdat,file='Pvalues_2023AUG03.txt',quote=F,sep='\t',row.names=F)

#check normality in case
# Function to check ANOVA assumptions
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

# Check assumptions for each axis
unique_axes <- unique(data$axis)
for (axis in unique_axes) {
  cat("Checking assumptions for", axis, "\n")
  check_anova_assumptions(data, axis)
}
