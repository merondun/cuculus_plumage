#### PCA from PLINK on
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/autosomes')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)

RUN = c('Autosomes.IF-GF-MM2-BP-ANN-AC2.HUNGARY-P50.10.1-LD')

#Read metadata
mdk = read.table('../EntireMetadata.txt',header=TRUE,comment.char='',sep='\t')
mdk$Species_Latin = factor(mdk$Species_Latin,levels=c('C. canorus','C. optatus','C. poliocephalus'))
#import eigenvectors
eigs = read.table(paste0(RUN,'.eigenvec'),header=FALSE) %>% select(-V2)
names(eigs) = c('ID',paste0('PC',seq(1,ncol(eigs)-1)))
eigs = eigs[1:5]
#only grab the first 10 axes
pcd = left_join(eigs,mdk)
#import eigenvalues
relval = read.table(paste0(RUN,'.eigenval'),header=FALSE)
relval = relval %>% mutate(VE = V1/sum(V1),
                           label = paste0('PC',row_number(),' ',signif(VE,3)*100,'%'))

#take a peek at first 2 axes
pc12 = pcd %>% ggplot(aes(x=PC1,y=PC2,fill=Plumage,shape=Species_Latin))+
  geom_point(size=3)+
  scale_shape_manual(values=pcd$Shape,breaks=pcd$Species_Latin)+
  scale_fill_manual(values=pcd$PlumageColor,breaks=pcd$Plumage)+
  xlab(paste0('PC1: ',round(relval$VE[[1]],3)*100,'%'))+ylab(paste0('PC2: ',round(relval$VE[[2]],3)*100,'%'))+
  theme_classic()+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(legend.position='none',
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.3, 'cm'))
pc12 

pdf('PCA_HUNGARY_2023AUG03.pdf',height=4,width=5)
pc12
dev.off()
