#### Plot randomized FST 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/divergence/out')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(gghalves)
library(boot)
library(ggdist)
library(LICORS)
library(broom)
library(purrr)
library(ggpubr)
library(gmodels)

fst = read.table('Randomized_Fst.dxy',header=TRUE)
names(fst) = c('chr','start','end','mid','sites','piA','piB','dxy','Fst','Group')
#drop NA dxy values, missing data for one clade
fstd = fst %>% drop_na(dxy,piA,piB)
#at this point, all NAs for FST are FST=0, since there is no diversity at that site
nas = fstd %>% filter(is.na(Fst)) #check that's true
fstd = fstd %>% replace_na(list(Fst=0))
#add ceiling to fst since for plotting negative fst isn't sensible for plotting
fstd = fstd %>% mutate(Fst = threshold(Fst,min=0,max=1))

#analyze via stats
fs = fstd %>% select(chr,Fst,Group) %>% as_tibble
fs$Fst = as.numeric(fs$Fst)
fs$chr = gsub('chr_','',fs$chr)
fs = fs %>% mutate(t1 = gsub('rs','sub@',Group)) %>% 
  separate(t1,into=c('trash','NumberSubsampled'),sep='@')
fs = fs %>% mutate(Replicate = ifelse(Group == 'HH_HG','True Color Comparison','Randomized'))
fs = fs %>% mutate(Treatment = ifelse(grepl('rs',Group),'Sequential_Rufous_Removal',
                                      ifelse(grepl('rf',Group),'Full_Sample_Shuffling',
                                             ifelse(grepl('r',Group),'Subsampled_Shuffling','Color_Comparison'))))
fsx = fs %>% group_by(chr,Group,NumberSubsampled,Replicate,Treatment) %>% 
  summarize(m = ci(Fst)[1],
            lwr = ci(Fst)[2],
            upr = ci(Fst)[3])

#plot subsampled first by itself
fss = fsx %>% filter(grepl('Seque|Color',Treatment)) %>% mutate_at('NumberSubsampled',as.numeric) %>% replace_na(list(NumberSubsampled = 0))
fsp = fss %>% ggplot(aes(x=NumberSubsampled,y=m,ymin=lwr,ymax=upr,col=Replicate,group=Group))+
  geom_point(stroke=1,shape=16,size=2,position=position_dodge(width=0.5)) +
  geom_errorbar(size=1,width=0.5,position=position_dodge(width=0.5)) +
  ggtitle('Sequentially Shuffling n Rufous as Grey')+ylab('Mean and 95% CI FST')+
  facet_grid(.~chr,scales='free')+xlab('Number Rufous and Grey Swapped')+
  theme_bw()
fsp

#plot both the full n = 20 comparison 
fsr = fsx %>% filter(grepl('Full|Color',Treatment))
frp = fsr %>% ggplot(aes(x=chr,y=m,ymin=lwr,ymax=upr,col=Replicate,group=Group))+
  geom_point(stroke=1,shape=16,size=2,position=position_dodge(width=0.5)) +
  geom_errorbar(size=1,width=0.5,position=position_dodge(width=0.5)) +
  ggtitle('Full Sample Shuffling')+
  facet_grid(.~chr,scales='free')+xlab('Shuffled Iteration')+ylab('Mean and 95% CI FST')+
  theme_bw()
frp

#plot both the full n = 20 comparison and the n = 5 for both groups 
fsv = fsx %>% filter(grepl('Sub|Color',Treatment))
fvp = fsv %>% ggplot(aes(x=chr,y=m,ymin=lwr,ymax=upr,col=Replicate,group=Group))+
  geom_point(stroke=1,shape=16,size=2,position=position_dodge(width=0.5)) +
  geom_errorbar(size=1,width=0.5,position=position_dodge(width=0.5)) +
  ggtitle('Subsampled Shuffling (n=5)')+
  facet_grid(.~chr,scales='free')+xlab('Shuffled Iteration')+ylab('Mean and 95% CI FST')+
  theme_bw()
fvp

ggarrange(fsp,frp,fvp,common.legend=TRUE,nrow=3,ncol=1)
pdf('FST_Sensitivity_WChromosome_2023JUNE14.pdf',height=7,width=5)
ggarrange(fsp,frp,common.legend=TRUE,nrow=2,ncol=1)
dev.off()
