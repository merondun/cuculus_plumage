#### Divergence Dating
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/martin_dxy_singletons/out')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(gghalves)
library(boot)
library(ggdist)
library(LICORS)
library(gmodels)

#read in data 
fst = read.table('Sensitivity_DXY.input',header=TRUE) #or start here 
fst = fst %>% mutate(chrnum = ifelse(chr == 'chr_Z',40,ifelse(chr == 'chr_W',41,ifelse(chr== 'chr_MT',42,chr))),
                     chrnum = gsub('chr_','',chrnum),
                     chr = gsub('chr_','',chr)) %>% 
  mutate_at('chrnum',as.numeric) %>% filter(chrnum < 21 | chrnum > 39)  %>%
  dplyr::rename(dxy=DXY,Fst=FST)
fst %>% head

### First, calculate average pi / dxy by chromosome, and then calculate average Da by chromosome 
dac = fst %>% group_by(chr,Group,chrnum) %>% summarize(piA = mean(piA,na.rm=TRUE),
                                                                      piB = mean(piB,na.rm=TRUE),
                                                                      dxy = mean(dxy,na.rm=TRUE),
                                                                      Fst = mean(Fst,na.rm=TRUE),
                                                                      Da = dxy - ((piA+piB)/2),
                                                                      Da = ifelse(Da < 0 , 0, Da)) #with undiverged clades sometimes Da negative, change to 0

#if W, divide mu by 2, if Z divide by 3/4, if autosome leave alone 
mu = 1.01e-08
da = dac %>% group_by(chr,Group,chrnum) %>% 
  mutate(time = ifelse(chr == 'W', Da / (2*mu*(1/2)),
                       ifelse(chr == 'Z', Da / (2*mu),
                              Da / (2*mu))))
#add some aesthetic changes 
dak = da %>% mutate(AvZ = ifelse(chr == 'W','W',ifelse(chr=='Z','Z','Autosome')),
                   Shape = ifelse(Group == 'OG_CP' | Group == 'CP_CG','Outgroup',
                                  ifelse(Group == 'OG_CG' | Group == 'OH_CH','canorus_optatus','Plumage')))
#only keep sensible pairiwse group comparisons 
dak$Group = factor(dak$Group,levels=c('HH_HG','CG_CH','OG_OH','OG_CG','OH_CH'))
chord = dak %>% ungroup %>% select(chr,chrnum) %>% unique %>% arrange(chrnum)
dak$chr = factor(dak$chr,levels=chord$chr)
dak = dak %>% drop_na(Group) %>% ungroup

#only plot species comparison
species = dak %>% filter((Group == 'OG_CG' | Group == 'OH_CH') & AvZ == 'Autosome')
speciestimes = species %>% group_by(Group) %>% summarize(mean = ci(time)[1],
                                                         lwr = ci(time)[2],
                                                         upr = ci(time)[3])
speciestimes
# Group   mean    lwr    upr
# <fct>  <dbl>  <dbl>  <dbl>
#   1 OG_CG 50988. 46481. 55495.
# 2 OH_CH 45107. 40590. 49624.

dtimes = ggplot(species, aes(x = Group, y = threshold(time,max=10e4),fill=Group)) +
  geom_violin(width = .5, alpha = 0.2,draw_quantiles = c(0.25,0.75)) +
  geom_label(data=speciestimes,
             aes(x=Group,col=Group,group=Group,y=95000,label=paste0(round(mean/1000,0),'K (',round(lwr/1000,0),'K - ',round(upr/1000),'K)')),
             position=position_dodge(width=0.5),fill='white',size=2) + 
  labs(x = "", y = "Time (Generations)") +
  theme_bw() +
  scale_y_continuous(
    limits = c(0,10e4),
    breaks = c(0,2.5e4,5e4,7.5e4,10e4),
    labels = c('0','25K','50K','75K','100K')) +
  scale_fill_manual(values=c('black','orange')) + 
  scale_color_manual(values=c('black','orange')) + 
  theme(legend.position='none')
dtimes
pdf('~/merondun/cuculus_plumage/divergence_dating/Time_Estimates_Autosomes-Singletons_2024JAN12.pdf',height=4,width=2.5)
dtimes
dev.off()

