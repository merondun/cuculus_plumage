#### SNP summaries
setwd('/dss/dsshome1/lxc07/di39dux/merondun/cuculus_plumage/general')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(ggpubr)

#read in data
counts = read.table('SNP_Numbers.txt',header=TRUE)
#assign nice names, create a numeric chr variable for ordering
c1 = counts %>% 
  pivot_longer(!c(CHR))  %>% 
  mutate(Filter = ifelse(grepl('acfilt',name),'Without Singletons',
                                           ifelse(grepl('bpfilter',name),'Excluding SNPs within 5bp',
                                                  ifelse(grepl('genofilter',name),'FORMAT Filters',
                                                         ifelse(grepl('infofilter',name),'INFO Filters',
                                                                ifelse(grepl('site',name),'Total Sites',
                                                                       ifelse(grepl('snps',name),'Raw SNPs','NO')))))),
         CHR = gsub('chr_','',CHR),
         chord = ifelse(CHR == 'Z',40,ifelse(CHR=='W',41,ifelse(CHR == 'MT',42,CHR)))) %>% 
  mutate_at('chord',as.numeric)
#aesthetic ordering 
chrlevel = c1 %>% select(CHR,chord) %>% unique %>% arrange(chord)
level = data.frame(lv = unique(c1$Filter)) #it's already in order somehow, very nice! 
c1$Filter = factor(c1$Filter,levels=level$lv)
c1$CHR = factor(c1$CHR,levels=chrlevel$CHR)
cols = viridis(5)

#calculate average across all samples, then calculate % difference 
alls = c1 %>% filter(Filter == 'Total Sites') %>% dplyr::rename(total = value) %>% select(-c(name,Filter))
c2 = left_join(c1,alls)
c2 = c2 %>% mutate(Pct_Diff = value / total * 100)

#plot 
c3 = c2 %>%  
  filter(Filter != 'Total Sites') %>% 
  ggplot(aes(y=CHR,x=Pct_Diff,col=Filter,shape=Filter))+
  geom_label(data = c2 %>% filter(Filter == 'Excluding SNPs within 5bp'),aes(y=CHR,x=14,label=paste0(round(value/1000,0),'K')),size=2,fill=cols[4],inherit.aes = FALSE)+
  geom_label(data = c2 %>% filter(Filter == 'Without Singletons'),aes(y=CHR,x=15,label=paste0(round(value/1000,0),'K')),size=2,fill=cols[5],inherit.aes = FALSE)+
  geom_point(size=3)+
  scale_x_continuous(limits = c(0,15),breaks = c(0,5,10),
                     labels = c('0','5%','10%')) +
  ylab('')+xlab('Percent SNPs Retained After Filter') +
  scale_fill_viridis(discrete=TRUE)+
  guides(color = guide_legend(nrow = 2))+
  scale_color_viridis(discrete=TRUE)+
  theme_bw()+
  theme(legend.position='top',legend.key.size = unit(0.05,'cm'))

#pdf('SNP_Counts_2023JUNE16.pdf',height=7,width=6)
png('SNP_Counts_2023JUNE16.png',height=7,width=6,units='in',res=600)
c3
dev.off()

#And SNP quality metrics
stats = read.table('SNP_Stats.txt',header=TRUE)
#assign nice names, classify as Autosome / Z / W, etc.
s1 = stats %>% 
  dplyr::rename(CHR = CHROM) %>% 
  pivot_longer(!c(CHR,ID)) %>% 
  mutate(Filter = ifelse(grepl('nREF',name),'Number Reference Alleles',
                         ifelse(grepl('nALT',name),'Number Alternate Alleles',
                                ifelse(grepl('nHET',name),'Number Heterozygous',
                                       ifelse(grepl('nTs',name),'Number Transitions',
                                              ifelse(grepl('nTv',name),'Number Transversions',
                                                     ifelse(grepl('avgDP',name),'Average Depth',
                                                            ifelse(grepl('Singletons',name),'Number Singletons',
                                                                   ifelse(grepl('Missing_Sites',name),'Number Missing Genotypes','NO')))))))),
         CHR = gsub('chr_','',CHR),
         chord = ifelse(CHR == 'Z',40,ifelse(CHR=='W',41,ifelse(CHR == 'MT',42,CHR))),
         AvZ = ifelse(CHR == 'Z','Z',ifelse(CHR=='W','W',ifelse(CHR == 'MT','MT','Autosome'))),
         Species = ifelse(grepl('_CO_',ID),'C. optatus',
                          ifelse(grepl('_CC_',ID),'C. canorus',
                                 ifelse(grepl('_CM_',ID),'C. micropterus',
                                        ifelse(grepl('_CP_',ID),'C. poliocephalus','NO'))))) %>% 
         mutate_at('chord',as.numeric)
s2 = s1 %>% group_by(AvZ,ID,Species,Filter) %>% summarize(mean = mean(value))
s2$AvZ = factor(s2$AvZ,levels=c('Autosome','Z','W','MT'))

#Calculate TsTv, and then also average for each variable to examine sample disparities
s3 = s2 %>% filter(Filter == 'Number Transitions' | Filter == 'Number Transversions') %>% 
  group_by(AvZ,ID,Species) %>% summarize(mean = mean[Filter == 'Number Transitions'] / mean[Filter == 'Number Transversions'],
                                 Filter = 'TsTv') 
s4 = rbind(s3,s2 %>% filter(Filter != 'Number Transitions' & Filter != 'Number Transversions'))
avg = s4 %>% ungroup %>% group_by(AvZ,Filter,Species) %>% summarize(avg = mean(mean))
s5 = left_join(s4,avg)

#Calculate percent difference from the mean 
s5 = s5 %>% mutate(Pct_Diff = (mean - avg) / avg * 100)

#plot 
RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)
s6 = s5 %>% filter(AvZ != 'MT' & !grepl('_CP_|_CM_|252|276|277',ID)) %>% 
  filter(!grepl('Reference|Alternate',Filter)) %>% 
  ggplot(aes(y=ID,x=Pct_Diff,col=Filter,shape=Filter))+
  geom_vline(xintercept=0,lty=2)+
  geom_point(alpha=0.9)+
  facet_grid(Species~AvZ,scales='free',space = 'free_y')+
  scale_color_viridis(discrete=TRUE,option='turbo')+
  theme_bw()+xlab('Percent Difference')+ylab('')+
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position='top',axis.text.y = element_text(size = 5))
  
#pdf('SNP_Stats_2023JUNE16.pdf',height=8,width=6)
png('SNP_Stats_2023JUNE16.png',height=8,width=8,units='in',res=600)
s6
dev.off()


