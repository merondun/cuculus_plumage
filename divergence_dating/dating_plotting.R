#### Divergence Dating
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/divergence/out')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(gghalves)
library(boot)
library(ggdist)

#read in data 
fst = read.table('Sensitivity_DXY.input',header=TRUE) #or start here 
fst = fst %>% mutate(chrnum = ifelse(chr == 'chr_Z',40,ifelse(chr == 'chr_W',41,ifelse(chr== 'chr_MT',42,chr))),
                     chrnum = gsub('chr_','',chrnum),
                     chr = gsub('chr_','',chr)) %>% 
  mutate_at('chrnum',as.numeric) %>% filter(chrnum < 21 | chrnum > 39)
fst %>% head

### First, calculate average pi / dxy by chromosome, and then calculate average Da by chromosome 
dac = fst %>% group_by(chr,Group,chrnum,Iteration,Mask) %>% summarize(piA = mean(piA,na.rm=TRUE),
                                                                      piB = mean(piB,na.rm=TRUE),
                                                                      dxy = mean(dxy,na.rm=TRUE),
                                                                      Fst = mean(Fst,na.rm=TRUE),
                                                                      Da = dxy - ((piA+piB)/2),
                                                                      Da = ifelse(Da < 0 , 0, Da)) #with undiverged clades sometimes Da negative, change to 0

#if W, divide mu by 2, if Z divide by 3/4, if autosome leave alone 
da = dac %>% group_by(chr,Group,chrnum,Iteration,Mask) %>% 
  mutate(time = ifelse(chr == 'W', Da / (2*9.08e-09*(1/2)),
                       ifelse(chr == 'Z', Da / (2*9.08e-09),
                              Da / (2*9.08e-09))))
#add some aesthetic changes 
da = da %>% mutate(AvZ = ifelse(chr == 'W','W',ifelse(chr=='Z','Z','Autosome')),
                   Shape = ifelse(Group == 'OG_CP' | Group == 'CP_CG','Outgroup',
                                  ifelse(Group == 'OH_CH' | Group == 'OG_CG','canorus_optatus','Plumage')))
#only keep sensible pairiwse group comparisons 
dak = da %>% filter(chr != 'MT')
dak$Group = factor(dak$Group,levels=c('HG_HH','CG_CH','OG_OH','CG_OG','CH_OH'))
chord = dak %>% ungroup %>% select(chr,chrnum) %>% unique %>% arrange(chrnum)
dak$chr = factor(dak$chr,levels=chord$chr)
dak = dak %>% drop_na(Group) %>% ungroup

#plot chromosome estimates 
chrplot = dak %>% ungroup %>% 
  ggplot(aes(x=time, y=chr,col=Group,shape=Mask,size=Mask)) + 
  geom_jitter(width = 0.25)+
  scale_shape_manual(values=c(6,3))+
  scale_size_manual(values=c(3,1))+
  xlab('Generations')+ylab('')+
  scale_fill_viridis(discrete=TRUE)+
  scale_color_viridis(discrete=TRUE)+
  scale_x_continuous(breaks = c(0,25000,50000,75000,100000),
                     labels = c('0','25K','50K','75K','100K')) +
  theme_bw()
pdf('../Time_Estimates_Sensitivity_2023JUNE12.pdf',height=6,width=6)
chrplot
dev.off()

### Bootstrapping 95% CIs 
set.seed(123) # Set seed for reproducibility
#function for calculating Da 
calculate_Da <- function(data, indices) {
  resampled_data <- data[indices, ]
  piA_mean <- mean(resampled_data$piA, na.rm = TRUE)
  piB_mean <- mean(resampled_data$piB, na.rm = TRUE)
  dxy_mean <- mean(resampled_data$dxy, na.rm = TRUE)
  Da <- dxy_mean - ((piA_mean + piB_mean) / 2)
  return(Da)
}

#subset 50KB data, we will have 11 iterations with 0 = All samples 
bootdat = fst %>% filter(Window == 50000) %>% mutate(Iteration = ifelse(Mask == 'All',0,Iteration))
bootres = NULL
for (it in unique(bootdat$Iteration)){
  cat('Working on iteration: ',it,'\n')
  
  bd = bootdat %>% filter(Iteration == it)
  # Apply bootstrapping for each 'chr' and store the results in a list
  bootstrapped_results <- bd %>% 
    group_by(AvZ,Group) %>% 
    group_split() %>% 
    lapply(function(df) {
      boot(df, calculate_Da, R = 1000) # R is the number of bootstrap replicates; adjust as needed
    })
  
  # Combine the results into a single data frame
  bootstrapped_results_df <- lapply(bootstrapped_results, function(x) {
    data.frame(AvZ = unique(x$data$AvZ),
               Group = unique(x$data$Group),
               Da = x$t0,
               Iteration = it,
               Da_lower = quantile(x$t, 0.025),
               Da_upper = quantile(x$t, 0.975))
  }) %>% bind_rows()
  
  bootres = rbind(bootres,bootstrapped_results_df)
}

#keep only sensible comparisons 
write.table(bootres,file='Bootstrapped_Da.txt',quote=F,sep='\t',row.names=F)
bootres = read.table('Bootstrapped_Da.txt',header=TRUE)
bdak = bootres
bdak$Group = factor(bdak$Group,levels=c('HG_HH','CG_CH','OG_OH','CG_OG','CH_OH'))
bdak = bdak %>% drop_na(Group) %>% ungroup
bdak = bdak %>% 
  mutate(Da = ifelse(Da < 0 , 0, Da), #set Da and upper and lower to 0 if Da is les than 0 
         Da_lower = ifelse(Da_lower < 0 , 0, Da_lower),
         Da_upper = ifelse(Da_upper < 0 , 0, Da_upper),
         time = ifelse(AvZ == 'W', Da / (2*9.08e-09*(1/2)),
                       ifelse(AvZ == 'Z', Da / (2*9.08e-09),
                              Da / (2*9.08e-09))),
         time_lower = ifelse(AvZ == 'W', Da_lower /(2*9.08e-09*(1/2)),
                             ifelse(AvZ == 'Z', Da_lower / (2*9.08e-09),
                                    Da_lower / (2*9.08e-09))),
         time_upper = ifelse(AvZ == 'W', Da_upper /(2*9.08e-09*(1/2)),
                             ifelse(AvZ == 'Z', Da_upper / (2*9.08e-09),
                                    Da_upper / (2*9.08e-09))))
bdak$AvZ = factor(bdak$AvZ,levels=c('Autosome','Z','W'))
bdak = bdak %>% mutate(Mask = ifelse(Iteration == 0,'All','Subset'),
                       Facet = ifelse(Group == 'HG_HH' | Group == 'CG_CH' | Group == 'OG_OH','Plumage Comparison','Species Comparison'))
#plot 
dtimes = ggplot(bdak, aes(x = AvZ, y = time,fill=Group,col=Group)) +
  #geom_boxplot(width = .5, alpha = 0.3) +
  geom_point(data = bdak %>% filter(Mask == 'All'),aes(col=Group),stroke=1,shape=4,size=3,position=position_dodge(width=0.5)) +
  geom_errorbar(data = bdak %>% filter(Mask == 'All'),aes(col=Group,ymin=time_lower,ymax=time_upper),size=1,width=0.5,position=position_dodge(width=0.5)) +
  #stat_halfeye(adjust = .5,width = .25,.width = 0,justification = -1.25, point_colour = NA,alpha = 0.8,normalize='groups')+
  labs(x = "",
       y = "Time (Generations) \nbootstrapped mean with subsampling") +
  theme_bw(base_size=14) +
  scale_y_continuous(
    limits = c(0,10e4),
    breaks = c(0,2.5e4,5e4,7.5e4,10e4),
    labels = c('0','25K','50K','75K','100K')) +
  scale_shape_manual(values=c(8,3))+
  facet_grid(.~Facet,)+
  scale_fill_viridis(discrete=TRUE) + 
  scale_color_viridis(discrete=TRUE) + 
  theme(legend.position='top')

pdf('../Divergence_Times-Simple_2023JUNE12.pdf',height=5,width=7)
dtimes
dev.off()

bdak %>% filter((Group =='CG_CH' | Group == 'OG_OH' ) & AvZ == 'W')

### Show the time estimates for each window instead, or across window sizes 
dak = fst
dak$Group = factor(dak$Group,levels=c('HG_HH','CG_CH','OG_OH','CG_OG','CH_OH'))
chord = dak %>% ungroup %>% select(chr,chrnum) %>% unique %>% arrange(chrnum)
dak$chr = factor(dak$chr,levels=chord$chr)
dak = dak %>% drop_na(Group) %>% ungroup
dak = dak %>% 
  mutate(Da = dxy - ((piA + piB) / 2),
         Da = ifelse(Da < 0 , 0, Da),
         time = ifelse(chr == 'W', Da / (2*9.08e-09*(1/2)),
                       ifelse(chr == 'Z', Da / (2*9.08e-09),
                              Da / (2*9.08e-09))))
dak$AvZ = factor(dak$AvZ,levels=c('Autosome','Z','W'))
write.table(bdak,file='../Da_Estimates_Generations_Results_2023JUNE12.txt',quote=F,sep='\t',row.names=F)

#plot window boxplot and density time estimates 
divs = dak %>% 
  mutate(time = ifelse(time > 1.5e6, 1.5e6,time)) %>%  #add ceilingat 2M 
  ggplot(aes(y=time, x=AvZ,fill=Group))  + 
  stat_halfeye(adjust = .5,width = .25,.width = 0,justification = -1.5, point_colour = NA,alpha = 0.5,normalize='groups')+
  geom_boxplot(width = .5,outlier.shape = NA, alpha = 0.3) +
  ylab('Time')+xlab('')+
  #scale_y_continuous(limits = c(0,1.5e6),breaks = c(100000, 250000, 500000, 1000000, 1500000),
  #                   labels = c("100K", "250K", "500K", "1M", "1.5M")) +
  scale_fill_viridis(discrete=TRUE)+
  scale_color_viridis(discrete=TRUE)+
  theme_bw()+
  coord_flip()
divs
#pdf('../Divergence_Dating_250KB.pdf',height=6,width=7)
png('../Divergence_Dating_50KB_2023JUNE12.png',units='in',res=600,height=6,width=7)
divs
dev.off()
