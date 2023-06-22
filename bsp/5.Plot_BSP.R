#### Plot BSP
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/bsp/plot')
.libPaths('~/miniconda3/envs/R/lib/R/library')
library(tidyverse)
.libPaths('~/miniconda3/envs/mtdna/lib/R/library')
library(RColorBrewer)
library(scales)
bsp = read.table('BSP.indat',header=FALSE)
names(bsp) = c('Time','Mean','Median','Upper','Lower','Group')
bspp = bsp %>% ggplot(aes(x=Time,col=Group))+
  geom_line(aes(y=Mean))+
  geom_line(aes(y=Upper),lty=2)+
  geom_line(aes(y=Lower),lty=2)+
  theme_bw()+
  scale_color_manual(values=brewer.pal(6,'Set2'))+
  facet_wrap(Group~.,scales='free')+
  scale_y_log10(expression(paste(N[e], " (M)")),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = function(l) {trans = round(l/100000,2) ; paste0(trans, "M")})+
  scale_x_log10('Years (M)',
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = function(l) {trans = round(l/100000,2) ; paste0(trans, "M")})
bspp

png('BSP.png',height=6,width=8,units='in',res=600)
bspp
dev.off()
