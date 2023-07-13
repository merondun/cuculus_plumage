#### Plot BSP
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/bsp/2023JUNE_GENES/nex')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(RColorBrewer)
library(scales)
bsp = read.table('BSP_Input.txt',header=TRUE) %>% as_tibble 
names(bsp) = c('Time','Mean','Median','Upper','Lower','Group')
bsp = bsp %>% mutate(Clade = ifelse(grepl('CG',Group),'canorus grey',
                                    ifelse(grepl('CH',Group),'canorus hepatic',
                                           ifelse(grepl('OG',Group),'optatus grey',
                                                  ifelse(grepl('OH',Group),'optatus hepatic',
                                                         ifelse(grepl('HH',Group),'hungarian hepatic','hungarian grey'))))),
                     Sympatry = ifelse(grepl('HH|HG',Group),'Hungarian','All'))

#Summary stats of the Mean across all time points 
bsp %>% group_by(Group,Clade,Sympatry) %>% summarize(mean = mean(Mean,na.rm=TRUE))

#Contemporary
bsp %>% group_by(Group,Clade,Sympatry) %>% filter(Time == 0)

#Plot
bspp = bsp %>% ggplot(aes(x=Time,col=Clade))+
  geom_ribbon(aes(fill=Clade,ymin=Lower,ymax=Upper),alpha=0.2,col=NA)+
  geom_line(aes(y=Mean),lwd=1)+
  theme_bw()+
  scale_color_manual(values=brewer.pal(6,'Dark2'))+
  scale_fill_manual(values=brewer.pal(6,'Dark2'))+
  facet_grid(Sympatry~.,scales='free')+
  scale_y_log10(expression(paste(N[e], " (M)")),n.breaks=5,
                #breaks = trans_breaks("log10", function(x) 10^x),
                labels = function(l) {trans = round(l/1000000,2) ; paste0(trans, "M")})+
  scale_x_log10('Generations (g)',n.breaks=5,
                #breaks = trans_breaks("log10", function(x) 10^x),
                labels = function(l) {trans = round(l/1000,2) ; paste0(trans, "K")})
bspp

png('BSP_2023JULY3.png',height=4,width=4.5,units='in',res=600)
pdf('BSP_2023JUNE26.pdf',height=2.5,width=3.5)
bspp
dev.off()

