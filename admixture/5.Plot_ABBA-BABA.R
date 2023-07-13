#### GW ABBA-BABA from https://github.com/simonhmartin/tutorials/tree/master/ABBA_BABA_whole_genome
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/ABBA/')
library(tidyverse)
library(viridis)
library(data.table)
source("~/modules/genomics_general/jackknife.R")

D.stat <- function(p1, p2, p3) {
  ABBA <- (1 - p1) * p2 * p3
  BABA <- p1 * (1 - p2) * p3
  (sum(ABBA, na.rm=T) - sum(BABA, na.rm=T)) / (sum(ABBA, na.rm=T) + sum(BABA, na.rm=T))
}

groups = c('CH','CG')
AvZs = c('W','Z','Autosome')
res = NULL

for (group in groups){
  freq_table_full = read.table(paste0(group,'.ABBA.txt'),header=FALSE)
  names(freq_table_full) = c('scaffold','position','OG','OH',group,'CP')
  freq_table_full = freq_table_full %>% mutate(AvZ = ifelse(scaffold == 'chr_Z','Z',ifelse(scaffold == 'chr_W','W','Autosome')))
  
  for (comp in AvZs) {
    cat('Working on group : ',group,' for compartment: ',comp,'\n')
    freq_table = freq_table_full %>% filter(AvZ == comp)
    
    P1 = "OG"
    P2 = "OH"
    P3 = group
    
    # Calculate D statistic
    D = D.stat(freq_table[,P1], freq_table[,P2], freq_table[,P3])
    
    # Divide for jackknifing
    block_indices <- get.block.indices(block_size=1e6,
                                       positions=freq_table$position,
                                       chromosomes=freq_table$scaffold)
    
    n_blocks <- length(block_indices)
    
    print(paste("Genome divided into", n_blocks, "blocks."))
    
    # Jackknife it
    D_jackknife <- block.jackknife(block_indices=block_indices,
                                   FUN=D.stat,
                                   freq_table[,P1], freq_table[,P2], freq_table[,P3])
    print(paste("D jackknife mean =", round(D_jackknife$mean,4)))
    
    D_Z <- D_jackknife$mean / D_jackknife$standard_error
    print(paste("D Z score = ", round(D_Z,3)))
    
    # Merge and bind results
    r = data.frame(Group = group, AvZ = comp, D = D, Djack = D_jackknife, Z = D_Z)
    res = rbind(r,res)
  }
}

write.table(res,file='Jackknife_Dstatistics_2023JUNE26.txt',quote=F,sep='\t',row.names=F)

res = read.table('Jackknife_Dstatistics_2023JUNE26.txt',header=TRUE)
ap = res %>% 
  #filter(AvZ != 'W') %>% 
  ggplot(aes(x=AvZ,y=Djack.mean,ymin=Djack.mean-Djack.standard_deviation,
                        ymax=Djack.mean+Djack.standard_deviation,
                        col=Group,label=paste0(signif(Djack.mean,3),'\n (Z = ',signif(Z,3),')')))+
  geom_point(position=position_dodge(width=0.75))+
  geom_text(aes(y=Djack.mean+Djack.standard_deviation + 0.06),position=position_dodge(width=0.75),size=1.75)+
  geom_errorbar(width=0.25,position=position_dodge(width=0.75))+
  xlab('')+ylab('D-Statistic')+
  scale_color_manual(guide = 'none',values=c('#000000','#e78f41'))+
  geom_hline(yintercept=0,lty=2)+
  theme_bw()

pdf('Dstatistics_2023JUNE26.pdf',height=5,width=4)
ap
dev.off()

