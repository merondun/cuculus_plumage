#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#### GW ABBA-BABA from https://github.com/simonhmartin/tutorials/tree/master/ABBA_BABA_whole_genome
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/ABBA/2024JAN')
library(tidyverse)
library(viridis)
library(data.table)
source("~/modules/genomics_general/jackknife.R")

#submit with group e.g. HUN_Parity
pop = args[1]

D.stat <- function(p1, p2, p3) {
  ABBA <- (1 - p1) * p2 * p3
  BABA <- p1 * (1 - p2) * p3
  (sum(ABBA, na.rm=T) - sum(BABA, na.rm=T)) / (sum(ABBA, na.rm=T) + sum(BABA, na.rm=T))
}

groups = c('CH','CG')
res = NULL

freq_table = read_tsv(paste0('ABBA/',pop,'_2024JAN17.txt')) %>% na.omit %>% as.data.frame

for (group in groups){
  cat('Working on group : ',group,'\n')

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
  r = data.frame(Group = group, Djack = D_jackknife, Z = D_Z)
  res = rbind(r,res)

}

#admixture proportion (f) using P3a and P3b "average proportion of foreign ancestry in any given haploid genome"
abba = function(p1, p2, p3) (1 - p1) * p2 * p3
baba = function(p1, p2, p3) p1 * (1 - p2) * p3
P1 = 'OG'
P2 = 'OH'
P3a = 'CG'
P3b = 'CH'

ABBA_1_2_3a = abba(freq_table[,P1], freq_table[,P2], freq_table[,P3a])
BABA_1_2_3a = baba(freq_table[,P1], freq_table[,P2], freq_table[,P3a])

ABBA_1_3b_3a = abba(freq_table[,P1], freq_table[,P3b], freq_table[,P3a])
BABA_1_3b_3a = baba(freq_table[,P1], freq_table[,P3b], freq_table[,P3a])

f = (sum(ABBA_1_2_3a) - sum(BABA_1_2_3a))/
  (sum(ABBA_1_3b_3a) - sum(BABA_1_3b_3a))
f

#bind d stats and f
res$f = f
res$pop = pop

write.table(res,file=paste0('Jackknife_Dstatistics_',pop,'2024JAN12.txt'),quote=F,sep='\t',row.names=F)


#### Plot ABBA-BABA D-statistics 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/ABBA/2024JAN')
library(tidyverse)
library(viridis)
library(data.table)
source("~/modules/genomics_general/jackknife.R")

res = read.table('~/merondun/cuculus_plumage/admixture/Jackknife_Dstatistics_2024JAN17.txt',header=TRUE)

ap = res %>% 
  filter(pop == 'HUN_Parity') %>% 
  ggplot(aes(x=Group,y=Djack.mean,ymin=Djack.mean-Djack.standard_deviation,
             ymax=Djack.mean+Djack.standard_deviation,
             col=Group,label=paste0(signif(Djack.mean,3),'\n (Z = ',signif(Z,3),')')))+
  geom_point(position=position_dodge(width=0.75))+
  geom_label(aes(y=Djack.mean+Djack.standard_deviation + 0.03),position=position_dodge(width=0.75),size=1.75)+
  geom_errorbar(width=0.25,position=position_dodge(width=0.75))+
  xlab('')+ylab('D-Statistic')+
  scale_color_manual(values=c('black','orange'))+
  geom_hline(yintercept=0,lty=2)+
  ylim(c(-0.06,0.14))+
  theme_bw()
ap

pdf('~/merondun/cuculus_plumage/admixture/Dstatistics_2024JAN17.pdf',height=2,width=2.5)
ap
dev.off()
