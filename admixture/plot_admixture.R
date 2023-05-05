### Plot ADMIXTURE
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/admixture')
library(tidyverse)
library(viridis)

#RUN='Autosomes.IF-GF-MM2-BP-ANN-AC2-LD.HUNGARY'
RUN='Autosomes.IF-GF-MM2-BP-ANN-AC2-LDSTRICT.AnalysisBIngroup'
qs = list.files('.',paste0(RUN,'.*.Q$')) #find all files ending in .Q
famfile = paste0('../autosomes/',RUN,'.fam') #specify the .fam file 
samps = read.table(famfile,header=FALSE) %>% select(V2) %>% dplyr::rename(ID=V2)
qdat = NULL
for (q in qs){
  cat('Melting K: ',q,'\n')
  qf = read.table(q) #read in file
  k = ncol(qf) #k is simply the total number of columns
  names(qf) = paste0('K',seq(1,k)) #add 'KX' for all columns
  qfm = cbind(samps,qf) #bind with samples, ensure you're using the FAM file used for admixture! 
  qfm$MaxK = k   
  meltq = qfm %>% pivot_longer(!c(ID,MaxK),names_to = 'K',values_to = 'Q')  #melt it 
  qdat = rbind(qdat,meltq)
}

#read in metadata 
md = read.table('../EntireMetadata.txt',header=TRUE,comment.char = '',sep='\t')
#only keep data analysis set A+B. Since I align samples with the tree order, the left join with ord will remove individuals not in the tree anyways!
qdat = left_join(qdat,md) %>% mutate(IDtree = paste0(IDNumber,'_',Species,'_',Plumage))

#I have a specified order already since I want to align samples in the same order as a phylogenetic tree, this is just a vector of IDs 
ord = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/mltree/Chr_W_TreeNoChinaNoToesMask.order',header=TRUE)
qdat = left_join(ord %>% dplyr::rename(IDtree = ID),qdat)
qdat = qdat %>% mutate(IDtree = fct_reorder(IDtree,desc(TreeOrder)),
                       Kclust = paste0('K',MaxK))
qdat$Kclust <- factor(qdat$Kclust, levels = paste0("K", 2:10))

#without order
#qdat = qdat %>% mutate(IDtree = fct_reorder(ID,desc(Plumage)))

#Now plot
k2plot =
  qdat %>% filter(MaxK < 7) %>%  #I only want to show K2-4
  ggplot(aes(x = factor(IDtree), y = Q, fill = factor(K))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(Kclust~ ., switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "",y = "Ancestry Coefficient") +
  scale_y_continuous(n.breaks = 3) +
  scale_fill_viridis(discrete=TRUE,option='viridis')+ 
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    #axis.text.x = element_blank(),
    axis.text.x = element_text(angle=90),
    panel.grid = element_blank(),
    legend.position='bottom')
k2plot

pdf('Admixture_NoChinaNoToesMask_LDStrict-K6.pdf',height=4,width=8)
#png('Admixture_HUNGARY.png',height=2,width=7,units='in',res=600)
k2plot
dev.off()


