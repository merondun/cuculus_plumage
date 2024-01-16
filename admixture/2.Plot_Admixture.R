### Plot ADMIXTURE
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/admixture')
library(tidyverse)
library(viridis)

#First, read in q matrices 
RUN='Autosomes.IF-GF-MM2-BP-ANN-AC2.TransSpecies'
qs = list.files('.',paste0(RUN,'.*.Q$')) #find all files ending in .Q
famfile = paste0('../autosomes/',RUN,'.fam') #specify the .fam file 
samps = read.table(famfile,header=FALSE) %>% select(V2) %>% dplyr::rename(ID=V2) %>% 
  mutate(ID = gsub('_F_.*','_F',ID)) #with plink you will projbably get a double ID, I use this to strip the second ID, double check
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
qdat = left_join(qdat,md) %>% mutate(IDtree = paste0(IDNumber,'_',Species,'_',Plumage))
#save to github for plotting 
write.table(qdat,file='/dss/dsshome1/lxc07/di39dux/merondun/cuculus_plumage/admixture/Qcoefficients.txt')

####Skip to hear if starting from github
setwd('/dss/dsshome1/lxc07/di39dux/merondun/cuculus_plumage/admixture/')
qdat = read.table('Qcoefficients.txt',header=TRUE)
#I have a specified order already since I want to align samples in the same order as a phylogenetic tree, this is just a vector of IDs 
ord = read.table('../trees/max-missing-20_SetB/Chr_W_TreeWithMIC-SetB_ROOT.order',header=TRUE)
qdat = left_join(ord %>% dplyr::rename(IDtree = ID),qdat) %>% drop_na(Q)
qdat$IDtree = factor(qdat$IDtree,levels=ord$ID)

#without order
#qdat = qdat %>% mutate(IDtree = fct_reorder(ID,desc(Plumage)))

#Prep
qdat = qdat %>% mutate(Kclust = paste0('K',MaxK))
qdat$Kclust <- factor(qdat$Kclust, levels = paste0("K", 2:10))

#summarize canorus and their optatus within their respective K2 clusters
qdat %>% filter(Kclust == 'K2' & K == 'K1' & Species == 'CC') %>% summarize(mean = mean(Q),min=min(Q),max=max(Q))
qdat %>% filter(Kclust == 'K2' & K == 'K2' & Species == 'CO') %>% summarize(mean = mean(Q),min=min(Q),max=max(Q))

#assign each K to a species, and then calculate the maximum Q observed from the alternate species at that that K , for all Ks 
qd = qdat %>% mutate(Species = ifelse(grepl('CC',Species),'C. canorus','C. optatus'))
max_species_per_K <- qd %>%
  group_by(Kclust, K) %>%
  filter(Q == max(Q)) 
max_other_species_per_K <- qd %>%
  anti_join(max_species_per_K, by = c("Kclust", "K", "Species")) %>%
  group_by(Kclust, K) %>%
  filter(Q == max(Q))
kvar = ggplot() +
  geom_point(data = max_species_per_K, aes(y = K, x = Q, color = Species)) +
  geom_point(data = max_other_species_per_K, aes(y = K, x = Q, color = Species)) +
  facet_grid( Kclust ~ .,scales='free',space='free') +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme_bw()+xlab('Maximum Observed Ancestry Coefficient (Q)')+
  theme(legend.position='top')

png('Admixture_WithMIC-SetC_MM2-KInterspecies.png',height=8,width=7,units='in',res=600)
kvar
dev.off()

#Now plot
k2plot =
  qdat %>% filter(MaxK < 11) %>%  #I only want to show K2-4
  mutate(Group = ifelse(grepl('CG',Group),'C. canorus grey',ifelse(grepl('CH',Group),'C. canorus rufous',
                                                                   ifelse(grepl('OG',Group),'C. optatus grey',ifelse(grepl('OH',Group),'C. optatus rufous','NO'))))) %>% 
  ggplot(aes(x = factor(IDtree), y = Q, fill = factor(K))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(Kclust~ Group, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "",y = "Ancestry Coefficient") +
  scale_y_continuous(n.breaks = 3) +
  scale_fill_viridis('K',discrete=TRUE,option='turbo')+ 
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    #axis.text.x = element_blank(),
    axis.text.x = element_text(size=5,angle=90,vjust=0,hjust=0),
    panel.grid = element_blank(),
    strip.text = element_text(size = 5),
    legend.position='top')
k2plot

png('Admixture_TransSpecies-K10-Facet.png',height=6,width=8,units='in',res=600)
k2plot
dev.off()

#Now plot main figure 
fig3b =
  qdat %>% filter(MaxK < 6) %>%  #I only want to show K2-4
  ggplot(aes(x = factor(IDtree), y = Q, fill = factor(K))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(Kclust~., switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "",y = "Ancestry Coefficient") +
  scale_y_continuous(n.breaks = 3) +
  scale_fill_viridis('K',discrete=TRUE)+ 
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    #axis.text.x = element_blank(),
    axis.text.x = element_text(size=5,angle=90,vjust=0,hjust=0),
    panel.grid = element_blank(),
    strip.text = element_text(size = 5),
    legend.position='top')
fig3b

pdf('Admixture_TransSpecies-K6.pdf',height=3,width=6)
fig3b
dev.off()