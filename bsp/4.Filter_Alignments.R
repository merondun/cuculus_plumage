### Identify reliable sites 
.libPaths('~/mambaforge/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/bsp/2023JUNE_GENES')
library(ggpubr)
library(tidyverse)
library(viridis)
library(LICORS)
library(karyoploteR)
library(zoo)
library(RColorBrewer)

p = read.table('output.tsv',header=TRUE)
p = p %>% filter(sequence_id != 'Reference') %>% dplyr::rename(ID=sequence_id)
p$file_name = gsub('.fa','',p$file_name)
md = read.table('~/merondun/cuculus_plumage/EntireMetadata.txt',header=TRUE,comment.char = '',sep='\t')
pt = left_join(p,md) #ensure that the 'ID' from the vcf matches the 'ID' from your metadata
sample_pi = pt %>% 
  mutate(ID = fct_relevel(ID,Plumage)) %>% 
  ggplot(aes(y=ID,x=bs_distance,fill=Plumage))+
  geom_boxplot(show.legend=FALSE)+
  facet_grid(Group~.,scales='free')+
  theme_bw()+
  scale_fill_viridis(discrete=TRUE)

length_pi = pt %>% 
  select(file_name,num_polymorphic,sequence_length) %>%
  unique %>% 
  ggplot(aes(x=num_polymorphic,y=sequence_length))+
  geom_point()+
  theme_bw()

#Filter genes without any polymorphic sites 
lowgenes = pt %>% filter(num_polymorphic < 10) %>% select(file_name) %>% unique
pf = pt %>% filter(!file_name %in% lowgenes$file_name)
length(unique(pf$file_name))
pf %>% 
  select(file_name,num_polymorphic,sequence_length) %>%
  unique %>% 
  ggplot(aes(x=num_polymorphic,y=sequence_length))+
  geom_point()+
  theme_bw()

#what total sequence is left
pf %>% 
  select(file_name,num_polymorphic,sequence_length) %>%
  unique %>% summarize(total_length = sum(sequence_length))
#still too large, 3.8MB, let's remove some of the largest and smallest genes 

#Filter very large genes and very small genes 
edgegenes = pf %>% filter(sequence_length > 75000 | sequence_length < 10000) %>% select(file_name) %>% unique
pf2 = pf %>% filter(!file_name %in% edgegenes$file_name)
length(unique(pf2$file_name))
pf2 %>% 
  select(file_name,num_polymorphic,sequence_length) %>%
  unique %>% summarize(total_length = sum(sequence_length))
pf2 %>% 
  select(file_name,num_polymorphic,sequence_length) %>%
  unique %>% 
  ggplot(aes(x=num_polymorphic,y=sequence_length))+
  geom_point()+
  theme_bw()

#1474285, let's try this 
write.table(unique(pf2$file_name),file='Retained_Genes_Between10k-75k_BSP.txt',quote=F,sep='\t',row.names=F,col.names=F)

