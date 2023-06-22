#### Full ml tree
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/bsp')
.libPaths('~/mambaforge/envs/mtdna/lib/R/library')
library(rhierbaps)
library(ggtree)
library(phytools)
library(ape)
library(openxlsx)
library(treeio)
library(ggnewscale)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

#Read in metadata
md = read.table('~/merondun/cuculus_plumage/EntireMetadata.txt',header=TRUE,comment.char='',sep='\t')
md = md %>% mutate(PlumageColor = gsub("#F7F7F7",'grey80',PlumageColor))
files = list.files('.','*contree')
iqtree = read.iqtree(files) 
gg <- ggtree(iqtree, layout = "ape") %<+% md
p = gg + 
      geom_tippoint(aes(fill = Plumage,shape=Species_Latin),size=3)+
      scale_fill_manual(values=gg$data$PlumageColor,breaks=gg$data$Plumage)+
      scale_color_manual(values=gg$data$PlumageColor,breaks=gg$data$Plumage)+
      scale_shape_manual(values=gg$data$Shape,breaks=gg$data$Species_Latin)+
      geom_nodepoint(mapping=aes(subset=(SH_aLRT > 99)),col='black',pch=9,size=2,show.legend=F)+
      theme(legend.position = 'none')
p
 
pdf('BSP_Tree.pdf',height=6,width=7)
p
dev.off()
