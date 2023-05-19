setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/withmic/mltree/')
.libPaths('~/mambaforge/envs/mtdna/lib/R/library')
library(ape)
library(tidyverse)

#Read in metadata
md = read.table('../EntireMetadata.txt',header=TRUE,comment.char='',sep='\t')
md = md %>% mutate(PlumageColor = gsub("#F7F7F7",'grey80',PlumageColor))

tree1 <- read.tree(list.files('.',paste0('chr_W.*__SetB.*contree')))
tree2 <- read.tree(list.files('.',paste0('chr_MT.*SetB.*contree')))
tree1$data = left_join(data.frame(ID=tree1$tip.label),md )
tree1$data = tree1$data %>% mutate(PID = paste0(IDNumber,'_',Species,'_',Plumage))
tree1$tip.label = tree1$data$PID
tree2$data = left_join(data.frame(ID=tree2$tip.label),md )
tree2$data = tree2$data %>% mutate(PID = paste0(IDNumber,'_',Species,'_',Plumage))
tree2$tip.label = tree2$data$PID

# create associations matrix 
a <- as.character(c(tree2$tip.label))
b <- as.character(c(tree2$tip.label)) 
association <- cbind(a, b)

pdf('Cophyloplot_chr_W-MT_SETB.pdf',height=10,width=20)
cophyloplot(tree1, tree2, assoc = association,length.line = 25, space = 150, gap = 5)
dev.off()