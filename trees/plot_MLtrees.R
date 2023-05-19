#### Full ml tree
setwd('/dss/dsshome1/lxc07/di39dux/merondun/cuculus_plumage/trees')
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
md = read.table('../EntireMetadata.txt',header=TRUE,comment.char='',sep='\t')
md = md %>% mutate(PlumageColor = gsub("#F7F7F7",'grey80',PlumageColor))

chrs = c('chr_W','chr_MT')
#chrs = 'chr_W'
#Symbology bars
counter = 0 
for (chr in chrs) {
  files = list.files('.',paste0(chr, '.*WithChinaWithMicropterusMask__SetB.*contree'))
  for (tree in files) {
    cat('Making tree for ',tree,'\n')
    counter = counter + 1 
    iqtree = read.iqtree(tree) 
    iqtree2 = phytools::reroot(as.phylo(iqtree),interactive=TRUE)
    iqtree3 = drop.tip(iqtree2, '327_CO_SCC_RUS_F')
    label = gsub('.min4.*','',tree)
    gg <- ggtree(iqtree3, layout = "dendrogram") %<+% md
    gg$data = gg$data %>% mutate(Plumage = ifelse(Species == 'CP' | (Species == 'CM' & is.na(Plumage)),'OUT',Plumage),
                                 PlumageColor = ifelse(Species == 'CP' | Species == 'CM','grey60',PlumageColor))
    gg$data$label = ifelse(gg$data$isTip == TRUE,paste0(gg$data$IDNumber,'_',gg$data$Species,'_',gg$data$Plumage),gg$data$label)
    #without short
    m <- MRCA(gg, '387_CP_OUT')
    y <- groupClade(gg$data, m)
    y$SHalrt = as.numeric(y$label)
    p = ggtree(
      #iqtree, 
      y,
      aes(linetype = group),
      layout='dendrogram')+ 
      geom_tippoint(aes(fill = Plumage,shape=Species_Latin),size=3)+
      geom_tiplab(size=2,aes(col=Plumage,lwd=Species),align = TRUE,angle=90,hjust=1)+
      scale_fill_manual(values=gg$data$PlumageColor,breaks=gg$data$Plumage)+
      scale_color_manual(values=gg$data$PlumageColor,breaks=gg$data$Plumage)+
      scale_shape_manual(values=gg$data$Shape,breaks=gg$data$Species_Latin)+
      geom_nodepoint(mapping=aes(subset=(SHalrt > 99)),col='black',pch=9,size=2,show.legend=F)+
      ggtitle(label)+
      theme(legend.position = 'none')
    p
    #shorten branch length
    p$data[p$data$node %in% c(m), "x"] <- mean(p$data$x)
    #p = p + geom_treescale(x = -0.5, y = 15,linesize = 1)
    p = p + geom_treescale(x = 0.15)
    p
    assign(paste0('p',counter),p)
    
  }
}
p1
p2
#chrW
ord = get_taxa_name(p1) %>% as.data.frame()
ord = ord %>% mutate(TreeOrder = row_number())
names(ord) = c('ID','TreeOrder')
write.table(ord,file='/dss/dsshome1/lxc07/di39dux/merondun/cuculus_plumage/trees/Chr_W_TreeWithMIC-SetB_ROOT.order',quote=F,sep='\t',row.names=F)

#chrMT
ord = get_taxa_name(p2) %>% as.data.frame()
ord = ord %>% mutate(TreeOrder = row_number())
names(ord) = c('ID','TreeOrder')
write.table(ord,file='/dss/dsshome1/lxc07/di39dux/merondun/cuculus_plumage/trees/Chr_MT_TreeWithMIC-SetB_ROOT.order',quote=F,sep='\t',row.names=F)

pdf('chr_W_WithMIC-SetB_ROOT_Halloween.pdf',height=4,width=8)
p1
dev.off()

pdf('chr_MT_WithMIC-SetB_ROOT_Halloween.pdf',height=4,width=8)
p2
dev.off()