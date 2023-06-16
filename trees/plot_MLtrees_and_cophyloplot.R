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

chrs = c('chr_W','chr_MT','chr_Z','chr_6')
#chrs = 'chr_W'
#Symbology bars
counter = 0 
for (chr in chrs) {
  files = list.files('max-missing-05_SetB',paste0(chr, '.*contree'),full.names = TRUE)
  for (tree in files) {
    cat('Making tree for ',tree,'\n')
    counter = counter + 1 
    iqtree = read.iqtree(tree) 
    iqtree2 = phytools::reroot(as.phylo(iqtree),interactive=TRUE)
    iqtree3 = drop.tip(iqtree2, '327_CO_SCC_RUS_F')
    label = gsub('.min4.*','',tree)
    gg <- ggtree(iqtree3, layout = "dendrogram") %<+% md
    gg$data = gg$data %>% mutate(Plumage = ifelse(Species == 'CP' | Species == 'CM','OUT',Plumage),
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
      #geom_tiplab(size=2,aes(col=Plumage,lwd=Species),align = TRUE,angle=90,hjust=1)+
      scale_fill_manual(values=gg$data$PlumageColor,breaks=gg$data$Plumage)+
      scale_color_manual(values=gg$data$PlumageColor,breaks=gg$data$Plumage)+
      scale_shape_manual(values=gg$data$Shape,breaks=gg$data$Species_Latin)+
      geom_nodepoint(mapping=aes(subset=(SHalrt > 99)),col='black',pch=9,size=2,show.legend=F)+
      ggtitle(chr)+
      theme(legend.position = 'none',plot.title = element_text(hjust = 0.5))
    p
    #shorten branch length
    p$data[p$data$node %in% c(m), "x"] <- mean(p$data$x)
    #p = p + geom_treescale(x = -0.5, y = 15,linesize = 1)
    p = p + geom_treescale()
    p
    assign(paste0('p',counter),p)
    
  }
}
p1
p2
p3
p4
#chrW
ord = get_taxa_name(p1) %>% as.data.frame()
ord = ord %>% mutate(TreeOrder = row_number())
names(ord) = c('ID','TreeOrder')
write.table(ord,file='/dss/dsshome1/lxc07/di39dux/merondun/cuculus_plumage/trees/Chr_W_TreeWithMIC-NoChina-MM05-SetC_ROOT.order',quote=F,sep='\t',row.names=F)

#chrMT
ord = get_taxa_name(p2) %>% as.data.frame()
ord = ord %>% mutate(TreeOrder = row_number())
names(ord) = c('ID','TreeOrder')
write.table(ord,file='/dss/dsshome1/lxc07/di39dux/merondun/cuculus_plumage/trees/Chr_MT_TreeWithMIC-NoChina-MM05-SetC_ROOT.order',quote=F,sep='\t',row.names=F)

pdf('chr_W_WithMIC-NoChina-MM05-SetC_ROOT_Halloween.pdf',height=4,width=8)
p1
dev.off()

pdf('chr_MT_WithMIC-NoChina-MM05-SetC_ROOT_Halloween.pdf',height=4,width=8)
p2
dev.off()

pdf('AllChromosomes-MM05-SetB_ROOT_Halloween.pdf',height=6,width=9)
ggarrange(p1,p2,p3,p4,nrow=2,ncol=2)
dev.off()

#cophyloplot - plot the trees facing each other with links 
files = list.files('.',paste0('.*_WithMicropterus-NoChina-Mask-MM05__SetC.*contree'))
tree1 <- read.tree(files[2])
tree2 <- read.tree(files[1])
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

pdf('Cophyloplot_chr_W-MT_WithMic-NoChina-MM05_SETC.pdf',height=10,width=20)
cophyloplot(tree1, tree2, assoc = association,length.line = 25, space = 150, gap = 5)
dev.off()

