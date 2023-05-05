#### Full ml tree
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
chrs = 'chr_W'
#Symbology bars
counter = 0 
for (chr in chrs) {
  files = list.files('.',paste0(chr, '.*WithChinaMask.*contree'))
  for (tree in files) {
    cat('Making tree for ',tree,'\n')
    counter = counter + 1 
    iqtree = read.iqtree(tree) 
    label = gsub('.min4.*','',tree)
    #iqtree = read.iqtree(list.files('.',paste0(chr, '.NoChina\\..*contree')))
    gg <- ggtree(iqtree, layout = "dendrogram") %<+% md
    gg$data$label = ifelse(gg$data$isTip == TRUE,paste0(gg$data$IDNumber,'_',gg$data$Species,'_',gg$data$Plumage),gg$data$label)
    
    m <- MRCA(gg, '387_CP_mixed')
    y <- groupClade(gg$data, m)
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
      #geom_nodelab(mapping=aes(subset=(GeneConcordance >= 25),label=GeneConcordance),col='black',size=2,show.legend=F,vjust=-1,hjust=2)+
      geom_nodepoint(mapping=aes(subset=(SH_aLRT >= 80 & UFboot > 95)),col='black',pch=9,size=2,show.legend=F)+
      #scale_x_continuous(expand = expansion(mult=0.05)) +
      ggtitle(label)+
      #geom_treescale(y = 50) +
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
ord = get_taxa_name(p) %>% as.data.frame()
ord = ord %>% mutate(TreeOrder = row_number())
names(ord) = c('ID','TreeOrder')
write.table(ord,file='/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/mltree/Chr_W_TreeWithChinaMask.order',quote=F,sep='\t',row.names=F)

pdf('chr_W_Tree_WithChinaMask_Halloween.pdf',height=5,width=8)
p
dev.off()

