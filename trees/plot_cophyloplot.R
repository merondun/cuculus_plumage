library(ape)
tree1 <- read.tree("chr_W_WithChinaMask.contree")
tree2 <- read.tree("chr_MT_WithChinaMask.contree")
tree1$data = left_join(data.frame(ID=tree1$tip.label),md )
tree1$data = tree1$data %>% mutate(PID = paste0(IDNumber,'_',Species,'_',Plumage))
tree1$tip.label = tree1$data$PID
tree2$data = left_join(data.frame(ID=tree2$tip.label),md )
tree2$data = tree2$data %>% mutate(PID = paste0(IDNumber,'_',Species,'_',Plumage))
tree2$tip.label = tree2$data$PID

# create associations matrix 
a <- as.character(c(tree1$tip.label))
b <- as.character(c(tree1$tip.label)) 
association <- cbind(a, b)

pdf('Cophyloplot_chr_W-MT.pdf',height=10,width=20)
cophyloplot(tree1, tree2, assoc = association,length.line = 2, space = 75, gap = 7)
dev.off()

