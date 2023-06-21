# Set the working directory and load necessary R libraries
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

# Read in metadata file and replace specific color in 'PlumageColor' column
md = read.table('../EntireMetadata.txt',header=TRUE,comment.char='',sep='\t')
md = md %>% mutate(PlumageColor = gsub("#F7F7F7",'grey80',PlumageColor))

# Define a vector of chromosomes for the subsequent analyses
chrs = c('chr_W','chr_MT','chr_Z','chr_6')

# Initialize a counter for the loop
counter = 0

# Loop over each chromosome
for (chr in chrs) {
  # Get list of files for the current chromosome
  files = list.files('max-missing-05_SetB',paste0(chr, '.*contree'),full.names = TRUE)
  
  # Loop over each file (tree)
  for (tree in files) {
	# Print the current tree being processed
	cat('Making tree for ',tree,'\n')

	# Increment the counter by 1
	counter = counter + 1

	# Read the tree using iqtree function from treeio package
	iqtree = read.iqtree(tree)

	# Reroot the tree interactively using the phytools package
	iqtree2 = phytools::reroot(as.phylo(iqtree),interactive=TRUE)

	# Remove a specific tip (leaf) from the tree 
	iqtree3 = drop.tip(iqtree2, '327_CO_SCC_RUS_F')

	# Clean up the tree label
	label = gsub('.min4.*','',tree)

	# Create a ggtree object from the iqtree3 object and add metadata to it
	gg <- ggtree(iqtree3, layout = "dendrogram") %<+% md

	# Modify the ggtree data to change 'Plumage' and 'PlumageColor' values for certain species
	gg$data = gg$data %>% mutate(Plumage = ifelse(Species == 'CP' | Species == 'CM','OUT',Plumage),
	                             PlumageColor = ifelse(Species == 'CP' | Species == 'CM','grey60',PlumageColor))

	# Modify the label column for tip nodes
	gg$data$label = ifelse(gg$data$isTip == TRUE,paste0(gg$data$IDNumber,'_',gg$data$Species,'_',gg$data$Plumage),gg$data$label)

	# Find the most recent common ancestor (MRCA) of '387_CP_OUT'
	m <- MRCA(gg, '387_CP_OUT')

	# Group the clades based on the MRCA
	y <- groupClade(gg$data, m)

	# Convert SHalrt labels to numeric
	y$SHalrt = as.numeric(y$label)

	# Create the ggtree plot with specific aesthetics and configurations
	p = ggtree(y, aes(linetype = group), layout='dendrogram')+
	    geom_tippoint(aes(fill = Plumage,shape=Species_Latin),size=3)+
	    scale_fill_manual(values=gg$data$PlumageColor,breaks=gg$data$Plumage)+
	    scale_color_manual(values=gg$data$PlumageColor,breaks=gg$data$Plumage)+
	    scale_shape_manual(values=gg$data$Shape,breaks=gg$data$Species_Latin)+
	    geom_nodepoint(mapping=aes(subset=(SHalrt > 99)),col='black',pch=9,size=2,show.legend=F)+
	    ggtitle(chr)+
	    theme(legend.position = 'none',plot.title = element_text(hjust = 0.5))

	# Shorten branch length for the node identified as MRCA
	p$data[p$data$node %in% c(m), "x"] <- mean(p$data$x)

	# Add a tree scale to the plot
	p = p + geom_treescale()

	# Assign the plot to a dynamically generated variable name (e.g., p1, p2, etc.)
	assign(paste0('p',counter),p)

  }
}

p1
p2
p3
p4

# Save the ordering of tip labels on each tree into separate files: chrW
ord = get_taxa_name(p1) %>% as.data.frame()
ord = ord %>% mutate(TreeOrder = row_number())
names(ord) = c('ID','TreeOrder')
write.table(ord,file='/dss/dsshome1/lxc07/di39dux/merondun/cuculus_plumage/trees/Chr_W_TreeWithMIC-NoChina-MM05-SetC_ROOT.order',quote=F,sep='\t',row.names=F)

# Save the ordering of tip labels on each tree into separate files: mtDNA
ord = get_taxa_name(p2) %>% as.data.frame()
ord = ord %>% mutate(TreeOrder = row_number())
names(ord) = c('ID','TreeOrder')
write.table(ord,file='/dss/dsshome1/lxc07/di39dux/merondun/cuculus_plumage/trees/Chr_MT_TreeWithMIC-NoChina-MM05-SetC_ROOT.order',quote=F,sep='\t',row.names=F)

# Save the plot for 'chr_W' as a PDF
pdf('chr_W_WithMIC-NoChina-MM05-SetC_ROOT_Halloween.pdf',height=4,width=8)
# Display the plot for 'chr_W'
p1
# Turn off the current device driver (closes the PDF file)
dev.off()

# Save the plot for 'chr_MT' as a PDF
pdf('chr_MT_WithMIC-NoChina-MM05-SetC_ROOT_Halloween.pdf',height=4,width=8)
p2
dev.off()

# Save the combined plot of all chromosomes as a PDF
pdf('AllChromosomes-MM05-SetB_ROOT_Halloween.pdf',height=6,width=9)
ggarrange(p1,p2,p3,p4,nrow=2,ncol=2)
dev.off()

# Cophyloplot creation
# Get list of all  tree files in the current directory
files = list.files('.',paste0('.*_WithMicropterus-NoChina-Mask-MM05__SetC.*contree'))
# Read in the second and first trees from the files
tree1 <- read.tree(files[2])
tree2 <- read.tree(files[1])

# Join the metadata with the trees, create PID from IDNumber, Species and Plumage, and update tip.label with PID
tree1$data = left_join(data.frame(ID=tree1$tip.label),md )
tree1$data = tree1$data %>% mutate(PID = paste0(IDNumber,'_',Species,'_',Plumage))
tree1$tip.label = tree1$data$PID

tree2$data = left_join(data.frame(ID=tree2$tip.label),md )
tree2$data = tree2$data %>% mutate(PID = paste0(IDNumber,'_',Species,'_',Plumage))
tree2$tip.label = tree2$data$PID

# Create an associations matrix
a <- as.character(c(tree2$tip.label))
b <- as.character(c(tree2$tip.label))
association <- cbind(a, b)

# Save the cophyloplot as a PDF
pdf('Cophyloplot_chr_W-MT_WithMic-NoChina-MM05_SETC.pdf',height=10,width=20)
# Generate a cophyloplot (trees facing each other with links)
cophyloplot(tree1, tree2, assoc = association,length.line = 25, space = 150, gap = 5)
dev.off()

