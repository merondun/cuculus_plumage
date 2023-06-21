#### Plot fst and dxy between groups, show variation in window size, plot heatmaps, and plot genome-wide variation 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/divergence/out')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(karyoploteR)
library(viridis)
library(LICORS)
library(ggplot2)
library(openxlsx)
library(scales)
library(ggpubr)
library(RColorBrewer)
library(gmodels)

# Load data and rename the columns
fst = read.table('Sensitivity_DXY.txt',header=TRUE)
names(fst) = c('chr','start','end','mid','sites','piA','piB','dxy','Fst','Group','Iteration', 'Window','Mask')

# Remove rows that have missing values (NA) in 'dxy', 'piA', 'piB' columns
fstd = fst %>% drop_na(dxy,piA,piB)

# Replace missing values in the 'Fst' column with 0, as these indicate no diversity at the site
fstd = fstd %>% replace_na(list(Fst=0))

# Threshold 'Fst' values to be between 0 and 1
fstd = fstd %>% mutate(Fst = threshold(Fst,min=0,max=1))

# Split the 'Group' column into two separate columns 'GroupA' and 'GroupB' based on underscore
fs = fstd %>%
  separate(Group,into=c('GroupA','GroupB'),remove=FALSE,sep='_')

# Sort the groups alphabetically and join them back with underscore
fs1 = fs %>%  rowwise() %>% mutate(pair = sort(c(GroupA,GroupB)) %>% paste(collapse = "_"))

# Drop the old 'GroupA' and 'GroupB' columns and rename the 'pair' column to 'Group', then split the 'Group' column again
fs2 = fs1 %>% select(!c(GroupA,GroupB,Group)) %>% dplyr::rename(Group = pair) %>% separate(Group,into=c('GroupA','GroupB'),remove=FALSE,sep='_')

# Add a new column 'AvZ' to classify chromosomes into 'W', 'Z', 'MT' and 'Autosome'
d1 = fs2 %>%
  mutate(AvZ = ifelse(chr == 'chr_W','W',ifelse(chr == 'chr_Z','Z',ifelse(chr == 'chr_MT','MT','Autosome'))))

# Exclude non-masked versions of 'chr_Z' and 'chr_W' from the data
d2 = d1 %>% filter(!((chr == "chr_Z" | chr == "chr_W") & (Mask == "Subset" | Mask == "All")))

# Clean up the 'Mask' column values
d3 = d2 %>% mutate(Mask = gsub('SubsetMask','Subset',Mask),
                   Mask = gsub('AllMask','All',Mask))

# Retain only one iteration for the 'All' mask
d4 = d3 %>% filter(!(Mask == "All" & Iteration != "1"))

# Write out the cleaned and manipulated data to a new file
write.table(d4,file='Sensitivity_DXY.input',quote=F,sep='\t',row.names=F)

# Read back the cleaned data
d4 = read.table('Sensitivity_DXY.input',header=TRUE)

# Reorder the levels of 'Group' factor
d4$Group = factor(d4$Group,levels=c('HG_HH','CG_CH','OG_OH','CG_OG','CH_OH'))

# Calculate mean, lower and upper confidence intervals for 'Fst' and 'dxy'
data = d4 %>% select(Group,Iteration,Window,Mask,AvZ,dxy,Fst) %>%
  pivot_longer(!c(Group,Iteration,Window,Mask,AvZ)) %>%
  group_by(Group,Iteration,Window,Mask,AvZ,name) %>%
  summarize(mean = ci(value)[1],
            lo = ci(value)[2],
            hi = ci(value)[3])

# Generate plot for 'Fst'
fstsen = data %>% filter(name == 'Fst') %>%
  ggplot(aes(x=AvZ,y=mean,ymin=lo,ymax=hi,col=Group))+
  geom_point(position=position_dodge(width=0.5))+
  geom_errorbar(position=position_dodge(width=0.5))+
  scale_color_viridis(discrete=TRUE,option='turbo')+
  facet_grid(Window~Mask,scales='free')+ggtitle('FST')+
  theme_bw()

# Generate plot for 'dxy'
dxysen = data %>% filter(name == 'dxy') %>%
  ggplot(aes(x=AvZ,y=mean,ymin=lo,ymax=hi,col=Group))+
  geom_point(position=position_dodge(width=0.5))+
  geom_errorbar(position=position_dodge(width=0.5))+
  scale_color_viridis(discrete=TRUE,option='turbo')+
  facet_grid(Window~Mask,scales='free')+ggtitle('DXY')+
  theme_bw()

# Save the plots into a pdf file
pdf('../Sensitivity_Window-FST-DXY.pdf',height=8,width=6)
ggarrange(fstsen,dxysen,common.legend = TRUE,nrow=2,ncol=1)
dev.off()

#Plot Overall estimates, so now ignore the iterations and just summarize all 
wins = c(5000,50000,500000)
vars = c('Fst','dxy')
overall = d4 %>% 
  filter(Mask == 'All') %>% 
  select(Group,Window,AvZ,dxy,Fst) %>% 
  pivot_longer(!c(Group,Window,AvZ)) %>% 
  group_by(Group,Window,AvZ,name) %>% 
  summarize(mean = mean(value,na.rm=TRUE),
            sd = sd(value,na.rm=TRUE))
for (win in wins) {
  for (va in vars) {
    cat('Plotting: ',va,' in ',win,'kb windows \n')
    #plot
    pp = ggplot(overall %>% filter(name == va & Window == win), aes(x = Group, y = AvZ, fill = mean,label=paste0(signif(mean,2),' ± ',signif(sd,2)))) +
      geom_tile() +
      geom_text(size=2)+ ylab('')+xlab('')+
      scale_fill_continuous(va,low='yellow',high='red')+
      theme_bw() + ggtitle(paste0(va,': ',win/1000,'KB windows')) +
      theme(strip.background = element_blank(),
            strip.text = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position='bottom')
    assign(paste0(va,win),pp)
  }
}
pdf('../Sensitivity_Window-FST-DXY_heatmap.pdf',height=12,width=10)
ggarrange(dxy5000,Fst5000,dxy50000,Fst50000,dxy500000,Fst500000,nrow=3,ncol=2)
dev.off()

# Set up genome
genome <- read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa.bed',header=FALSE) %>% filter(str_detect(V1,'scaff|MT',negate=T)) %>% arrange(desc(V2))
names(genome) <- c('chr','start','end')
genome$chr <- gsub('chr_','',genome$chr)

# Only chromosomes > 2mb 
genome = subset(genome,end > 2000000)
cuckooG <- makeGRangesFromDataFrame(genome,seqnames.field = 'chr',start.field = 'start',end.field = 'end')

# For sex chromosome-only view make mtDNA larger for visualization 
wgenome = rbind(genome %>% filter(chr == 'W' | chr == 'Z'),data.frame(chr = 'MT',start=-2000000,end=20000))
WG = makeGRangesFromDataFrame(wgenome,seqnames.field = 'chr',start.field = 'start',end.field = 'end')

# Stagger chromosome names 
evens <- ((1:length(genome$chr))%%2)==0
chr.names <- genome$chr
even.names <- chr.names
even.names[!evens] <- ""
odd.names <- chr.names
odd.names[evens] <- ""
odd.names[length(odd.names) - 1] <- "W"
odd.names = gsub('29','',odd.names) #micros are too busy, remove 
odd.names = gsub('29','',odd.names) #micros are too busy, remove 
odd.names = gsub('25','',odd.names) #micros are too busy, remove 

# Only plot 50kb and all samples 
divs = d4 %>% filter(Window == 50000 & Mask == 'All') %>% mutate(chr = gsub('chr_','',chr))

#subset each raw dataset into specific colors
gb1 = genome[seq(1, nrow(genome), 2), ][[1]]
gb2 = genome[seq(2, nrow(genome), 2), ][[1]]
div_gb1 <- divs[grepl(paste(gb1, collapse="$|"), divs$chr), ]
div_gb1$Color <- "grey60"
div_gb2 <- divs[grepl(paste(gb2, collapse="$|"), divs$chr), ]
div_gb2$Color <- "grey10"
divs_use <- rbind(div_gb1,div_gb2)

# Define the variable to examine 
var = 'dxy'
#pdf('../Karyoplot_Legend.pdf',height=5,width=9)
#pdf('../DXY_Genome-wide.pdf',height=5,width=9)

# Use this section for whole genome
pp <- getDefaultPlotParams(plot.type=4)
pp$leftmargin <- 0.1
png(paste0('../',var,'_Genome-wide_50kb.png'),units='in',res=600,height=5,width=9)
kp <- plotKaryotype(plot.type=4, genome = cuckooG,plot.params = pp,labels.plotter = NULL)
kpAddChromosomeNames(kp, chr.names = odd.names,cex=0.8)
kpAddBaseNumbers(kp,tick.dist = 25000000,cex=0.4)

# Use this section for sex chromosomes 
# png(paste0('../',var,'_Sex_50kb.png'),units='in',res=600,height=5,width=9)
# kp <- plotKaryotype(plot.type=4, genome = WG)
# kpAddBaseNumbers(kp,tick.dist = 25000000,cex=0.4)
counter = 0 
tracks = 5

# Loop through all the groups, plot the variable indicated above, add 1 to the counter so it gets added to the next panel 
vars = c('HG_HH','CG_CH','OG_OH','CG_OG','CH_OH')
for (varz in vars) {
  if (var == 'dxy') {
    cat('Working with dxy: ',varz,'\n')
    dat = divs_use %>% filter(Group == varz) %>% mutate(chr = gsub('chr_','',chr)) %>% dplyr::rename(value=dxy)
    lower = 0
    upper = 0.02
  } else if (var == 'maxAF') {
    cat('Working with dAF: ',varz,'\n')
    dat = divs_use %>% filter(Group == varz) %>% mutate(chr = gsub('chr_','',chr)) %>% dplyr::rename(value=maxAF)
    lower = 0
    upper = 1
  } else {
    cat('Working with fst: ',varz,'\n')
    dat = divs_use %>% filter(Group == varz) %>% mutate(chr = gsub('chr_','',chr)) %>% dplyr::rename(value=Fst)
    lower = 0
    upper = 1
  }
  counter <- counter+1
  cat('Current Autotrack: ',counter,' for ',varz,'\n')
  at <- autotrack(current.track = counter, total.tracks = tracks);
  at$r1 <- at$r1-0.03

  #add lines at 0,1 and 0.5
  kpAbline(kp,h=c(0,1),lty=1,r0=at$r0,r1=at$r1,cex=0.6,col='grey60')
  kpAbline(kp,h=0.5,lty=2,r0=at$r0,r1=at$r1,cex=0.6,col='grey60')

  #add points 
  kpPoints(kp,chr=dat$chr,x=dat$start,y=dat$value,ymin=lower,ymax=upper,r0=at$r0,r1=at$r1,col=dat$Color)
  #kpLines(kp,chr=dat$chr,x=dat$start,y=dat$dxy,ymin=min(dat$dxy),ymax=max(dat$dxy),r0=at$r0,r1=at$r1,alpha=0.5,col=dat$Color)

  # Add labels
  kpAxis(kp,ymin=lower,ymax=upper,r0=at$r0,r1=at$r1,cex=0.6)
  kpAddLabels(kp,cex=0.8,labels = varz,r0=at$r0+.01, r1=at$r1,col="black",srt=0,label.margin = 0.05)
}
dev.off()

