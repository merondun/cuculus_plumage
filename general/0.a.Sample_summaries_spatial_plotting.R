.libPaths('~/mambaforge/envs/r/lib/R/library')
setwd('/dss/dsshome1/lxc07/di39dux/merondun/cuculus_plumage/general')
library(openxlsx)
library(tidyverse)
library(sf)
library(ggspatial)
library(spThin)
library(viridis)
library(factoextra)
library(ggpubr)
library(forcats)
library(RColorBrewer)
library(gghalves)

mdk = read.table('../EntireMetadata.txt',header=TRUE,comment.char='',sep='\t')
mdk = mdk %>% group_by(Species) %>% 
  mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
  { cbind(., pca = prcomp(.[, c("Longitude", "Latitude")], scale. = TRUE)$x[,1]) } %>%
  rename(Distance = pca) 
mdk = mdk %>% mutate(Species_Latin = gsub('C. canorus canorus','C. canorus',Species_Latin))
mdk$Species_Latin = factor(mdk$Species_Latin,levels=c('C. canorus','C. optatus','C. poliocephalus'))
mdk = mdk %>% dplyr::rename(Plumage = AdultPlumage)
md = mdk
q1 = md %>% 
  group_by(Species_Latin,Plumage) %>% count() %>% 
  ggplot(aes(x=Species_Latin,fill=Plumage,y=n,label=n))+
  geom_text(position=position_dodge2(width=1,preserve = 'single'),vjust=-.5)+
  theme_bw()+
  geom_bar(stat='identity',width=1,position=position_dodge2(width=1,preserve = 'single'),col='black')+
  xlab('')+ylab('Count')+
  scale_fill_manual(values=md$PlumageColor,breaks=md$Plumage)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))+
  scale_y_continuous(expand = expansion(0.1))
q1

#boxplot
q2 = md %>% 
  ggplot(aes(y=MedianCoverage, x=Plumage,fill=Plumage))  + 
  ggdist::stat_halfeye(adjust = .5,width = .6,.width = 0,justification = -.2, point_colour = NA,alpha = 0.5)+
  geom_boxplot(width = .15,outlier.shape = NA, alpha = 0.3) +
  gghalves::geom_half_point(aes(col=Plumage),side='l',range_scale = .4,alpha = .3)+
  ylab('Median Coverage')+xlab('')+
  scale_fill_manual(values=md$PlumageColor,breaks=md$Plumage)+
  scale_color_manual(values=md$PlumageColor,breaks=md$Plumage)+
  facet_grid(.~Species_Latin,scales='free')+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position = 'none')
q2

pdf('Data_Summary.pdf',height=4,width=9)
png('Data_Summary.png',units='in',res=600,height=4,width=9)
ggarrange(q1,q2,common.legend = TRUE,widths=c(0.3,0.7))
dev.off()

md %>% 
  ggplot(aes(x=Source,fill=Plumage))+
  geom_bar()+
  scale_fill_manual(values=md$PlumageColor,breaks=md$Plumage)+
  facet_grid(Species_Latin~.,scales='free')+
  theme_bw()


####### Plot All Samples ######
set.seed(111)
md = read.table('../EntireMetadata.txt',header=TRUE,sep='\t',comment.char = '')
md = md %>% mutate(LatJit = jitter(Latitude,amount =1),
                   LonJit = jitter(Longitude,amount=1))
#if you want a custom order, add a field with it
ord = read.table('/dss/dsshome1/lxc07/di39dux/merondun/cuculus_plumage/trees/Chr_W_TreeWithMIC-NoChina-SetC_ROOT.order',header=TRUE) %>% arrange(desc(TreeOrder)) %>% dplyr::rename(IDtree = ID)
md = left_join(ord ,md %>% mutate(IDtree = gsub('307_CM_grey','307_CM_OUT',IDtree)))
md = md %>% arrange(desc(TreeOrder))

world <- map_data("world")
sites <- st_as_sf(md, coords = c("LonJit", "LatJit"), 
                  crs = 4326, agr = "constant") 

#add labels
# Calculate the y-coordinate for the labels based on the map extent
y_coordinate <- max(md$Latitude) + 10

# Calculate the minimum and maximum x-coordinates for the labels
min_x <- min(md$Longitude) - 2
max_x <- max(md$Longitude) + 2
label_positions <- seq(min_x, max_x, length.out = nrow(md)) # Calculate the number of evenly spaced x-coordinates for the labels

# Create the new data frame with label positions and connecting line endpoints
label_data <- md %>% ungroup %>%
  mutate(Longitude_label = label_positions,
         Latitude_label = y_coordinate)

#Plot with labels 
hsp = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') + #add basemap
  geom_sf(data = sites , 
          aes(fill=Plumage,shape=Species_Latin),
          size=2,show.legend = T,col='grey20',alpha=0.9) +  #add the geometry for the points, I use pch 21-25 so the points have outlines
  scale_shape_manual(values=sites$Shape,breaks=sites$Species_Latin)+ #custom shape encoded from metadata
  scale_fill_manual(values=sites$JPlumageColor,breaks=sites$Plumage)+ #custom fill encoded from metadata
  scale_color_manual(values=sites$JPlumageColor,breaks=sites$Plumage)+ #custom color encoded from metadata
  geom_segment(data = label_data %>% filter(Plumage == 'rufous' & Country != 'HUN'),  #only show segments for some hepatic individuals, but exclude the hungarian ones because it's too many lines
               aes(x = Longitude_label, y = Latitude_label, xend = LonJit, yend = LatJit, col=Plumage), linetype = 3) + #color by plumage
  geom_rect(data = label_data, inherit.aes = FALSE,
            aes(xmin = min(Longitude_label)-2, xmax = max(Longitude_label)+2, ymin = Latitude_label, ymax = max(md$Latitude) + 25),
            fill = "white", color = NA) +  #add a white rectangle base so that the map doesn't show
  geom_text(data = label_data, aes(x = Longitude_label, y = Latitude_label, label = IDtree, col=Plumage), angle = 90, hjust = -.05,size=1.75) + #add the labels on top
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(md$Longitude)-5, max(md$Longitude)+5), 
           ylim = c(min(md$Latitude)-5, max(md$Latitude)+25), expand = FALSE)+ #expand the map boundaries so that the labels show up, probably an easier way but I just add a gray rectangle in inkscape afterwards in between the labels and the map 
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+ #make the ocean blue
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+ 
  annotation_scale(line_width=0.5,location='tl')+  #add scale bar (requires ggspatial)
  guides(fill=guide_legend(override.aes=list(shape=21))) #ensure the legend has a fill-able shape

hsp

pdf('Plumage_Spatial_HalloweenData_WithMIC-NoChina_SETC_1DEGREEJITTER.pdf',height=4,width=7)
#png('Plumage_Spatial.png',units='in',res=600,height=5,width=9)
hsp
dev.off()

####### Plot Hungary ######
world <- map_data("world")
pt = md %>% filter(Country == 'HUN')
pt = pt %>% mutate(LatJit = jitter(Latitude,amount =0.25),
                   LonJit = jitter(Longitude,amount=0.25)) 
sites <- st_as_sf(pt, coords = c("LonJit", "LatJit"), 
                  crs = 4326, agr = "constant") 
dat = sites
#or single plot
hsp <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') +
  geom_sf(data = dat, 
          aes(fill=Plumage,shape=Species_Latin),
          size=1.5,show.legend = T) +
  scale_shape_manual(values=dat$Shape,breaks=dat$Species_Latin)+
  scale_fill_manual(values=dat$PlumageColor,breaks=dat$Plumage)+
  scale_alpha_manual(values=c(rep(1,22),0.5),guide='none')+
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(dat$Longitude)-1, max(dat$Longitude)+1), 
           ylim = c(min(dat$Latitude)-.5, max(dat$Latitude)+.5), expand = FALSE)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5,location='br')+
  guides(fill=guide_legend(override.aes=list(shape=21)))

hsp

pdf('Plumage_Spatial_HUNGARY_0.25JITTER.pdf',height=2,width=3)
hsp
dev.off()

####### Plot NE China ######
world <- map_data("world")
pt = md %>% filter(Longitude > 110 & Latitude < 55)
pt = pt %>% mutate(LatJit = jitter(Latitude,amount =0.5),
                   LonJit = jitter(Longitude,amount=0.5)) 
sites <- st_as_sf(pt, coords = c("LonJit", "LatJit"), 
                  crs = 4326, agr = "constant") 
dat = sites
#or single plot
hsp <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),col='grey90',fill='white') +
  geom_sf(data = dat, 
          aes(fill=Plumage,shape=Species_Latin),
          size=1.5,show.legend = T,alpha=0.9) +
  scale_shape_manual(values=dat$Shape,breaks=dat$Species_Latin)+
  scale_fill_manual(values=dat$PlumageColor,breaks=dat$Plumage)+
  scale_alpha_manual(values=c(rep(1,22),0.5),guide='none')+
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(dat$Longitude)-1, max(dat$Longitude)+1), 
           ylim = c(min(dat$Latitude)-.5, max(dat$Latitude)+.5), expand = FALSE)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5,location='br')+
  guides(fill=guide_legend(override.aes=list(shape=21)))

hsp

pdf('Plumage_Spatial_NECHINA_0.5JITTER.pdf',height=2,width=3)
hsp
dev.off()

