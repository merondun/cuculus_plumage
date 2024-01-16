#### dA Ratios
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/plumage/toes/martin_dxy_singletons/out')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(gghalves)
library(boot)
library(ggdist)

#read in data 
m3 = read_tsv('Sensitivity_DXY.input') #or start here 
m3 = m3 %>% mutate(chrnum = ifelse(chr == 'Z',40,ifelse(chr == 'W',41,ifelse(chr== 'MT',42,chr)))) %>% 
  mutate_at('chrnum',as.numeric) %>% filter(chrnum < 21 | chrnum > 39)

m3 = m3 %>% mutate(RatioGroup = ifelse(Group == 'CG_CH' | Group == 'OG_OH','Plumage',
                                       ifelse(Group == 'OG_CG' | Group == 'OH_CH','Species','DROP')))

rat = m3 %>% filter(RatioGroup != 'DROP')
rat = rat %>% mutate(Da = DXY - ((piA+piB)/2),
                     Da = ifelse(Da <= 0 , 1e-5, round(Da,5)),
                     site = paste0(chr,'__',start))

# calculate ratio plumage/species Da
set.seed(123)

# Compute ratio for a given data sample
compute_ratio <- function(data) {
  da_plumage <- mean(data$Da[data$RatioGroup == 'Plumage'], na.rm = TRUE)
  da_species <- mean(data$Da[data$RatioGroup == 'Species'], na.rm = TRUE)
  return(da_plumage / da_species)
}

# Bootstrapping function
boot_function <- function(data, indices) {
  resampled_data <- data[indices,]
  return(compute_ratio(resampled_data))
}

# Function to apply bootstrapping per group
bootstrap_per_group <- function(group_data) {
  boot_results <- boot(group_data, boot_function, R = 1000)
  mean_ratio <- mean(boot_results$t)
  lower_bound <- quantile(boot_results$t, 0.025)
  upper_bound <- quantile(boot_results$t, 0.975)
  return(data.frame(mean = mean_ratio, lower = lower_bound, upper = upper_bound))
}

# Grouping by AvZ and applying the bootstrap function
AvZ = data %>%
  group_by(AvZ) %>%
  do(bootstrap_per_group(.))

#save it, also can start here 
write_tsv(AvZ,file='Bootstrapped_Compartment_AvZ.txt')
AvZ = read_tsv('Bootstrapped_Compartment_AvZ.txt')

bds = AvZ %>%
  filter(AvZ != 'Z') %>% 
  ggplot(aes(x=AvZ,y=mean,ymin=lower,ymax=upper,col=AvZ,
             label=paste0(round(mean,2),' (',round(lower,2),' - ',round(upper,2),')')))+
  geom_errorbar(width=0.25)+
  geom_point(size=3)+
  geom_label(aes(y=mean+1.5),size=2)+
  geom_hline(yintercept=1,lty=2)+
  xlab('')+ylab("Da Ratio: Plumage / Species \n Bootstrapped mean (95% CI)")+
  scale_color_manual(values=viridis(2,option='turbo'))+
  theme_bw()+
  theme(legend.position='none')
bds

pdf('~/merondun/cuculus_plumage/divergence_dating/DaRatio_PlumageSpecies-SINGLETONS_2023AUG19.pdf',height=4,width=3)
bds
dev.off()

AvZ
# A tibble: 3 × 4
# AvZ        mean  lower  upper
# <chr>     <dbl>  <dbl>  <dbl>
#   1 Autosome 0.0475 0.0468 0.0482
# 2 W        8.59   8.02   9.20  
# 3 Z        0.0249 0.0233 0.0264
