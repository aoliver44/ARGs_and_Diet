###################################################################
# File: Supp_figure_dheas.R                             
#                                                                 
# Purpose: Examine abundance of DHEAS across ARG clusters                                 
#                                                                 
#                                                                 
# Author: A.Oliver				                                        
# Date: 1/20/21						                                        
#                                                                 
# Inputs (1):                                                     
# (1) Be able to source merging_features_new_data script          
#      successfully 
# (2) clinical_range_dheas.csv - clinical range of dheas
#                                                                 
#                     #####################                       
#                                                                 
# Outputs (1):                                                        
# (1) supplemental figure 3 about DHEAS abundance
#                                                                 
# Usage: Run the entire script without changes.                   
###################################################################

#########################
## Supp. Figure 3 - DHEAS
#########################

library(tidyverse)
library(ggpubr)
setwd("/home/datasets/new_datasets/")
#source(file = "/home/scripts/merging_features_new_data.R")
## Figure 3A - DHEAS Raw
## make sure you generate files from merging_features_new_data.R
boxplots <- abx_cluster_features %>% select(., subject_id, cluster, dheas_bd1) %>% drop_na() 
boxplots <- boxplots %>% mutate(., cluster = ifelse(cluster == "soft", 0,
                                                                    ifelse(cluster == "normal", 1, 
                                                                           ifelse(cluster == "hard", 2, cluster))))

boxplots$cluster <- factor(x = boxplots$cluster, levels = c("low", "medium", "high"), ordered = T)
dheas_boxplot <- boxplots %>%
  ggplot() + aes(x = as.factor(cluster), y = as.numeric(dheas_bd1)) + 
  geom_boxplot(aes(fill = as.factor(cluster)),outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.2) +
  scale_fill_jama() + theme_bw(base_size = 20) +
  labs(x = "AMR Cluster", y = "DHEA-S (mcg/dl)") + theme(legend.position = "none") +
  stat_compare_means(method = "kruskal") 

## Figure 3A - DHEAS clinical range
clin_range <- read.table(file = "clinical_range_dheas.csv", sep = ",", header = T)
clin_range$cluster <- factor(clin_range$cluster, levels = c("low", "medium", "high", 
                                                            "men_low_bound", "men_high_bound", 
                                                            "women_low_bound", 
                                                            "women_high_bound"), ordered = T)
clin_range_plot <- clin_range %>% ggplot() + aes(x = age, y = dheas_bd1, linetype = factor(cluster)) +
  geom_point(aes(shape = as.factor(sex), color = cluster), size = 2.5) +
  geom_line() +
  scale_shape_manual(values = c(0,17)) + scale_color_jama() +
  scale_linetype_manual(values = c("blank", "blank", "blank", "solid", "dashed", "solid", "dashed")) + theme_bw()
