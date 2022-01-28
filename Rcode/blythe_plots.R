###################################################################
# File: blythe_plots.R                                               
#                                                                 
# Purpose: Plot POLR regression results (smoothed scatter plots
#          of significant results from POLR analysis)                                
#                                                                 
#                                                                 
# Author: A.Oliver				                                        
# Date: 1/27/21						                                        
#                                                                 
# Inputs (1):                                                     
# (1) FL100_merged_norm_final.csv
# (2) abx_cluster_andrew.csv
#
#                     #####################                       
#                                                                 
# Outputs (1):                                                        
# (1) Figure with 4 panels, of 4 significant results from POLR
#      analysis.
#                                                                 
# Usage: Run the entire script without changes.                   
##################################################################

#################
## BLYTHE PLOTS
#################

library(tidyverse)
library(reshape2)
library(ggsci)
library(ggplot2)
library(reshape2)
setwd("/home/datasets/")

## read in data
amr_genes <- read.csv("amr_genes/FL100_merged_norm_final.csv", check.names = F)
abx_cluster <- read_delim("from_andrew/abx_cluster_andrew.csv", delim = ",")[2:3]

## group AMR by MEGID
amr_genes_MEGID <- amr_genes %>% 
  select(., -c("Gene", "Mechanism", "Group")) %>%
  group_by(., MEGID) %>% 
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames(., var = "MEGID") %>%
  t() %>% as.data.frame() %>% clean_names() %>%
  rownames_to_column(., var = "subject_id")

amr_genes_mech <- amr_genes %>% 
  select(., -c("Gene", "MEGID", "Group")) %>%
  group_by(., Mechanism) %>% 
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames(., var = "Mechanism") %>%
  t() %>% as.data.frame() %>% clean_names() %>%
  rownames_to_column(., var = "subject_id")

## merge with relevant dietary data
## (run merging features new data script first to genenerate for_directed_hypot
## -hesis_testing dataframe)
blythe_data <- merge(abx_cluster, ffq_data, by = "subject_id")
blythe_data <- blythe_data %>%
  select(., dt_kcal, dt_fibe, dt_fiber_sol, subject_id, cluster)

meg_1039_data <- amr_genes_MEGID %>%
  select(., subject_id, meg_1039)
meg_1039_data <- merge(meg_1039_data, blythe_data, by = "subject_id")

figure_3ab <- meg_1039_data %>% select(., dt_kcal, dt_fiber_sol, meg_1039, cluster, subject_id) %>%
  melt(., id.vars = c("subject_id", "cluster", "meg_1039")) %>%
  ggplot() + aes(x = value, y = log(meg_1039 + 0.1)) +
  geom_point(aes(color = cluster)) +
  geom_smooth(se = F, color = "red") + 
  facet_wrap(. ~ variable, scales = "free") +
  theme_bw() +
  scale_color_jama() + theme(panel.grid.minor = element_blank())

meg_1488_data <- amr_genes_mech %>%
  select(., subject_id, multi_metal_resistance)
meg_1488_data <- merge(meg_1488_data, blythe_data, by = "subject_id")

figure_3cd <- meg_1488_data %>% select(., dt_fibe, dt_fiber_sol, multi_metal_resistance, cluster, subject_id) %>%
  melt(., id.vars = c("subject_id", "cluster", "multi_metal_resistance")) %>%
  ggplot() + aes(x = value, y = log(multi_metal_resistance + 0.1)) +
  geom_point(aes(color = cluster)) +
  geom_smooth(se = F, color = "red") + 
  facet_wrap(. ~ variable, scales = "free") +
  theme_bw() +
  scale_color_jama() + theme(panel.grid.minor = element_blank())

