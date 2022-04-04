###################################################################
# File: generate_figures.R                                               
#                                                                 
# Purpose: Generate most of the figures in the manuscript                                
#                                                                 
#                                                                 
# Author: A.Oliver				                                        
# Date: 1/27/21						                                        
#                                                                 
# Inputs (1):                                                     
# (1) merging_features_new_data.R
# (2) amr_genes.R
# (3) alpha_diversity.R
# (4) directed_hypothesis_testing.R
# (5) blythe_plots.R
# (6) plot_shap_values.R
#
#                     #####################                       
#                                                                 
# Outputs (1):                                                        
# (1) ~8 figures
#                                                                 
# Usage: Run the entire script without changes.                   
##################################################################

library(gtsummary)
library(tidyverse)
library(patchwork)
library(reshape2)
library(ggsci)
##########
## Table 1
##########

source(file = "/home/Rcode/merging_features_new_data.R")
abx_cluster$cluster <- factor(abx_cluster$cluster, 
                              levels = c("low", "medium", "high"), ordered = T)
merge(abx_cluster, bmi_etc, by = "subject_id", all.x = T) %>%
  merge(., ethnicities, by = "subject_id", all = T) %>%
  merge(., country_origin, by = "subject_id", all = T) %>%
  filter(., cluster != "") %>% 
  select(., cluster, age, bmi_final, sex, ethnicity, birth_country) %>%
  mutate(., sex = ifelse(sex == 1, "male", "female")) %>%
  tbl_summary(., by = "cluster", 
              statistic = list(all_continuous() ~ "{mean} ({min}, {max})", 
                               all_categorical() ~ "{n} ({p}%)")) 

############
## Figure 1
############

source(file = "/home/Rcode/amr_genes.R")

fig_1_layout <- "
AAB
AAC"

figure_1a + figure_1b + figure_1c + 
  patchwork::plot_layout(design = fig_1_layout) + 
  patchwork::plot_annotation(tag_levels = 'A')

############
## Figure 2
############

source(file = "/home/Rcode/alpha_diversity.R")

fig_2_layout <- "
ABB"

figure_2a + figure_2b +
  patchwork::plot_layout(design = fig_2_layout) + 
  patchwork::plot_annotation(tag_levels = 'A')

############
## Table 2
############

source(file = "/home/Rcode/directed_hypothesis_testing.R")


############
## Figure 3
############

source(file = "/home/Rcode/blythe_plots.R")

fig_3_layout <- "
ABB"

figure_3ab + figure_3cd +
  patchwork::plot_layout(design = fig_3_layout) + 
  patchwork::plot_annotation(tag_levels = 'A')

############
## Figure 4
############

## NOTE, figure 4b and 4d were generated using python
## shap bar outputs

#source(file = "/home/Rcode/alpha_diversity.R")
source(file = "/home/Rcode/plot_shap_values.R")

fig_4_layout <- "
BC"

## Note this is a small version of the larger figure
## which was edited in Affinity designer
figure_4a

figure_4c_1 + figure_4c_2 +
  patchwork::plot_layout(design = fig_4_layout) + 
  patchwork::plot_annotation(tag_levels = 'A')
  

########################
## Supplemental Figure 1
########################

## Table 1
#source(file = "/home/Rcode/merging_features_new_data.R")
abx_cluster$cluster <- factor(abx_cluster$cluster, 
                              levels = c("low", "medium", "high"), ordered = T)
merge(abx_cluster, bmi_etc, by = "subject_id", all.x = T) %>%
  merge(., ethnicities, by = "subject_id", all = T) %>%
  merge(., country_origin, by = "subject_id", all = T) %>%
  filter(., cluster != "") %>% 
  select(., cluster, age, bmi_final, sex, ethnicity, birth_country) %>%
  mutate(., sex = ifelse(sex == 1, "male", "female")) %>%
  tbl_summary(., by = "cluster", 
              statistic = list(all_continuous() ~ "{mean} ({min}, {max})", 
                               all_categorical() ~ "{n} ({p}%)")) 

## Table 2

abx_cluster_features_no_na %>%
  filter(., cluster != "") %>%
  mutate(., cluster = ifelse(cluster == 0, "low", 
                             ifelse(cluster == 1, "medium", "high"))) %>%
  select(., cluster, age, bmi_final, sex, ethnicity, birth_country) %>%
  mutate(., sex = ifelse(sex == 1, "male", "female")) %>%
  tbl_summary(., by = "cluster", 
              statistic = list(all_continuous() ~ "{mean} ({min}, {max})", 
                               all_categorical() ~ "{n} ({p}%)")) 

## Figure
total_boxplots <- merge(abx_cluster, bmi_etc, by = "subject_id", all.x = T) %>%
  merge(., ethnicities, by = "subject_id", all = T) %>%
  merge(., country_origin, by = "subject_id", all = T) %>%
  filter(., cluster != "") %>% 
  select(., cluster, age, bmi_final) %>% 
  mutate(., study = "Total_Cohort")

ml_boxplots <- abx_cluster_features_no_na %>%
  filter(., cluster != "") %>%
  mutate(., cluster = ifelse(cluster == 0, "low", 
                             ifelse(cluster == 1, "medium", "high"))) %>%
  select(., cluster, age, bmi_final) %>%
  mutate(., study = "ML_Cohort")

sup1_boxplots <- rbind(total_boxplots, ml_boxplots)
sup1_boxplots <- melt(sup1_boxplots, id.vars = c("cluster", "study"))
sup1_boxplots$cluster <- factor(sup1_boxplots$cluster, 
                              levels = c("low", "medium", "high"), ordered = T)
sup1_boxplots %>% ggplot() + aes(x = study, y = value) + 
  geom_boxplot(aes(fill = cluster)) +
  facet_grid(. ~ variable) + 
  ggsci::scale_fill_jama()

########################
## Supplemental Figure 2
########################

#source(file = "/home/Rcode/amr_genes.R")
source(file = "/home/Rcode/suppp_figure_cat-bat.R")

supp_1_layout <- "
AA
BB
CC"
sup_fig_1a_top + sup_fig_1a_bottom + supplemental_fig1 +
  patchwork::plot_layout(design = supp_1_layout) + 
  patchwork::plot_annotation(tag_levels = 'A')
  

########################
## Supplemental Figure 3
########################

source(file = "/home/Rcode/Supp_figure_lefse.R")
lefse_plot


########################
## Supplemental Figure 4
########################
#source(file = "/home/Rcode/blythe_plots.R")
supplemental_fig4a

########################
## Supplemental Figure 5
########################

source(file = "/home/Rcode/Supp_figure_dheas.R")

supp_3_layout <- "
AB"

dheas_boxplot + clin_range_plot +
  patchwork::plot_layout(design = supp_3_layout) + 
  patchwork::plot_annotation(tag_levels = 'A')

########################
## Supplemental Table 1 
########################

source(file = "/home/Rcode/data_dictionary.R")

View(paper_dictionary_clean)

########################
## Supplemental Table 2 
########################

## there is some slight variability in the features that get
## kept or co-correlated
View(cor_tmp)

