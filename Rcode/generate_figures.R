###################################################################
# File: figures.R                                               
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
library(kableExtra)

##########
## Table 1
##########

source(file = "/home/scripts/merging_features_new_data.R")
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

source(file = "/home/scripts/amr_genes.R")

fig_1_layout <- "
AAB
AAC"

figure_1a + figure_1b + figure_1c + 
  patchwork::plot_layout(design = fig_1_layout) + 
  patchwork::plot_annotation(tag_levels = 'A')

############
## Figure 2
############

source(file = "/home/scripts/alpha_diversity.R")

fig_2_layout <- "
ABB"

figure_2a + figure_2b +
  patchwork::plot_layout(design = fig_2_layout) + 
  patchwork::plot_annotation(tag_levels = 'A')

############
## Table 2
############

source(file = "/home/scripts/directed_hypothesis_testing.R")

directed_hypothesis_table %>% kbl() %>% kable_classic()


############
## Figure 3
############

source(file = "/home/scripts/blythe_plots.R")

fig_3_layout <- "
A
B"

figure_3ab + figure_3cd +
  patchwork::plot_layout(design = fig_3_layout) + 
  patchwork::plot_annotation(tag_levels = 'A')

############
## Figure 4
############

## note, figure 4b and 4d were generated using python

#source(file = "/home/scripts/alpha_diversity.R")
source(file = "/home/scripts/plot_shap_values.R")

fig_4_layout <- "
AABC
AABC"

figure_4a + figure_4c_1 + figure_4c_2 +
  patchwork::plot_layout(design = fig_4_layout) + 
  patchwork::plot_annotation(tag_levels = 'A')
  

  
########################
## Supplemental Figure 1
########################

#source(file = "/home/scripts/amr_genes.R")
supplemental_fig1

########################
## Supplemental Figure 2
########################

source(file = "/home/scripts/Supp_figure_lefse.R")
lefse_plot

########################
## Supplemental Figure 3
########################

source(file = "/home/scripts/Supp_figure_dheas.R")

supp_3_layout <- "
AB"

dheas_boxplot + clin_range_plot +
  patchwork::plot_layout(design = supp_3_layout) + 
  patchwork::plot_annotation(tag_levels = 'A')


