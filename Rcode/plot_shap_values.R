###################
## Plot SHAP Values
###################


library(tidyverse)
library(patchwork)
library(reshape2)
library(viridis)
setwd("/home/machine_learning/new_dataset_models/shap_values/")

## read in shap data - diet + lifestyle - bi (low vs medium)
diet_lifestyle_bi_shap <- read_delim(file = "diet-life_bi_class1_shapvals-NEW.csv")
diet_lifestyle_bi_shap <- diet_lifestyle_bi_shap %>% select(., -c(`...1`))
## read in raw data - diet + lifestyle - bi (low vs medium)
diet_lifestyle_bi_raw <- read_delim(file = "/home/datasets/new_datasets/output_for_ml/diet-life_input_bi.csv")
#diet_lifestyle_bi_raw <- diet_lifestyle_bi_raw %>% select(., -c(`...1`, `Unnamed: 0`))

#########################
## MODIFIED BEESWARM PLOT
#########################

## features to pull and display
## for diet lifestyle - both ASA and FFQ - BI
diet_lifestyle_bi_features <- c("pd_whole_tree_total", "dheas_bd1", 
                                "ur_epi_ug_gcreat", "dt_fiber_sol_per_kcal","per_kcal_fiber_tnfs", "dt_p225", 
                                "f15d1", "pd_whole_tree_fat", "pf_seafd_low", "plasma_lbp_bd1", "pf_organ", "hdl_bd1", "coumestrol", "avg_fibe_tnfs")


feature_plots_shap = list()
feature_plots_raw = list()
diet_lifestyle_bi_plot = list()
for (feature in diet_lifestyle_bi_features) {

## pull features and combine
diet_lifestyle_bi_shap_loop <- diet_lifestyle_bi_shap %>% select(., feature)
diet_lifestyle_bi_shap_loop <- melt(diet_lifestyle_bi_shap_loop)
diet_lifestyle_bi_shap_loop <- diet_lifestyle_bi_shap_loop %>%
  rename(., "shap_category" = "variable", "shap_values" = "value")
diet_lifestyle_bi_raw_loop <- diet_lifestyle_bi_raw %>% select(., feature)
diet_lifestyle_bi_raw_loop <- melt(diet_lifestyle_bi_raw_loop)
diet_lifestyle_bi_raw_loop <- diet_lifestyle_bi_raw_loop %>%
  rename(., "raw_category" = "variable", "raw_values" = "value")
diet_lifestyle_bi_plot[[feature]] <- cbind(diet_lifestyle_bi_shap_loop, diet_lifestyle_bi_raw_loop)
diet_lifestyle_bi_plot[[feature]]$cluster <- as.factor(diet_lifestyle_bi_raw$cluster)

print(feature)

midpoint <- mean(diet_lifestyle_bi_plot[[feature]]$raw_values)
low <- min(diet_lifestyle_bi_plot[[feature]]$raw_values)
high <- max(diet_lifestyle_bi_plot[[feature]]$raw_values)
print(midpoint)
print(low)
print(high)
feature_plots_shap[[feature]] <- print(ggplot(diet_lifestyle_bi_plot[[feature]], aes(x = shap_category, y = shap_values)) +
                                        geom_point(position = position_jitterdodge(jitter.width = 1), aes(color = log(raw_values)), size = 0.8) + 
                                         theme_bw() +
                                         coord_flip(ylim = c(-0.09, 0.09)) +
                                         labs(x= "", y = "", title = feature) +
                                         theme_minimal() +
                                         theme(axis.text.y = element_blank(), legend.position = "none", 
                                               plot.margin=grid::unit(c(0,0,0,4.8), "cm"), axis.text.x = element_blank(),
                                               axis.ticks = element_blank(), panel.grid.minor = element_blank(),
                                               panel.background = element_blank(), panel.grid.major.y = element_line(colour = "black"),
                                               panel.grid.major = element_line(colour = "black"), plot.title = element_text(vjust = -76.5, hjust = -0.3)) + 
                                         scale_color_viridis(option = "C"))
                                         #scale_color_gradient2(low="#0092F7", mid = "#AC40A1", high="#F90050",midpoint = midpoint, limits=c(low, high)))

diet_lifestyle_bi_plot[[feature]] <- diet_lifestyle_bi_plot[[feature]] %>% mutate(., cluster = ifelse(cluster == 1, "medium", "low"))
diet_lifestyle_bi_plot[[feature]]$cluster <- factor(x = diet_lifestyle_bi_plot[[feature]]$cluster, levels = c("medium", "low"), ordered = T)
feature_plots_raw[[feature]] <- print(ggplot(diet_lifestyle_bi_plot[[feature]], aes(x = as.factor(cluster), y = raw_values)) +
                                        #geom_violin(aes(fill = as.factor(cluster)), alpha = 0.6) + 
                                        #geom_point(position = position_jitter(width = .1), alpha = 0.2) + 
                                        geom_boxplot(width=0.8, aes(fill = as.factor(cluster)), alpha = 1, outlier.alpha = 0) + 
                                        scale_fill_manual(values = c("#DF8F44", "#374E55","white")) +
                                        theme_bw() +
                                        theme(axis.text.y = element_blank(), 
                                              axis.ticks.y = element_blank(),
                                              legend.position = "none", plot.margin=grid::unit(c(0,0,0,0), "mm"),
                                              panel.grid = element_blank(), panel.background = element_blank()) +
                                        labs(x = "", y = "") + scale_y_continuous(position = "right") + coord_flip())



}

fig_layout_final <- "
AAAAB
CCCCD
EEEEF
GGGGH
IIIIJ
KKKKL
MMMMN
OOOOP
QQQQR
SSSST
UUUUV
WWWWX
YYYYZ"



figure_4a <- feature_plots_shap$pd_whole_tree_total + feature_plots_raw$pd_whole_tree_total +
  feature_plots_shap$dheas_bd1 + feature_plots_raw$dheas_bd1 +
  feature_plots_shap$ur_epi_ug_gcreat + feature_plots_raw$ur_epi_ug_gcreat +
  feature_plots_shap$dt_fiber_sol_per_kcal + feature_plots_raw$dt_fiber_sol_per_kcal +
  feature_plots_shap$per_kcal_fiber_tnfs + feature_plots_raw$per_kcal_fiber_tnfs +
  feature_plots_shap$dt_p225 + feature_plots_raw$dt_p225 +
  feature_plots_shap$f15d1 + feature_plots_raw$f15d1 +
  feature_plots_shap$pd_whole_tree_fat + feature_plots_raw$pd_whole_tree_fat +
  feature_plots_shap$pf_seafd_low + feature_plots_raw$pf_seafd_low +
  feature_plots_shap$plasma_lbp_bd1 + feature_plots_raw$plasma_lbp_bd1 +
  feature_plots_shap$pf_organ + feature_plots_raw$pf_organ +
  feature_plots_shap$hdl_bd1 + feature_plots_raw$hdl_bd1 +
  feature_plots_shap$coumestrol + feature_plots_raw$coumestrol +
  #feature_plots_shap$avg_fibe_tnfs + feature_plots_raw$avg_fibe_tnfs +
  plot_layout(design = fig_layout_final)



