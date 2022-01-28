###################################################################
# File: amr_genes.R                                               
#                                                                 
# Purpose: Examine the composition and diversity of ARGs in 
#          FL100 cohort
#                                                                 
#                                                                 
# Author: A.Oliver				                                        
# Date: 1/20/21						                                        
#                                                                 
# Inputs (1):                                                     
# (1) FL100_merged_norm_final.csv
# (2) abx_cluster_andrew.csv
#                                                                 
#                     #####################                       
#                                                                 
# Outputs (1):                                                        
# (1) Bar plots of AMR genes
# (2) Diversity plot of AMR genes
# (3) AMR clusters based on percentiles
#                                                                 
# Usage: Run the entire script without changes.                   
##################################################################



#################
## AMR GENES
#################

library(tidyverse)
library(reshape2)
library(ggsci)
library(vegan)
library(bestNormalize)
library(janitor)
library(patchwork)
setwd("/home/datasets/")

## define function for CV
CV <- function(x){
  (sd(x)/mean(x))*100
}

## read in data
amr_genes <- read.csv("amr_genes/FL100_merged_norm_final.csv", check.names = F)
abx_cluster <- read_delim("from_andrew/abx_cluster_andrew.csv", delim = ",")[2:3]

## group AMR by mechanism
amr_genes_mech <- amr_genes %>% 
  select(., -c("MEGID", "Gene", "Group")) %>%
  group_by(., Mechanism) %>% 
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames(., var = "Mechanism") %>%
  t() %>% as.data.frame() %>% clean_names() %>%
  rownames_to_column(., var = "subject_id")
## group AMR by gene
amr_genes_mech <- amr_genes %>% 
  select(., -c("MEGID", "Mechanism", "Group")) %>%
  group_by(., Gene) %>% 
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames(., var = "Gene") %>%
  t() %>% as.data.frame() %>% clean_names() %>%
  rownames_to_column(., var = "subject_id")

########################
## Barplot for mechanism
########################

## group AMR by mechanism
amr_genes_mech <- amr_genes %>% 
  select(., -c("MEGID", "Gene", "Group")) %>%
  group_by(., Mechanism) %>% 
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames(., var = "Mechanism") %>%
  t() %>% as.data.frame() %>% clean_names() %>%
  rownames_to_column(., var = "subject_id")

## clean up after that, rename some columns
amr_genes_mech_melt <- melt(id.vars = "subject_id", data = amr_genes_mech) %>%
  rename(., "Mechanism" = "variable")
amr_genes_mech_melt$Mechanism <- as.character(amr_genes_mech_melt$Mechanism)

## get top amr mechanisms by cohort abundance
top_genes <- amr_genes_mech_melt %>% 
  group_by(., Mechanism) %>%
  summarise(., total = sum(value), ) %>%
  arrange(., desc(total)) %>%
  slice(., (1:9)) %>%
  droplevels()

## View the top mechanisms for statistics
amr_genes_mech_melt %>% 
  group_by(., Mechanism, subject_id) %>%
  summarise(., total = sum(value) ) %>%
  summarise(., average = mean(total), CV = CV(total), sd = sd(total)) %>%
  arrange(desc(average))

## figure out what top mechanisms to plot
amr_genes_mech_melt <- amr_genes_mech_melt %>%
  mutate(., plot = ifelse(Mechanism %in% top_genes$Mechanism, Mechanism, "zlow_abundant"))
amr_genes_mech_melt$plot <- as.factor(amr_genes_mech_melt$plot)
levels(amr_genes_mech_melt$plot)
amr_genes_mech_melt <- amr_genes_mech_melt %>% group_by(., subject_id) %>% 
  mutate(., total_gene_abundance = sum(value))

## get andrew clusters based on quartiles
## ex high is > 0.75 , med = 0.25-0.50, low < 0.25 OF TOTAL AMR ABUNDANCE
amr_genes_mech_melt <- amr_genes_mech_melt %>% 
  ungroup() %>% 
  mutate(., andrew_cluster = ifelse(total_gene_abundance < 49.85, "low", 
                                    ifelse(total_gene_abundance > 60.02, "high", "medium")))

## plot distribution of andrew-clusters (figure 1C)
andrew_cluster_volcano <- amr_genes_mech_melt %>% 
  select(subject_id, total_gene_abundance, andrew_cluster) %>% 
  distinct
andrew_cluster_volcano$andrew_cluster <- factor(andrew_cluster_volcano$andrew_cluster, 
                                                levels = c("low", "medium", "high"),
                                                ordered = TRUE)
figure_1c <- ggplot(data = andrew_cluster_volcano, aes(x = andrew_cluster, y = total_gene_abundance)) +
  geom_violin(aes(fill = as.factor(andrew_cluster))) + 
  geom_point(position = position_jitter(width = .1), alpha = 0.2) + 
  geom_boxplot(width=0.1, aes(fill = "white"), alpha = 0.7, outlier.alpha = 0) + 
  scale_fill_manual(values = c("#374E55", "#DF8F44", "#00A1D5", "white")) + 
  theme_classic() + theme(legend.position = "none") +
  labs(x = "Antibiotic Cluster", y = "Normalized Gene Abundance")

abx_cluster_andrew <- andrew_cluster_volcano %>% select(., -total_gene_abundance) %>% rename(., "cluster" = "andrew_cluster")
## Write main abundance file
write.csv(abx_cluster_andrew, file = "from_andrew/abx_cluster_andrew.csv", sep = ",", quote = F)

## plot
## relative abundance 
amr_relative_abundance <- ggplot(data = amr_genes_mech_melt, aes(x = as.factor(subject_id), weight = value, group = plot)) +
  geom_bar(position = position_fill(), aes(fill = plot), width=1) +
  scale_y_continuous(labels = scales::percent_format()) + 
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank()) + 
  scale_fill_d3()
## by total gene abundance (Figure 1 A)
figure_1a <- ggplot(data = amr_genes_mech_melt, aes(x = reorder(as.factor(subject_id), -total_gene_abundance), y = as.numeric(value), group = plot)) +
  geom_bar(aes(fill = plot), width=1, stat = "identity") +
  coord_cartesian(ylim = c(0, 125), expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank()) + 
  labs(y = "Normalized Gene Abundance", x = "Individual") +
  geom_hline(yintercept = 54.59) +
  geom_hline(yintercept = 49.84) +
  geom_hline(yintercept = 60.03) +
  scale_fill_d3()

## plot the low abundant (Figure 1B)
low_abundant_amr <- amr_genes_mech_melt %>% filter(., plot == "zlow_abundant") %>% select(., -plot)
top_low_abundant_genes <- low_abundant_amr %>% 
  filter(., andrew_cluster == "high") %>%
  group_by(., Mechanism) %>%
  summarise(., total = sum(value)) %>%
  arrange(., desc(total)) %>%
  slice(., (1:9)) %>%
  droplevels() 

low_abundant_amr <- low_abundant_amr %>%
  mutate(., plot = ifelse(Mechanism %in% top_low_abundant_genes$Mechanism, Mechanism, "zlow_abundant"))

low_abundant_amr$plot <- as.factor(low_abundant_amr$plot)

figure_1b <- ggplot(data = subset(low_abundant_amr, andrew_cluster == "high"), aes(x = reorder(as.factor(subject_id), -total_gene_abundance), y = as.numeric(value), group = plot)) +
  geom_bar(aes(fill = plot), width=1, stat = "identity") +
  coord_cartesian(ylim = c(0, 40), expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank()) + 
  labs(y = "Normalized Gene Abundance", x = "Low abundant category in\n High-ARG individuals") +
  geom_vline(xintercept = "6067") +
  geom_vline(xintercept = "6053") +
  #geom_hline(yintercept = 60.03) +
  scale_fill_jco()

#################
## PLOT FIGURE 1
#################

# fig_1_layout <- "
# AAB
# AAC"
# 
# figure_1a + figure_1b + figure_1c + 
#   patchwork::plot_layout(design = fig_1_layout) + 
#   patchwork::plot_annotation(tag_levels = 'A')

### beta diversity
## Permanovas

## group AMR by gene
amr_genes_mech <- amr_genes %>% 
  select(., -c("MEGID", "Mechanism", "Group")) %>%
  group_by(., Gene) %>% 
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames(., var = "Gene") %>%
  t() %>% as.data.frame() %>% clean_names() %>%
  rownames_to_column(., var = "subject_id")

set.seed(999)
amr_genes_mech_beta <- amr_genes_mech %>% column_to_rownames(., var = "subject_id")
amr_genes_mech_beta_permanova <- merge(abx_cluster, amr_genes_mech_beta, by.x = "subject_id", by.y = "row.names")
cluster_permanova <- adonis(amr_genes_mech_beta_permanova[,3:NCOL(amr_genes_mech_beta_permanova)] ~ amr_genes_mech_beta_permanova$cluster, data = amr_genes_mech_beta_permanova, permutations = 999, parallel = 4, method = "euclidean")
cluster_permanova

## check permanova assumptions - FAIL
tmp <- amr_genes_mech_beta_permanova %>% 
  column_to_rownames(., var = "subject_id") %>% 
  select(., 2:390)
tmp_dist <- vegdist(tmp, method = "bray")
tmp <- betadisper(tmp_dist, group = amr_genes_mech_beta_permanova$cluster, type = "centroid")
anova(tmp)
plot(tmp)

## vis of MDS of the data
set.seed(seed = 999)
beta.mds <- metaMDS(amr_genes_mech_beta, distance="bray", k=2)
stressplot(beta.mds)
beta.mds$stress

sites <- as.data.frame(scores(beta.mds, display = "sites"))
species <- as.data.frame(scores(beta.mds, display = "species"))

nmds.sites <- merge(sites, abx_cluster, by.x = "row.names", by.y = "subject_id")
nmds.sites$cluster <- factor(nmds.sites$cluster, levels = c("low", "medium", "high"), ordered = T)

amr_nmds <- ggplot() + 
  geom_point(data = nmds.sites, aes(NMDS1, NMDS2, color = cluster), size = 2, alpha = 0.8) + 
  theme_bw() + theme(panel.border = element_rect(colour = "black", size = 1.5), 
                     panel.grid.minor = element_blank(),
                     legend.position = "none",
                     panel.background = element_blank()) + scale_color_jama()

beta.mds$stress

################################
### alpha diversity - Supp Fig 1
################################

## read in data
amr_genes <- read.csv("amr_genes/FL100_merged_norm_final.csv", check.names = F)
abx_cluster <- read_delim("from_andrew/abx_cluster_andrew.csv", delim = ",")[2:3]

## group AMR by group
amr_genes_mech <- amr_genes %>% 
  select(., -c("MEGID", "Mechanism", "Group")) %>%
  group_by(., Gene) %>% 
  summarise(across(where(is.numeric), sum)) %>%
  column_to_rownames(., var = "Gene") %>%
  t() %>% as.data.frame() %>% clean_names() %>%
  rownames_to_column(., var = "subject_id")

## melt
amr_genes_mech_melt <- melt(id.vars = "subject_id", data = amr_genes_mech) %>%
  rename(., "Gene" = "variable")
amr_genes_mech_melt$Gene <- as.character(amr_genes_mech_melt$Gene)
amr_genes_mech_melt <- amr_genes_mech_melt %>% group_by(., subject_id) %>% 
  mutate(., total_gene_abundance = sum(value))

## get andrew clusters based on quartiles
## ex high is > 0.75 , med = 0.25-0.50, low < 0.25 OF TOTAL AMR ABUNDANCE
amr_genes_mech_melt <- amr_genes_mech_melt %>% 
  ungroup() %>% 
  mutate(., andrew_cluster = ifelse(total_gene_abundance < 49.85, "low", 
                                    ifelse(total_gene_abundance > 60.02, "high", "medium")))


amr_genes_mech_beta <- amr_genes_mech %>% column_to_rownames(., var = "subject_id")

## calculate alpha diveristy
shannon_div <- diversity(amr_genes_mech_beta, "shannon") 
richness <- specnumber(amr_genes_mech_beta)
Evenness <- shannon_div/log(richness)
all_diversity <- cbind(Evenness, richness, shannon_div)

## merge with Nutritional Metadata
## read in metadata clusters
abx_clusters <- read_csv("../datasets/from_andrew/abx_cluster_andrew.csv")[2:3]
alpha_diversity_data <- merge(abx_clusters, all_diversity, by.x = "subject_id", by.y = "row.names")

## This is for ggplot. Change the measure vars for what y-axis you want to plot
alpha_diversity_data_melt <- melt(data = alpha_diversity_data, id.vars = c("subject_id", "cluster"), measure.vars = c("Evenness", "richness", "shannon_div"))
alpha_diversity_data_melt$cluster <- factor(alpha_diversity_data_melt$cluster, levels = c("low", "medium", "high"), ordered = T)

## make single alpha diversity box plots
supplemental_fig1 <- ggplot(data = alpha_diversity_data_melt) +
  aes(x = cluster, y = value, fill = cluster) +
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) +
  ##geom_jitter(width = 0.15, alpha = 0.2) +
  labs(x = '') +
  theme_bw(base_size = 14) + theme(legend.position = "none") + facet_wrap(.~variable, scales = "free") + scale_fill_jama()

## test alpha diversity
alpha_diversity_data_melt_test <- subset(alpha_diversity_data_melt, alpha_diversity_data_melt$variable == "shannon_div")
leveneTest(alpha_diversity_data_melt_test$value, group = alpha_diversity_data_melt_test$cluster)
mod <- aov(value ~ cluster, data = alpha_diversity_data_melt_test)
resids <- residuals(mod)
shapiro.test(resids)
## transform
bestNormalize(alpha_diversity_data_melt_test$value)
alpha_diversity_data_melt_test$value_boxcox <- (bestNormalize::boxcox(alpha_diversity_data_melt_test$value))$x.t
mod <- aov(value_boxcox ~ cluster, data = alpha_diversity_data_melt_test)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(alpha_diversity_data_melt_test$value_boxcox, group = alpha_diversity_data_melt_test$cluster)
summary(mod)
TukeyHSD(x = mod)
## describe
alpha_diversity_data_melt_test %>% group_by(., cluster) %>% summarise(., average = mean(value),
                                  
                                                                                                           sd = sd(value))
##################################
## WRITE AMR ABUNDANCE FILE FOR ML 
##################################
## this analysis is playground stuff

amr_genes_mech_ml <- amr_genes_mech %>% column_to_rownames(., var = "subject_id")
amr_genes_mech_ml <- merge(abx_cluster, amr_genes_mech_beta, by.x = "subject_id", by.y = "row.names")
amr_genes_mech_ml <- amr_genes_mech_ml %>%
  mutate(., cluster = ifelse(cluster == "low", 0, 
                             ifelse(cluster == "medium", 1, 2)))

abx_cluster_genes_post_process <- preprocess_data(dataset = amr_genes_mech_ml,
                                                     method = NULL,
                                                     outcome_colname = "cluster", 
                                                     collapse_corr_feats = T, 
                                                     remove_var = "zv")

pre_corr_clean_genes <- abx_cluster_genes_post_process$dat_transformed
## get the correlated features to a dataframe
cor_tmp_genes <- mikropml:::group_correlated_features(abx_cluster_genes_post_process$dat_transformed, 
                                                corr_thresh = 0.80)
cor_tmp_genes <- as.data.frame(cor_tmp_genes)
cor_tmp_genes <- cor_tmp_genes %>% separate(., col = cor_tmp_genes, into = c("keep", "co-correlated"), sep = "\\|", extra = "merge")
## filter out only 1 of the co-correlated groups, results in 146 features (including subject ID)
corr_clean_genes <- pre_corr_clean_genes %>% select(., cor_tmp_genes$keep)
## write to file
corr_clean_genes_write <- corr_clean_genes %>% select(., -subject_id)
## All clusters - 187 samples 
#write.csv(corr_clean_genes_write, file = "output_for_ML/AMR_genes_input.csv", quote = F, row.names = FALSE)
#write.csv(cor_tmp_genes, file = "output_for_ML/co-correlated-features_amr_genes.csv", quote = F, row.names = FALSE)


  