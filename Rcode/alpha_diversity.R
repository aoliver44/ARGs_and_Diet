######################
## Microbiome analyses
######################

library(vegan)
library(ggplot2)
#library(tidyverse)
library(nlme)
library(reshape2)
library(janitor)

## calculate alpha diveristy
mpa_rare <- readRDS("/home/Kracken/mpa_rare_5833371_perm_5.rds")
shannon_div <- diversity(mpa_rare, "shannon") 
richness <- specnumber(mpa_rare)
Evenness <- shannon_div/log(richness)
all_diversity <- cbind(Evenness, richness, shannon_div)

## merge with Nutritional Metadata
## read in metadata clusters
abx_clusters <- read_csv("/home/datasets/new_datasets/abx_cluster_andrew.csv")
abx_clusters <- abx_clusters %>% select(., subject_id, cluster)
alpha_diversity_data <- merge(abx_clusters, all_diversity, by.x = "subject_id", by.y = "row.names")

## This is for ggplot. Change the measure vars for what y-axis you want to plot
alpha_diversity_data_melt <- melt(data = alpha_diversity_data, id.vars = c("subject_id", "cluster"), measure.vars = c("Evenness", "richness", "shannon_div"))
alpha_diversity_data_melt$cluster <- factor(alpha_diversity_data_melt$cluster, levels = c("low", "medium", "high"), ordered = T)

## make single alpha diversity box plots
figure_2a <- alpha_diversity_data_melt %>% filter(., variable == "shannon_div") %>%
  ggplot() +
  aes(x = cluster, y = value, fill = cluster) +
  geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(), alpha = 0.3) +
  ##geom_jitter(width = 0.15, alpha = 0.2) +
  labs(x = '') +
  theme_bw(base_size = 14) + 
  theme(legend.position = "none") + 
  #facet_wrap(.~variable, scales = "free") + 
  scale_fill_jama()

## test alpha diversity
## normality of residuals
alpha_mod <- aov(alpha_diversity_data$shannon_div ~ alpha_diversity_data$cluster)
shannon_residuals <- residuals(alpha_mod)
shapiro.test(shannon_residuals) # assumption not met
car::leveneTest(alpha_diversity_data$shannon_div ~ alpha_diversity_data$cluster) # assumption met
kruskal.test(alpha_diversity_data$shannon_div ~ alpha_diversity_data$cluster)
dunn.test::dunn.test(alpha_diversity_data$shannon_div, g = alpha_diversity_data$cluster, method = "bh")

#### beta diversity
## Permanovas
set.seed(999)
beta_diversity_data <- merge(abx_clusters, mpa_rare, by.x = "subject_id", by.y = "row.names")
cluster_permanova <- adonis(beta_diversity_data[,3:NCOL(beta_diversity_data)] ~ beta_diversity_data$cluster, data = beta_diversity_data, permutations = 999, parallel = 4, method = "bray")
cluster_permanova

## check permanova assumptions
tmp <- beta_diversity_data %>% column_to_rownames(., var = "subject_id") %>% select(., 2:8943)
tmp_dist <- vegdist(tmp, method = "bray")
tmp <- betadisper(tmp_dist, group = beta_diversity_data$cluster, type = "median")
anova(tmp)

## vis of MDS of the data
set.seed(seed = 999)
beta.mds <- metaMDS(mpa_rare, distance="bray", k=3)
stressplot(beta.mds)

sites <- as.data.frame(scores(beta.mds, display = "sites"))
species <- as.data.frame(scores(beta.mds, display = "species"))

nmds.sites <- merge(sites, abx_clusters, by.x = "row.names", by.y = "subject_id")
nmds.sites$cluster <- factor(nmds.sites$cluster, levels = c("low", "medium", "high"), ordered = T)

figure_2b <- ggplot() + 
  geom_point(data = nmds.sites, aes(NMDS1, NMDS2, color = cluster), size = 2, alpha = 0.8) + 
  theme_bw() + theme(panel.border = element_rect(colour = "black", size = 1.5), 
                     panel.grid.minor = element_blank(),
                     legend.position = "none",
                     panel.background = element_blank()) + scale_color_jama()

beta.mds$stress

#### Vector plot from ML results
## Permanovas
set.seed(999)
vector_data <- read_delim(file = "/home/datasets/output-for-ML/microbiome_family.csv")
vector_data_species <- vector_data %>% select(., -cluster)
## vis of MDS of the data
set.seed(seed = 999)
beta.mds <- metaMDS(vector_data, distance="bray", k=3)
stressplot(beta.mds)
beta.mds$stress

sites <- as.data.frame(scores(beta.mds, display = "sites"))
species <- as.data.frame(scores(beta.mds, display = "species"))

nmds.sites <- cbind(vector_data$cluster, sites)
nmds.sites$cluster <- factor(nmds.sites$`vector_data$cluster`, levels = c("0", "1", "2"), ordered = T)

ggplot() + 
  geom_point(data = nmds.sites, aes(NMDS1, NMDS2, color = cluster), size = 2, alpha = 0.8) + 
  theme_bw() + theme(panel.border = element_rect(colour = "black", size = 1.5), 
                     panel.grid.minor = element_blank(),
                     legend.position = "none",
                     panel.background = element_blank()) + scale_color_jama()

# vectorize the random forest selected taxa
ml_predictors <- c("f__f__CAG-508", "f__f__Enterobacteriaceae", "f__f__Streptococcaceae", "f__f__Treponemataceae", "f__f__CAG-302", "f__f__Lachnospiraceae")

set.seed(123)
vf <- envfit(sites, vector_data_species, perm = 99)
vf_scores <- as.data.frame(scores(vf, display = "vectors"))
vf_scores <- vf_scores[rownames(vf_scores) %in% ml_predictors, ]
vf_scores$species <- rownames(vf_scores)
figure_4c_1 <- ggplot(data = nmds.sites, aes(x = NMDS1, y = NMDS2)) + 
  labs(x = "NMDS1", y = "NMDS2", shape = "Site") +
  geom_point(inherit.aes = F, aes(x = NMDS1, y = NMDS2, color = cluster), size = 2.5, alpha = 0.3) + 
  geom_segment(data = vf_scores, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") +
  geom_text(data = vf_scores, aes(x = NMDS1, y = NMDS2, label = species), size = 3, position = position_dodge(width = 1)) +
  theme_bw(base_size = 14)+ theme(panel.grid = element_blank()) + scale_color_jama()


ml_predictors_box <- c("cluster", "f__f__CAG-508", "f__f__Enterobacteriaceae", "f__f__Streptococcaceae", 
                       "f__f__Treponemataceae", "f__f__CAG-302", "f__f__Lachnospiraceae")

figure_4c_2 <- vector_data %>% select(., any_of(ml_predictors_box)) %>% melt(., id.vars = "cluster") %>%
  ggplot() + aes(x = as.factor(cluster), y = log(value), group = as.factor(cluster)) + geom_boxplot(aes(fill = as.factor(cluster)), outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.2) +
  facet_wrap(.~ variable, ncol = 2, scales = "free_y") + theme_bw() + scale_fill_jama() + theme(panel.background = element_blank(), 
                                                                                                panel.grid.minor.y = element_blank(),
                                                                                                panel.grid.major.x = element_blank())
  
