## CAT- BAT Vector plots

## Purpose: ordinate AMR data and overlay vectors of CAT - BAT data

## setwd
setwd("/home")

## load packages
library(tidyverse)
library(vegan)
library(ggsci)

## load in the cat-bat files
cat_bat_family_otu <- readr::read_delim(file = "/home/cat_bat/cat_bat_family.csv")
cat_bat_family_otu <- cat_bat_family_otu %>% select(., -bacterium)

## load in the AMR files
## read in data
amr_genes <- read.csv("/home/datasets/amr_genes/FL100_merged_norm_final.csv", check.names = F)
abx_cluster <- read_delim("/home/datasets/from_andrew/abx_cluster_andrew.csv", delim = ",")[2:3]
abx_cluster <- abx_cluster %>% filter(., subject_id %in% cat_bat_family_otu$subject_id)
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
  rownames_to_column(., var = "subject_id") %>%
  filter(., subject_id %in% cat_bat_family_otu$subject_id)

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

## figure out what top mechanisms to plot
amr_genes_mech_melt <- amr_genes_mech_melt %>%
  mutate(., plot = ifelse(Mechanism %in% top_genes$Mechanism, Mechanism, "zlow_abundant"))
amr_genes_mech_melt$plot <- as.factor(amr_genes_mech_melt$plot)
levels(amr_genes_mech_melt$plot)
amr_genes_mech_melt <- amr_genes_mech_melt %>% group_by(., subject_id) %>% 
  mutate(., total_gene_abundance = sum(value))

sup_fig_1a_top <- ggplot(data = amr_genes_mech_melt, aes(x = reorder(as.factor(subject_id), -total_gene_abundance), y = as.numeric(value), group = plot)) +
  geom_bar(aes(fill = plot), width=1, stat = "identity") +
  coord_cartesian(ylim = c(0, 125), expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(axis.text.x.bottom = element_text(angle = 90), axis.ticks.x = element_blank(), panel.grid = element_blank()) + 
  labs(y = "Normalized Gene Abundance", x = "Individual") +
  scale_fill_d3()


## cat bat
order <- levels(reorder(as.factor(amr_genes_mech_melt$subject_id), -amr_genes_mech_melt$total_gene_abundance))
cat_bat_melt <- melt(id.vars = "subject_id", data = cat_bat_family_otu) %>%
  rename(., "Family" = "variable")
cat_bat_melt$Family <- as.character(cat_bat_melt$Family)

## get top amr mechanisms by cohort abundance
top_genes <- cat_bat_melt %>% 
  group_by(., Family) %>%
  summarise(., total = sum(value)) %>%
  arrange(., desc(total)) %>%
  slice(., (1:9)) %>%
  droplevels()

## figure out what top mechanisms to plot
cat_bat_melt <- cat_bat_melt %>%
  mutate(., plot = ifelse(Family %in% top_genes$Family, Family, "zlow_abundant"))
cat_bat_melt$plot <- as.factor(cat_bat_melt$plot)
levels(cat_bat_melt$plot)
cat_bat_melt <- cat_bat_melt %>% group_by(., subject_id) %>% 
  mutate(., total_gene_abundance = sum(value)) 
cat_bat_melt$subject_id <- factor(cat_bat_melt$subject_id, levels=order, ordered = T)
sup_fig_1a_bottom <- ggplot(data = cat_bat_melt, aes(x = as.factor(subject_id), weight = value, group = plot)) +
  geom_bar(aes(fill = plot), width=1) +
  labs(y = "contig count", x = "Individual") +
  #scale_y_continuous(labels = scales::percent_format()) + 
  #coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(axis.text.x.bottom = element_text(angle = 90), axis.ticks.x = element_blank(), panel.grid = element_blank()) + 
  ggsci::scale_fill_futurama()


