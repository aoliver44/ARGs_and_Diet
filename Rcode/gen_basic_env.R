###################################################################
# File: gen_basic_env.R                                           #
#                                                                 #
# Purpose: Generate rarefied microbiome datasets and microbiome   #
#          data grouped by family and genus levels                #
#                                                                 #
#                                                                 #
# Author: A.Oliver				                                        #
# Date: 1/24/21						                                        #
#                                                                 #
# Inputs (2):                                                     #
# (1) MPA formated combined Kraken outputs (combined_MPA.txt)     #
# (2) ABX cluster metadata by sample                              #
#                                                                 #                                                                #
#                     #####################                       #
#                                                                 #
# Outputs (2):                                                    #    
# (1) mpa_rare_5833371_perm_5.rds R data structure rarefied table #
# (1) mpa_rare_5833371_perm_5_family.csv at family level          #
#                                                                 #
# Usage: Run the entire script without changes.                   #
###################################################################


#######################
#### GEN BASIC ENV ####
#######################

set.seed(999)
# load libraries:
library(vegan)
library(spaa)
library(tidyverse)
library(EcolUtils)
library(reshape2)

setwd("/home/Kracken/")

###############
## READ IN DATA 
###############

# read in kracken (converted to mpa format) data
mpa_raw <- read_delim("combined_MPA.txt", delim = "\t")
# remove _report from the col names
names(mpa_raw) <- gsub("_report", "", names(mpa_raw))
colnames(mpa_raw)[1] <- "Classification"

# get rid of everything that is not a species
mpa_species <- mpa_raw %>% filter(., grepl("s__",Classification)) %>%
  column_to_rownames(., var = "Classification") %>% t()

# read in metadata clusters
abx_clusters <- read_csv("../datasets/from_andrew/abx_cluster_andrew.csv")[2:3]

##############
## RAREFACTION
##############

# get rid of all the species with <1 reads (kracken was run to keep 0 hits in)
mpa_species <- mpa_species[,colSums(mpa_species) >= 1 ]

# Check read distribution
min_reads <- min(rowSums(mpa_species))
max_reads <- max(rowSums(mpa_species))
barplot(sort(rowSums(mpa_species)), main = "Number of reads per sample", 
        sub = paste("Min read depth:",  min_reads, "\nMax read depth:", 
                    max_reads,sep=" "))
# Rarefy
#mpa_rare <- rrarefy.perm(mpa_species, sample = 5833371, n = 5, round.out = T)
#saveRDS(mpa_rare, file = "mpa_rare_5833371_perm_5.rds", version = 2)
mpa_rare <- readRDS("mpa_rare_5833371_perm_5.rds")

# Check read distribution post rarefaction
min_reads_post_rare <- min(rowSums(mpa_rare))
max_reads_post_rare <- max(rowSums(mpa_rare))
barplot(sort(rowSums(mpa_rare)), main = "Number of reads per sample post rarefaction", 
        sub = paste("Min read depth:",  min_reads_post_rare, "\nMax read depth:", 
                    max_reads_post_rare,sep=" "))

# check goods coverage:
# thanks to John Quensen for his function https://rdrr.io/github/jfq3/QsRutils/man/goods.html
goods <-
  function(com){
    no.seqs <- rowSums(com)
    sing <- com==1
    no.sing <- apply(sing, 1, sum)
    goods <- 100*(1-no.sing/no.seqs)
    goods.sum <- cbind(no.sing, no.seqs, goods)
    goods.sum <- as.data.frame(goods.sum)
    return(goods.sum)
  }

goods_analysis <- goods(mpa_rare) # every sample has like 99.99% coverage
##################
## ALPHA DIVERSITY
##################

# tidy up data
shannon_div <- diversity(mpa_rare, "shannon") 
richness <- specnumber(mpa_rare)
Evenness <- shannon_div/log(richness)
all_diversity <- as.data.frame(cbind(Evenness, richness, shannon_div))

# Plot diversity:
ggplot(all_diversity, aes(x = rownames(all_diversity), y = shannon_div)) +
  geom_point()

#####################
## ENTEROBACTERIACEAE
#####################

# seperate taxa name by taxonomic level
mpa_melt <- as.data.frame(t(mpa_rare)) %>%
  rownames_to_column(., var = "taxa_name") %>%
  separate(., col = taxa_name, 
           into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), 
           sep = "\\|", 
           remove = T, 
           extra = "drop") %>%
  melt()

# group by individual, Enterobacteriaceae

enterobacteriaceae <- mpa_melt %>% 
  filter(., Family == "f__f__Enterobacteriaceae") %>%
  group_by(., Family, variable) %>%
  summarise(., total_enterobacteriaceae = sum(value))

# plot enterobacteriaceae

ggplot(data = enterobacteriaceae, aes(x = reorder(as.factor(variable), total_enterobacteriaceae), y = log(total_enterobacteriaceae))) +
  geom_bar(stat = "identity") + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90))

########################
## COMBINE WITH METADATA
########################

merged_alpha <- merge(abx_clusters, all_diversity, by.x = "subject_id", by.y = "row.names")
merged_all <- merge(merged_alpha, enterobacteriaceae, by.x = "subject_id", by.y = "variable")
merged_all$cluster <- factor(merged_all$cluster, levels = c("low", "medium", "high"), ordered = T)

# plot enterobacteriaceae by abx cluster
ggplot(data = merged_all, aes(x = cluster , y = log2(total_enterobacteriaceae))) +
  geom_boxplot(aes(fill = cluster), outlier.alpha = 0) + geom_point(position = position_jitter(width = 0.1)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Antibiotic Gene Resistance Cluster",
       y = "Enterobacteriaeace abundance \n(rarefied, log2 transformed)") +
  theme_bw()

# plot cluster by alpha diversity
ggplot(data = merged_all, aes(x = cluster, y = shannon_div)) +
  geom_boxplot(aes(fill = cluster), outlier.alpha = 0) + 
  geom_point(position = position_jitter(width = 0.1)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Antibiotic Gene Resistance Cluster",
       y = "Shannon Diversity (rarefied)") +
  theme_bw()

write.csv(merged_all, file = "alpha_entero.csv", quote = F, row.names = F)

########################################
## GROUP MICROBIOME DATA BY HIGHER TAXON
########################################

mpa_taxa_group <- as.data.frame(t(mpa_rare)) %>%
  rownames_to_column(., var = "taxa_name") %>%
  separate(., col = taxa_name, 
           into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), 
           sep = "\\|", 
           remove = T, 
           extra = "drop") %>% melt()

mpa_taxa_group <- mpa_taxa_group %>%
  select(., variable, Family, value) %>%
  group_by(., Family, variable) %>%
  summarise(., reads = sum(value))

mpa_taxa_group_otu <- dcast(formula = variable ~ Family, data = mpa_taxa_group)
mpa_taxa_group_otu <- merge(abx_clusters, mpa_taxa_group_otu, by.x = "subject_id", by.y = "variable")
mpa_taxa_group_otu <- mpa_taxa_group_otu %>% 
  select(., -subject_id) %>%
  dplyr::mutate(., cluster = ifelse(cluster == "low", 0, 
                                    ifelse(cluster == "medium", 1, 2))) %>% 
  drop_na() %>% 
  clean_names()

#write.csv(mpa_taxa_group_otu, file = "/home/datasets/from_andrew/mpa_rare_5833371_perm_5_family.csv", quote = F, row.names = F)



