############################
## MERGE ALL DATA - NEW DATA
############################

library(tidyverse)
library(Hmisc)
library(janitor)
library(mikropml)
library(naniar)
setwd("/home/datasets/new_datasets/")

features <- list()

## START WITH ANDREW'S CLUSTERS
abx_cluster <- read_delim("abx_cluster_andrew.csv", delim = ",") 
abx_cluster <- abx_cluster %>% select(., -`...1`) 
features <- append(x = features, values = colnames(abx_cluster))

## diversity of fiber foods in diet
fiber_alpha <- read_delim("qiime1_alphadiv_fiber.txt", delim = "\t")
fiber_alpha <- fiber_alpha %>% 
  rename(., "subject_id" = `...1`, "PD_whole_tree_fiber" = "PD_whole_tree") %>% 
  select(., subject_id, PD_whole_tree_fiber) %>% clean_names()
features <- append(x = features, values = colnames(fiber_alpha))

## diversity of fat foods in diet
fat_alpha <- read_delim("qiime1_alphadiv_fat.txt", delim = "\t")
fat_alpha <- fat_alpha %>% 
  rename(., "subject_id" = `...1`, "PD_whole_tree_fat" = "PD_whole_tree") %>% 
  select(., subject_id, PD_whole_tree_fat) %>% clean_names()
features <- append(x = features, values = colnames(fat_alpha))

## diversity of protein foods in diet
protein_alpha <- read_delim("qiime1_alphadiv_protein.txt", delim = "\t")
protein_alpha <- protein_alpha %>% 
  rename(., "subject_id" = `...1`, "PD_whole_tree_protein" = "PD_whole_tree") %>% 
  select(., subject_id, PD_whole_tree_protein) %>% clean_names()
features <- append(x = features, values = colnames(protein_alpha))

## diversity of carb foods in diet
carb_alpha <- read_delim("qiime1_alphadiv_carb.txt", delim = "\t")
carb_alpha <- carb_alpha %>% 
  rename(., "subject_id" = `...1`, "PD_whole_tree_carb" = "PD_whole_tree") %>% 
  select(., subject_id, PD_whole_tree_carb) %>% clean_names()
features <- append(x = features, values = colnames(carb_alpha))

## diversity of all foods in diet
total_food_alpha <- read_delim("qiime1_alphadiv_filt.txt", delim = "\t")
total_food_alpha <- total_food_alpha %>% 
  rename(., "subject_id" = `...1`, "PD_whole_tree_total" = "PD_whole_tree") %>% 
  select(., subject_id, PD_whole_tree_total) %>% clean_names()
features <- append(x = features, values = colnames(total_food_alpha))

## ffq from files sent to Blythe from danielle
ffq_data <- read_delim("FL100_FFQ_cleaned_all_dt.csv", delim = ",") %>% clean_names()
features <- append(x = features, values = colnames(ffq_data))

## lifestyle features
lifestyle_data <- read_delim("report_for_DL_7262021_updated_AO.csv", delim = ",")
lifestyle_data <- rename(lifestyle_data, "subject_id" = "Participant ID") %>% clean_names()
features <- append(x = features, values = colnames(lifestyle_data))

## stool data and fecal cytokines
stool_variables <- read_delim("FL100_stool_variables.txt", delim = "\t")
stool_variables <- stool_variables %>% 
  filter(., AfterV2 != 1) %>% 
  filter(., diff_time_hrs < 24) %>%
  select(., -c("visit2_date", "stool_collection", "diff_time_hrs", "stool_specimen", "AfterV2")) %>% 
  clean_names()
features <- append(x = features, values = colnames(stool_variables))

## country of origin
country_origin <- readxl::read_xlsx("CountryOrigin.xlsx")
country_origin <- rename(country_origin, "subject_id" = "subject_id") %>% 
  mutate(., birth_country = ifelse(birth_country == "United States", "united_states", "non_us")) %>% 
  clean_names()
features <- append(x = features, values = colnames(country_origin))

## more dietary features (like ASA24) from elizabeth
stoolc_model2 <- read_delim("stoolc_model2_features.csv", delim = ",")
stoolc_model2 <- rename(stoolc_model2, "subject_id" = "SubjectID") %>% clean_names() %>% relocate(., subject_id)
stoolc_model3 <- read_delim("stoolc_model3_features.csv", delim = ",")
stoolc_model3 <- rename(stoolc_model3, "subject_id" = "SubjectID") %>% clean_names() %>% relocate(., subject_id)
#redundant_vars <- intersect(colnames(stoolc_model2), colnames(stoolc_model3)) %>% discard(~ .x %in% c("subject_id"))
#stoolc_model2 <- stoolc_model2 %>%  select(., -any_of(redundant_vars)) 
stoolc_model23 <- merge(stoolc_model2, stoolc_model3, by = "subject_id") %>% select(., -fibe_total)
features <- append(x = features, values = colnames(stoolc_model23))

## inflammatory markers from Charles
inflammatory_markers <- read_delim("CRP_WBC_9102021.csv", delim = ",")
inflammatory_markers <- inflammatory_markers %>% clean_names()
features <- append(x = features, values = colnames(inflammatory_markers))

## average ASA24 fiber intake from Blythe (danielle sent her)
asa_fiber_avg <- read_delim("ASA24_average_fiber_summary_variables.txt", delim = "\t") %>% clean_names()
features <- append(x = features, values = colnames(asa_fiber_avg))

##  ASA24 HEI  from Blythe (danielle sent her)
asa_hei <- read_delim("FL100_HEI_n378.txt", delim = "\t") %>% clean_names()
features <- append(x = features, values = colnames(asa_hei))

## Dexa ethnicities sent by elizabeth
ethnicities <- read_delim("DEXA_ethnicities04272020.csv", delim = ",") %>% clean_names()
features <- append(x = features, values = colnames(ethnicities))

## BMI/sex/age sent to Blythe from danielle
bmi_etc <- read_delim("FL100_FINAL_Bins_agesexbmi_clean.csv", delim = ",") %>% clean_names()
features <- append(x = features, values = colnames(bmi_etc))

## vitals from danielle from redcap (email to Jules)
vitals <- read_delim("CTSC24532USDAWHNRCNu-VitalsPhysiologyforJ_DATA_2021-10-07_1622.csv", delim = ",")
vitals <- vitals %>% select(., -c(redcap_survey_identifier, vital_signs_timestamp, vital_comments, vital_signs_complete)) %>% clean_names()
features <- append(x = features, values = colnames(vitals))

## plasma lbd from danielle from redcap (email to Jules)
inflammation_markers_2 <- read_delim("CTSC24532USDAWHNRCNu-GIMarkers7Oct2021_DATA_2021-10-07_1627.csv", delim = ",")
inflammation_markers_2 <- inflammation_markers_2 %>% select(., subject_id, plasma_lbp_bd1) %>% clean_names()
features <- append(x = features, values = colnames(inflammation_markers_2))

## FFQ HEI index
ffq_hei <- read_delim("HEI FFQ_scores_12072021.csv", delim = ",")
ffq_hei <- ffq_hei %>% select(., subject_id, hei_ffq_totalscore) %>% clean_names()
features <- append(x = features, values = colnames(ffq_hei))

features <- unique(unlist(features))

abx_cluster_features <- merge(fiber_alpha, abx_cluster, by = "subject_id") %>%
  merge(fat_alpha, by = "subject_id") %>%
  merge(protein_alpha, by = "subject_id") %>%
  merge(carb_alpha, by = "subject_id") %>%
  merge(total_food_alpha, by = "subject_id") %>%
  merge(lifestyle_data, by = "subject_id") %>%
  merge(country_origin, by = "subject_id") %>%
  merge(stool_variables, by = "subject_id") %>%
  merge(ffq_data, by = "subject_id") %>%
  #merge(micro_alpha, by = "subject_id") %>%
  merge(inflammatory_markers, by = "subject_id") %>%
  merge(inflammation_markers_2, by = "subject_id") %>%
  merge(bmi_etc, by = "subject_id") %>%
  merge(asa_fiber_avg, by = "subject_id") %>%
  merge(asa_hei, by = "subject_id") %>%
  merge(stoolc_model23, by = "subject_id") %>%
  merge(vitals, by = "subject_id") %>%
  merge(ffq_hei, by = "subject_id") %>%
  merge(ethnicities, by = "subject_id")
  

colnames(abx_cluster_features) = gsub("\\.y", "", colnames(abx_cluster_features))
colnames(abx_cluster_features) = gsub("\\.x", "", colnames(abx_cluster_features))
abx_cluster_features <- subset(abx_cluster_features, select=which(!duplicated(names(abx_cluster_features)))) 


## clean up NAs
gg_miss_upset(abx_cluster_features) 
# clean up data based on results - drop some samples or features
abx_cluster_features <- abx_cluster_features %>%
  filter(., stool_consistency_class != "NA") %>%
  filter(., wbc_bd1 != "NA") %>%
  select(., -tempavg) %>%
  select(., -fecal_neopterin) %>%
  select(., -waist_hip)

# check the missing values again
gg_miss_upset(abx_cluster_features) 
abx_cluster_features <- abx_cluster_features %>%
  select(., -c(do_you_ever_use_tanning_beds_or_sun_lamps, dov1))
gg_miss_upset(abx_cluster_features) 

## output dataframe for directed hypothesis testing
abx_cluster_features$total_fiber_per_kcal <- ((abx_cluster_features$total_fiber / abx_cluster_features$dt_kcal) * 1000)
abx_cluster_features$dt_fiber_sol_per_kcal <- ((abx_cluster_features$dt_fiber_sol / abx_cluster_features$dt_kcal) * 1000)
abx_cluster_features$pf_mps_total_per_kcal <- (abx_cluster_features$pf_mps_total / abx_cluster_features$dt_kcal) * 1000
abx_cluster_features$pf_meat_per_kcal <- (abx_cluster_features$pf_meat / abx_cluster_features$dt_kcal) * 1000
for_directed_hypothesis_testing <- abx_cluster_features

# rename the cluster column!
abx_cluster_features_no_na <- abx_cluster_features %>% 
  dplyr::mutate(., cluster = ifelse(cluster == "low", 0, 
                                 ifelse(cluster == "medium", 1, 2))) %>% drop_na()

#write.csv(abx_cluster_features, file = "abx_cluster_features.csv", quote = F, row.names = F)
#write.csv(abx_cluster_features_no_na, file = "abx_cluster_features_no_na.csv", quote = F, row.names = F)

###########################################
## MIKROPML - COLLAPSE CORRELATED VARIABLES
###########################################

## Clean up and preprocess with mikropml
abx_cluster_features_pre_process <- abx_cluster_features_no_na %>% droplevels()
abx_cluster_features_post_process <- preprocess_data(dataset = abx_cluster_features_pre_process,
                                                     method = NULL,
                                                     outcome_colname = "cluster", 
                                                     collapse_corr_feats = T, 
                                                     remove_var = "zv")

pre_corr_clean <- abx_cluster_features_post_process$dat_transformed %>% 
  # pre-remove correlated variables to sex, age, bmi that end up in
  # the feature-set before these variables do (these are important covariates!)
  dplyr::select(., -bin_number, -age_cat, -bmi_cat, -bmi, -hipavg_cm, 
                -screen_endmenstr, -screen_contracept, -menstr_contracept)

## get the correlated features to a dataframe
cor_tmp <- mikropml:::group_correlated_features(pre_corr_clean, 
                                                corr_thresh = 0.80)
cor_tmp <- as.data.frame(cor_tmp)
cor_tmp <- cor_tmp %>% separate(., col = cor_tmp, into = c("keep", "co-correlated"), sep = "\\|", extra = "merge")
## filter out only 1 of the co-correlated groups, results in 146 features (including subject ID)
corr_clean <- pre_corr_clean %>% select(., cor_tmp$keep)
## write to file
corr_clean_write <- corr_clean %>% select(., -subject_id) %>% clean_names()
## All clusters - 187 samples 
write.csv(corr_clean_write, file = "output_for_ML/diet-life_input.csv", quote = F, row.names = FALSE)
write.csv(cor_tmp, file = "output_for_ML/co-correlated-features_spearman80_covariates.csv", quote = F, row.names = FALSE)
## Just low-medium clusters - 140 samples
corr_clean_bi <- corr_clean_write %>% filter(., cluster != 2)
write.csv(corr_clean_bi, file = "output_for_ML/diet-life_input_bi.csv", quote = F, row.names = F)
## Just low-high clusters - 92 samples
corr_clean_low_high <- corr_clean_write %>% filter(., cluster != 1) %>% 
  mutate(., cluster = ifelse(cluster == 2, 1, 0))
write.csv(corr_clean_low_high, file = "output_for_ML/diet-life_input_low-high.csv", quote = F, row.names = F)

#############
## MICROBIOME
#############

## read in rarefied microbiome - ***family level***
microbiome_family <- read_delim("/home/datasets/from_andrew/mpa_rare_5833371_perm_5_family.csv", delim = ",")
microbiome_family <- mpa_taxa_group_otu
## get the correlated features to a dataframe
microbiome_features_post_process <- preprocess_data(dataset = microbiome_family,
                                                     method = NULL,
                                                     outcome_colname = "cluster", 
                                                     collapse_corr_feats = F, 
                                                     remove_var = "zv")
cor_tmp <- mikropml:::group_correlated_features(microbiome_features_post_process$dat_transformed, 
                                                corr_thresh = 0.80)
cor_tmp <- as.data.frame(cor_tmp)
cor_tmp <- cor_tmp %>% separate(., col = cor_tmp, into = c("keep", "co-correlated"), sep = "\\|", extra = "merge")
microbiome_family_clean <- microbiome_features_post_process$dat_transformed %>% select(., cor_tmp$keep)
write.csv(microbiome_family_clean, file = "/home/datasets/new_datasets/output_for_ml/microbiome_family.csv", quote = F, row.names = F)


## read in rarefied microbiome - ***genus level***
microbiome_genus <- read_delim("/home/datasets/from_andrew/mpa_rare_5833371_perm_5_genus.csv", delim = ",")
## get the correlated features to a dataframe
microbiome_features_post_process <- preprocess_data(dataset = microbiome_genus,
                                                    method = NULL,
                                                    outcome_colname = "cluster", 
                                                    collapse_corr_feats = F, 
                                                    remove_var = "zv")
cor_tmp <- mikropml:::group_correlated_features(microbiome_features_post_process$dat_transformed, 
                                                corr_thresh = 0.80)
cor_tmp <- as.data.frame(cor_tmp)
cor_tmp <- cor_tmp %>% separate(., col = cor_tmp, into = c("keep", "co-correlated"), sep = "\\|", extra = "merge")
microbiome_genus_clean <- microbiome_features_post_process$dat_transformed %>% select(., cor_tmp$keep)
write.csv(microbiome_genus_clean, file = "/home/datasets/new_datasets/output_for_ml/microbiome_genus.csv", quote = F, row.names = F)

