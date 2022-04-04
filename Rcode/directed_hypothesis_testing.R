###################################################################
# File: directed_hypothesis_testing.R                             #
#                                                                 #
# Purpose: After sourcing data wrangling script, assess directed  #
#          hypotheses statistics                                  #
#                                                                 #
#                                                                 #
# Author: A.Oliver				                                        #
# Date: 1/20/21						                                        #
#                                                                 #
# Inputs (1):                                                     #
# (1) Be able to source merging_features_new_data script          #
#      successfully                                               #
#                                                                 #
#                     #####################                       #
#                                                                 #
# Outputs (1):                                                    #    
# (1) data frame with directed hypotheses statistics              #
#                                                                 #
# Usage: Run the entire script without changes.                   #
###################################################################


#######################
## directed hypotheses
#######################

## load libraries
library(car)
library(patchwork)

## source or run merging_features_new_data.R to 
## generate merged features object (abx_cluster_features)
## abx_cluster_features is all the merged data, with NAs.

directed_hypotheses_features <- c("avg_fibe_tnfs", "per_kcal_fiber_tnfs", "total_fiber", "total_fiber_per_kcal",
                                  "dt_fiber_sol", "dt_fiber_sol_per_kcal", "dt_prot_animal", "dt_sfat", "dt_kcal",
                                  "pf_mps_total", "pf_mps_total_per_kcal", "pf_meat", "pf_meat_per_kcal", "hei_asa24_totalscore", "hei_ffq_totalscore",
                                  "fecal_ph", "bmi_final")

#source("/home/Rcode/merging_features_new_data.R")
## Select the features we care about and drop samples with NAs:
for_directed_hypothesis_testing <- abx_cluster_features %>%
  select(., subject_id, cluster, any_of(directed_hypotheses_features)) %>%
           drop_na() %>% droplevels()
## order the cluster factor
for_directed_hypothesis_testing$cluster <- factor(for_directed_hypothesis_testing$cluster, 
                                                  levels = c("low", "medium", "high"), ordered = TRUE)


## Extra function - SEM 
## thanks John https://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean
se <- function(x) sqrt(var(x)/length(x))
rm(Factor)
rm(data_source)
rm(resids)
rm(transformation_test)
#####################################
## 1. test average fiber intake ASA24
#####################################
## add factor information to dataframe
Factor <- c("Average fiber intake_tnfs, g")
data_source <- c("ASA24")
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$avg_fibe_tnfs, group = for_directed_hypothesis_testing$cluster)
## model function of avg fiber intake in response to AMR cluster
mod <- aov(avg_fibe_tnfs ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
set.seed(2022)
bestNormalize::bestNormalize(for_directed_hypothesis_testing$avg_fibe_tnfs)
## make new var that is now normalized based on previous data
for_directed_hypothesis_testing$avg_fibe_tnfs_norm <- (bestNormalize::sqrt_x(for_directed_hypothesis_testing$avg_fibe_tnfs))$x.t
## append information to dataframe
transformation_test <- c("sqrt(x), ANOVA, Tukey")
## re-model boxcox transformed avg fiber and AMR cluster
mod <- aov(avg_fibe_tnfs_norm ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
## double check homogeneity of variances assumption still met
leveneTest(for_directed_hypothesis_testing$avg_fibe_tnfs_norm, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- c(round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- c(round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- c(round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- c(round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(avg_fibe_tnfs),
                                                                       sd = sd(avg_fibe_tnfs),
                                                                       sem = se(avg_fibe_tnfs))
low_arg <- c(paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- c(paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- c(paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))

##################################################
## 2. Fiber, g/1000 kcal from food and supplements
##################################################
## add factor information to dataframe
Factor <- append(Factor, c("Fiber, g/1000 kcal from food and supplements"))
data_source <- append(data_source, c("ASA24"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$per_kcal_fiber_tnfs, group = for_directed_hypothesis_testing$cluster)
mod <- aov(per_kcal_fiber_tnfs ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
set.seed(2022)
bestNormalize::bestNormalize(for_directed_hypothesis_testing$per_kcal_fiber_tnfs)
for_directed_hypothesis_testing$per_kcal_fiber_tnfs_log <- (bestNormalize::yeojohnson(for_directed_hypothesis_testing$per_kcal_fiber_tnfs))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("yeojohnson, ANOVA, Tukey"))
mod <- aov(per_kcal_fiber_tnfs_log ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$per_kcal_fiber_tnfs_log, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(per_kcal_fiber_tnfs),
                                                                       sd = sd(per_kcal_fiber_tnfs),
                                                                       sem = se(per_kcal_fiber_tnfs))

low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))
############################
## 3. Total fiber intake, g
###########################
## add factor information to dataframe
Factor <- append(Factor, c("Total fiber intake, g"))
data_source <- append(data_source, c("FFQ"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$total_fiber, group = for_directed_hypothesis_testing$cluster)
mod <- aov(total_fiber ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
set.seed(2022)
bestNormalize::bestNormalize(for_directed_hypothesis_testing$total_fiber)
for_directed_hypothesis_testing$total_fiber_log <- (bestNormalize::yeojohnson(for_directed_hypothesis_testing$total_fiber))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("yeojohnson, ANOVA, Tukey"))
mod <- aov(total_fiber_log ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$total_fiber_log, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(total_fiber),
                                                                       sd = sd(total_fiber),
                                                                       sem = se(total_fiber))
low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))

##############################
## 4. Total fiber, g/1000 kcal 
#############################
## add factor information to dataframe
Factor <- append(Factor, c("Total fiber, g/1000 kcal"))
data_source <- append(data_source, c("FFQ"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$total_fiber_per_kcal, group = for_directed_hypothesis_testing$cluster)
mod <- aov(total_fiber_per_kcal ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
set.seed(2022)
bestNormalize::bestNormalize(for_directed_hypothesis_testing$total_fiber_per_kcal)
for_directed_hypothesis_testing$total_fiber_per_kcal_log <- (bestNormalize::arcsinh_x(for_directed_hypothesis_testing$total_fiber_per_kcal))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("arcsinh, ANOVA, Tukey"))
mod <- aov(total_fiber_per_kcal_log ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$total_fiber_per_kcal_log, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(total_fiber_per_kcal),
                                                                       sd = sd(total_fiber_per_kcal),
                                                                       sem = se(total_fiber_per_kcal))
low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))


##############################
## 5. Soluble fiber intake, g 
#############################
## add factor information to dataframe
Factor <- append(Factor, c("Soluble fiber intake, g"))
data_source <- append(data_source, c("FFQ"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$dt_fiber_sol, group = for_directed_hypothesis_testing$cluster)
mod <- aov(dt_fiber_sol ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
set.seed(2022)
bestNormalize::bestNormalize(for_directed_hypothesis_testing$dt_fiber_sol)
for_directed_hypothesis_testing$dt_fiber_sol_log <- (bestNormalize::yeojohnson(for_directed_hypothesis_testing$dt_fiber_sol))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("yeojohnson, ANOVA, Tukey"))
mod <- aov(dt_fiber_sol_log ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$dt_fiber_sol_log, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_fiber_sol),
                                                                       sd = sd(dt_fiber_sol),
                                                                       sem = se(dt_fiber_sol))
low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))

#####################################################
## 6. Soluble fiber intake, g/1000 kcal energy intake
#####################################################
## add factor information to dataframe
Factor <- append(Factor, c("Soluble fiber intake, g/1000 kcal energy intake"))
data_source <- append(data_source, c("FFQ"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$dt_fiber_sol_per_kcal, group = for_directed_hypothesis_testing$cluster)
mod <- aov(dt_fiber_sol_per_kcal ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
set.seed(2022)
bestNormalize::bestNormalize(for_directed_hypothesis_testing$dt_fiber_sol_per_kcal)
for_directed_hypothesis_testing$dt_fiber_sol_per_kcal_boxcox <- (bestNormalize::log_x(for_directed_hypothesis_testing$dt_fiber_sol_per_kcal))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("Log10(x), ANOVA, Tukey"))
mod <- aov(dt_fiber_sol_per_kcal_boxcox ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$dt_fiber_sol_per_kcal_boxcox, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_fiber_sol_per_kcal),
                                                                       sd = sd(dt_fiber_sol_per_kcal),
                                                                       sem = se(dt_fiber_sol_per_kcal))
low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))

##############################################
## 7. Habitual protein intake from animals, g
##############################################
## add factor information to dataframe
Factor <- append(Factor, c("Habitual protein intake from animals, g"))
data_source <- append(data_source, c("FFQ"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$dt_prot_animal, group = for_directed_hypothesis_testing$cluster)
mod <- aov(dt_prot_animal ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
set.seed(2022)
bestNormalize::bestNormalize(for_directed_hypothesis_testing$dt_prot_animal)
for_directed_hypothesis_testing$dt_prot_animal_log<- (bestNormalize::arcsinh_x(for_directed_hypothesis_testing$dt_prot_animal))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("arcsinh, ANOVA, Tukey"))
mod <- aov(dt_prot_animal_log ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$dt_prot_animal_log, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_prot_animal),
                                                                       sd = sd(dt_prot_animal),
                                                                       sem = se(dt_prot_animal))
low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))


#####################################
## 8. Saturated fatty acid intake, g
#####################################
## add factor information to dataframe
Factor <- append(Factor, c("Saturated fatty acid intake, g"))
data_source <- append(data_source, c("FFQ"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$dt_sfat, group = for_directed_hypothesis_testing$cluster)
mod <- aov(dt_sfat ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
set.seed(2022)
bestNormalize::bestNormalize(for_directed_hypothesis_testing$dt_sfat)
for_directed_hypothesis_testing$dt_sfat_boxcox<- (bestNormalize::boxcox(for_directed_hypothesis_testing$dt_sfat))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("BoxCox, ANOVA, Tukey"))
mod <- aov(dt_sfat_boxcox ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$dt_sfat_boxcox, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_sfat),
                                                                       sd = sd(dt_sfat),
                                                                       sem = se(dt_sfat))
low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))

###############################
## 9. Total energy intake, kcal
###############################
## add factor information to dataframe
Factor <- append(Factor, c("Total energy intake, kcal"))
data_source <- append(data_source, c("FFQ"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$dt_kcal, group = for_directed_hypothesis_testing$cluster)
mod <- aov(dt_kcal ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
set.seed(2022)
bestNormalize::bestNormalize(for_directed_hypothesis_testing$dt_kcal)
for_directed_hypothesis_testing$dt_kcal_boxcox<- (bestNormalize::arcsinh_x(for_directed_hypothesis_testing$dt_kcal))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("arcsinh, ANOVA, Tukey"))
mod <- aov(dt_kcal_boxcox ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$dt_kcal_boxcox, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_kcal),
                                                                       sd = sd(dt_kcal),
                                                                       sem = se(dt_kcal))
low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))

########################################################
## 10. Beef/pork consumption including cured meat, oz eq
########################################################
## add factor information to dataframe
Factor <- append(Factor, c("Beef/pork consumption including cured meat, oz eq"))
data_source <- append(data_source, c("FFQ"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$pf_mps_total, group = for_directed_hypothesis_testing$cluster)
mod <- aov(pf_mps_total ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
set.seed(2022)
bestNormalize::bestNormalize(for_directed_hypothesis_testing$pf_mps_total)
for_directed_hypothesis_testing$pf_mps_total_boxcox<- (bestNormalize::arcsinh_x(for_directed_hypothesis_testing$pf_mps_total))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("arcsinh, ANOVA, Tukey"))
mod <- aov(pf_mps_total_boxcox ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$pf_mps_total_boxcox, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(pf_mps_total),
                                                                       sd = sd(pf_mps_total),
                                                                       sem = se(pf_mps_total))
low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))

####################################################################
## 11. Beef/pork consumption including cured meat,  oz eq /1000 kcal
####################################################################
## add factor information to dataframe
Factor <- append(Factor, c("Beef/pork consumption including cured meat,  oz eq /1000 kcal"))
data_source <- append(data_source, c("FFQ"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$pf_mps_total_per_kcal, group = for_directed_hypothesis_testing$cluster)
mod <- aov(pf_mps_total_per_kcal ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
set.seed(2022)
bestNormalize::bestNormalize(for_directed_hypothesis_testing$pf_mps_total_per_kcal)
for_directed_hypothesis_testing$pf_mps_total_per_kcal_boxcox<- (bestNormalize::yeojohnson(for_directed_hypothesis_testing$pf_mps_total_per_kcal))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("yeojohnson, ANOVA, Tukey"))
mod <- aov(pf_mps_total_per_kcal_boxcox ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$pf_mps_total_per_kcal_boxcox, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(pf_mps_total_per_kcal),
                                                                       sd = sd(pf_mps_total_per_kcal),
                                                                       sem = se(pf_mps_total_per_kcal))
low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))

#########################################################
## 12. Beef/pork consumption excluding cured meat,  oz eq
#########################################################
## add factor information to dataframe
Factor <- append(Factor, c("Beef/pork consumption excluding cured meat,  oz eq"))
data_source <- append(data_source, c("FFQ"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$pf_meat, group = for_directed_hypothesis_testing$cluster)
mod <- aov(pf_meat ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
set.seed(2022)
bestNormalize::bestNormalize(for_directed_hypothesis_testing$pf_meat)
for_directed_hypothesis_testing$pf_meat_log<- (bestNormalize::boxcox(for_directed_hypothesis_testing$pf_meat))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("boxcox, ANOVA, Tukey"))
mod <- aov(pf_meat_log ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$pf_meat_log, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(pf_meat),
                                                                       sd = sd(pf_meat),
                                                                       sem = se(pf_meat))
low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))

####################################################################
## 13. Beef/pork consumption excluding cured meat,  oz eq /1000 kcal
####################################################################
## add factor information to dataframe
Factor <- append(Factor, c("Beef/pork consumption excluding cured meat,  oz eq /1000 kcal"))
data_source <- append(data_source, c("FFQ"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$pf_meat_per_kcal, group = for_directed_hypothesis_testing$cluster)
mod <- aov(pf_meat_per_kcal ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
set.seed(2022)
bestNormalize::bestNormalize(for_directed_hypothesis_testing$pf_meat_per_kcal)
for_directed_hypothesis_testing$pf_meat_per_kcal_boxcox<- (bestNormalize::yeojohnson(for_directed_hypothesis_testing$pf_meat_per_kcal))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("yeojohnson, ANOVA, Tukey"))
mod <- aov(pf_meat_per_kcal_boxcox ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$pf_meat_per_kcal_boxcox, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(pf_meat_per_kcal),
                                                                       sd = sd(pf_meat_per_kcal),
                                                                       sem = se(pf_meat_per_kcal))
low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))

#######################
## 14. HEI total score
#######################
## add factor information to dataframe
Factor <- append(Factor, c("HEI total score"))
data_source <- append(data_source, c("ASA24"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$hei_asa24_totalscore, group = for_directed_hypothesis_testing$cluster)
mod <- aov(hei_asa24_totalscore ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## many good options, no transformation performed
transformation_test <- append(transformation_test,c("None, ANOVA, Tukey"))
mod <- aov(hei_asa24_totalscore ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$hei_asa24_totalscore, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(hei_asa24_totalscore),
                                                                       sd = sd(hei_asa24_totalscore),
                                                                       sem = se(hei_asa24_totalscore))
low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))

#######################
## 15. Stool pH
#######################
## add factor information to dataframe
Factor <- append(Factor, c("Stool pH"))
data_source <- append(data_source, c("Not applicable"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$fecal_ph, group = for_directed_hypothesis_testing$cluster)
mod <- aov(fecal_ph ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## many good options, no transformation performed
transformation_test <- append(transformation_test,c("None, ANOVA, Tukey"))
mod <- aov(fecal_ph ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$fecal_ph, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(fecal_ph),
                                                                       sd = sd(fecal_ph),
                                                                       sem = se(fecal_ph))
low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))

#######################
## 16. BMI
#######################
## add factor information to dataframe
Factor <- append(Factor, c("BMI"))
data_source <- append(data_source, c("Not applicable"))
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$bmi_final, group = for_directed_hypothesis_testing$cluster)
mod <- aov(bmi_final ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
set.seed(2022)
bestNormalize::bestNormalize(for_directed_hypothesis_testing$bmi_final)
for_directed_hypothesis_testing$bmi_final_boxcox<- (bestNormalize::yeojohnson(for_directed_hypothesis_testing$bmi_final))$x.t
## many good options, no transformation performed
transformation_test <- append(transformation_test,c("yeojohnson, ANOVA, Tukey"))
mod <- aov(bmi_final_boxcox ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$bmi_final_boxcox, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
tmp2 <- summary(mod)
p_all <- append(p_all, round(tmp2[[1]][["Pr(>F)"]][1], digits = 3))
tmp3 <- TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, round(tmp3$cluster[1,4], digits = 3))
p_low_high_tukey <- append(p_low_high_tukey, round(tmp3$cluster[2,4], digits = 3))
p_med_high_tukey <- append(p_med_high_tukey, round(tmp3$cluster[3,4], digits = 3))
## describe
tmp4 <- for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(bmi_final),
                                                                       sd = sd(bmi_final),
                                                                       sem = se(bmi_final))
low_arg <- append(low_arg, paste0(round(tmp4[1,2],digits = 2), " +/- ", round(tmp4[1,3],digits = 2), " (", round(tmp4[1,4],digits = 2), ")"))
medium_arg <- append(medium_arg, paste0(round(tmp4[2,2],digits = 2), " +/- ", round(tmp4[2,3],digits = 2), " (", round(tmp4[2,4],digits = 2), ")"))
high_arg <- append(high_arg, paste0(round(tmp4[3,2],digits = 2), " +/- ", round(tmp4[3,3],digits = 2), " (", round(tmp4[3,4],digits = 2), ")"))

##################
## Complete the DF
##################

directed_hypothesis_table <- data.frame(Factor, data_source, 
           low_arg, medium_arg, high_arg, 
           transformation_test, p_low_med_tukey, p_low_high_tukey,p_med_high_tukey, p_all)

View(directed_hypothesis_table)



