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
#source("/home/scripts/merging_features_new_data.R")

## Select the features we care about and drop samples with NAs:
for_directed_hypothesis_testing <- abx_cluster_features %>%
  select(., subject_id, cluster, avg_fibe_tnfs, per_kcal_fiber_tnfs, total_fiber, total_fiber_per_kcal,
         dt_fiber_sol, dt_fiber_sol_per_kcal, dt_prot_animal, dt_sfat, dt_kcal,
         pf_mps_total, pf_mps_total_per_kcal, pf_meat, pf_meat_per_kcal, hei_asa24_totalscore, hei_ffq_totalscore,
         fecal_ph, bmi_final) %>% drop_na()
## order the cluster factor
for_directed_hypothesis_testing$cluster <- factor(for_directed_hypothesis_testing$cluster, 
                                                  levels = c("low", "medium", "high"), ordered = TRUE)


## Extra function - SEM 
## thanks John https://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean
se <- function(x) sqrt(var(x)/length(x))

#####################################
## 1. test average fiber intake ASA24
#####################################
## add factor information to dataframe
Factor <- c("Average fiber intake_tnfs, g")
data_source <- c("ASA24")
## test homogeneity of sample variances
leveneTest(for_directed_hypothesis_testing$avg_fibe_tnfs, group = for_directed_hypothesis_testing$cluster)
mod <- aov(avg_fibe_tnfs ~ cluster, data = for_directed_hypothesis_testing)
## test normality of residuals
resids <- residuals(mod)
shapiro.test(resids)
## transform to fit anova asumptions
bestNormalize::bestNormalize(for_directed_hypothesis_testing$avg_fibe_tnfs)
for_directed_hypothesis_testing$avg_fibe_tnfs_boxcox <- (bestNormalize::boxcox(for_directed_hypothesis_testing$avg_fibe_tnfs))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- c("BoxCox, ANOVA, Tukey")
mod <- aov(avg_fibe_tnfs_boxcox ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$avg_fibe_tnfs_boxcox, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- c("0.308")
TukeyHSD(x = mod)
p_low_med_tukey <- c("0.3318526")
p_low_high_tukey <- c("0.9484312")
p_med_high_tukey <- c("0.5660176")
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(avg_fibe_tnfs),
                                                                       sd = sd(avg_fibe_tnfs),
                                                                       sem = se(avg_fibe_tnfs))
low_arg <- c("26.9 +/- 10.8  (1.36)")
medium_arg <- c("25.1 +/- 13.5  (1.20)")
high_arg <- c("26.4 +/- 12.3  (1.63)")

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
bestNormalize::bestNormalize(for_directed_hypothesis_testing$per_kcal_fiber_tnfs)
for_directed_hypothesis_testing$per_kcal_fiber_tnfs_log <- (bestNormalize::log_x(for_directed_hypothesis_testing$per_kcal_fiber_tnfs))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("Log, ANOVA, Tukey"))
mod <- aov(per_kcal_fiber_tnfs_log ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$per_kcal_fiber_tnfs_log, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- append(p_all,c("0.0453"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.0719109"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.9616280"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.1623831"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(per_kcal_fiber_tnfs),
                                                                       sd = sd(per_kcal_fiber_tnfs),
                                                                       sem = se(per_kcal_fiber_tnfs))
low_arg <- append(low_arg, c("12.4 +/- 4.61  (0.581)"))
medium_arg <- append(medium_arg, c("10.8 +/- 4.55  (0.405)"))
high_arg<- append(high_arg, c("12.1 +/- 5.21  (0.690)"))

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
bestNormalize::bestNormalize(for_directed_hypothesis_testing$total_fiber)
for_directed_hypothesis_testing$total_fiber_log <- (bestNormalize::log_x(for_directed_hypothesis_testing$total_fiber))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("Log, ANOVA, Tukey"))
mod <- aov(total_fiber_log ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$total_fiber_log, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- append(p_all,c("0.126"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.7162577"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.4998148"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.1047019"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(total_fiber),
                                                                       sd = sd(total_fiber),
                                                                       sem = se(total_fiber))
low_arg <- append(low_arg, c("26.0 +/- 11.4  (1.44)"))
medium_arg <- append(medium_arg, c("25.0 +/- 12.1  (1.08)"))
high_arg<- append(high_arg, c("28.4 +/- 12.0  (1.58)"))




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
bestNormalize::bestNormalize(for_directed_hypothesis_testing$total_fiber_per_kcal)
for_directed_hypothesis_testing$total_fiber_per_kcal_log <- (bestNormalize::log_x(for_directed_hypothesis_testing$total_fiber_per_kcal))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("Log, ANOVA, Tukey"))
mod <- aov(total_fiber_per_kcal_log ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$total_fiber_per_kcal_log, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- append(p_all,c("0.0298"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.0318458"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.7461438"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.2417922"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(total_fiber_per_kcal),
                                                                       sd = sd(total_fiber_per_kcal),
                                                                       sem = se(total_fiber_per_kcal))
low_arg <- append(low_arg, c("13.1 +/- 3.92  (0.494)"))
medium_arg <- append(medium_arg, c("11.7 +/- 3.82  (0.340)"))
high_arg<- append(high_arg, c("12.6 +/- 3.77  (0.500)"))


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
bestNormalize::bestNormalize(for_directed_hypothesis_testing$dt_fiber_sol)
for_directed_hypothesis_testing$dt_fiber_sol_log <- (bestNormalize::log_x(for_directed_hypothesis_testing$dt_fiber_sol))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("Log, ANOVA, Tukey"))
mod <- aov(dt_fiber_sol_log ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$dt_fiber_sol_log, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- append(p_all,c("0.137"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.6068548"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.6301818"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.1213973"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_fiber_sol),
                                                                       sd = sd(dt_fiber_sol),
                                                                       sem = se(dt_fiber_sol))
low_arg <- append(low_arg, c("7.20 +/- 3.23  (0.406)"))
medium_arg <- append(medium_arg, c("6.86 +/- 3.72  (0.331)"))
high_arg<- append(high_arg, c("7.73 +/- 3.32  (0.440)"))

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
bestNormalize::bestNormalize(for_directed_hypothesis_testing$dt_fiber_sol_per_kcal)
for_directed_hypothesis_testing$dt_fiber_sol_per_kcal_boxcox <- (bestNormalize::boxcox(for_directed_hypothesis_testing$dt_fiber_sol_per_kcal))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("BoxCox, ANOVA, Tukey"))
mod <- aov(dt_fiber_sol_per_kcal_boxcox ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$dt_fiber_sol_per_kcal_boxcox, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- append(p_all,c("0.0154"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.0143785"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.5906558"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.2460671"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_fiber_sol_per_kcal),
                                                                       sd = sd(dt_fiber_sol_per_kcal),
                                                                       sem = se(dt_fiber_sol_per_kcal))
low_arg <- append(low_arg, c("3.61 +/- 1.07  (0.135)"))
medium_arg <- append(medium_arg, c("3.18 +/- 1.05  (0.0936)"))
high_arg<- append(high_arg, c("3.44 +/- 1.13  (0.150)"))

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
bestNormalize::bestNormalize(for_directed_hypothesis_testing$dt_prot_animal)
for_directed_hypothesis_testing$dt_prot_animal_log<- (bestNormalize::log_x(for_directed_hypothesis_testing$dt_prot_animal))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("Log, ANOVA, Tukey"))
mod <- aov(dt_prot_animal_log ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$dt_prot_animal_log, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- append(p_all,c("0.0186"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.0260474"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.0440254"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.9694855"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_prot_animal),
                                                                       sd = sd(dt_prot_animal),
                                                                       sem = se(dt_prot_animal))
low_arg <- append(low_arg, c("46.9 +/- 24.7  (3.11)"))
medium_arg <- append(medium_arg, c("54.7 +/- 25.4  (2.26)"))
high_arg<- append(high_arg, c("55.6 +/- 27.4  (3.62)"))


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
summary(mod)
p_all <- append(p_all,c("0.116"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.1920954"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.1323840"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.8579190"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_sfat),
                                                                       sd = sd(dt_sfat),
                                                                       sem = se(dt_sfat))
low_arg <- append(low_arg, c("26.2 +/- 11.3  (1.42)"))
medium_arg <- append(medium_arg, c("29.1 +/- 12.6  (1.12)"))
high_arg<- append(high_arg, c("30.1 +/- 12.5  (1.65)"))

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
bestNormalize::bestNormalize(for_directed_hypothesis_testing$dt_kcal)
for_directed_hypothesis_testing$dt_kcal_boxcox<- (bestNormalize::boxcox(for_directed_hypothesis_testing$dt_kcal))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("BoxCox, ANOVA, Tukey"))
mod <- aov(dt_kcal_boxcox ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$dt_kcal_boxcox, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- append(p_all,c("0.124"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.4767918"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.1020351"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.4382766"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_kcal),
                                                                       sd = sd(dt_kcal),
                                                                       sem = se(dt_kcal))
low_arg <- append(low_arg, c("2002 +/- 739  (93.1)"))
medium_arg <- append(medium_arg, c("2146 +/- 816  (72.7)"))
high_arg<- append(high_arg, c("2289 +/- 825  (109)"))

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
bestNormalize::bestNormalize(for_directed_hypothesis_testing$pf_mps_total)
for_directed_hypothesis_testing$pf_mps_total_boxcox<- (bestNormalize::boxcox(for_directed_hypothesis_testing$pf_mps_total))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("BoxCox, ANOVA, Tukey"))
mod <- aov(pf_mps_total_boxcox ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$pf_mps_total_boxcox, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- append(p_all,c("0.0262"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.0264200"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.0887268"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.9952044"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(pf_mps_total),
                                                                       sd = sd(pf_mps_total),
                                                                       sem = se(pf_mps_total))
low_arg <- append(low_arg, c("3.30 +/- 2.18  (0.274)"))
medium_arg <- append(medium_arg, c("4.08 +/- 2.30  (0.205)"))
high_arg<- append(high_arg, c("4.12 +/- 2.50  (0.331)"))

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
bestNormalize::bestNormalize(for_directed_hypothesis_testing$pf_mps_total_per_kcal)
for_directed_hypothesis_testing$pf_mps_total_per_kcal_boxcox<- (bestNormalize::boxcox(for_directed_hypothesis_testing$pf_mps_total_per_kcal))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("BoxCox, ANOVA, Tukey"))
mod <- aov(pf_mps_total_per_kcal_boxcox ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$pf_mps_total_per_kcal_boxcox, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- append(p_all,c("0.3015"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.0409767"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.4193234"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.6326168"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(pf_mps_total_per_kcal),
                                                                       sd = sd(pf_mps_total_per_kcal),
                                                                       sem = se(pf_mps_total_per_kcal))
low_arg <- append(low_arg, c("1.63 +/- 0.831  (0.105)"))
medium_arg <- append(medium_arg, c("1.92 +/- 0.781  (0.0696)"))
high_arg<- append(high_arg, c("1.81 +/- 0.823  (0.109)"))

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
bestNormalize::bestNormalize(for_directed_hypothesis_testing$pf_meat)
for_directed_hypothesis_testing$pf_meat_log<- (bestNormalize::log_x(for_directed_hypothesis_testing$pf_meat))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("Log, ANOVA, Tukey"))
mod <- aov(pf_meat_log ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$pf_meat_log, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- append(p_all,c("0.0366"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.0300194"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.1770071"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.9042890"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(pf_meat),
                                                                       sd = sd(pf_meat),
                                                                       sem = se(pf_meat))
low_arg <- append(low_arg, c("1.12 +/- 0.914  (0.115)"))
medium_arg <- append(medium_arg, c("1.44 +/- 1.14  (0.101)"))
high_arg<- append(high_arg, c("1.47 +/- 1.18  (0.157)"))

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
bestNormalize::bestNormalize(for_directed_hypothesis_testing$pf_meat_per_kcal)
for_directed_hypothesis_testing$pf_meat_per_kcal_boxcox<- (bestNormalize::boxcox(for_directed_hypothesis_testing$pf_meat_per_kcal))$x.t
## many good options, Box-Cox transformation performed
transformation_test <- append(transformation_test,c("BoxCox, ANOVA, Tukey"))
mod <- aov(pf_meat_per_kcal_boxcox ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$pf_meat_per_kcal_boxcox, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- append(p_all,c("0.0629"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.0508300"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.5070591"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.5783957"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(pf_meat_per_kcal),
                                                                       sd = sd(pf_meat_per_kcal),
                                                                       sem = se(pf_meat_per_kcal))
low_arg <- append(low_arg, c("0.534 +/- 0.321  (0.0404)"))
medium_arg <- append(medium_arg, c("0.669 +/- 0.423  (0.0377)"))
high_arg<- append(high_arg, c("0.612 +/- 0.372  (0.0492)"))

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
## transform to fit anova asumptions
bestNormalize::bestNormalize(for_directed_hypothesis_testing$hei_asa24_totalscore)
## many good options, no transformation performed
transformation_test <- append(transformation_test,c("None, ANOVA, Tukey"))
mod <- aov(hei_asa24_totalscore ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$hei_asa24_totalscore, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- append(p_all,c("0.343"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.3534857"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.9286009"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.6325084"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(hei_asa24_totalscore),
                                                                       sd = sd(hei_asa24_totalscore),
                                                                       sem = se(hei_asa24_totalscore))
low_arg <- append(low_arg, c("64.5 +/- 13.3  (1.68)"))
medium_arg <- append(medium_arg, c("61.6 +/- 14.1  (1.26)"))
high_arg<- append(high_arg, c("63.6 +/- 12.0  (1.59)"))

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
## transform to fit anova asumptions
bestNormalize::bestNormalize(for_directed_hypothesis_testing$fecal_ph)
## many good options, no transformation performed
transformation_test <- append(transformation_test,c("None, ANOVA, Tukey"))
mod <- aov(fecal_ph ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$fecal_ph, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- append(p_all,c("0.149"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.4712824"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.1251837"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.5095283"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(fecal_ph),
                                                                       sd = sd(fecal_ph),
                                                                       sem = se(fecal_ph))
low_arg <- append(low_arg, c("7.08 +/- 0.543  (0.0683)"))
medium_arg <- append(medium_arg, c("6.98 +/- 0.614  (0.0547)"))
high_arg<- append(high_arg, c("6.87 +/- 0.523  (0.0693)"))

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
bestNormalize::bestNormalize(for_directed_hypothesis_testing$bmi_final)
for_directed_hypothesis_testing$bmi_final_boxcox<- (bestNormalize::boxcox(for_directed_hypothesis_testing$bmi_final))$x.t
## many good options, no transformation performed
transformation_test <- append(transformation_test,c("BoxCox, ANOVA, Tukey"))
mod <- aov(bmi_final_boxcox ~ cluster, data = for_directed_hypothesis_testing)
## re-test those assumptions
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$bmi_final_boxcox, group = for_directed_hypothesis_testing$cluster)
## test the anova and do post-hoc tests
summary(mod)
p_all <- append(p_all,c("0.3506"))
TukeyHSD(x = mod)
p_low_med_tukey <- append(p_low_med_tukey, c("0.8133816"))
p_low_high_tukey <- append(p_low_high_tukey, c("0.2344705"))
p_med_high_tukey <- append(p_med_high_tukey, c("0.4101401"))
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(bmi_final),
                                                                       sd = sd(bmi_final),
                                                                       sem = se(bmi_final))
low_arg <- append(low_arg, c("27.7 +/- 5.42  (0.683)"))
medium_arg <- append(medium_arg, c("27.2 +/- 5.01  (0.446)"))
high_arg<- append(high_arg, c("26.1 +/- 4.08  (0.540)"))

##################
## Complete the DF
##################

directed_hypothesis_table <- data.frame(Factor, data_source, 
           low_arg, medium_arg, high_arg, 
           transformation_test, p_low_med_tukey, p_low_high_tukey,p_med_high_tukey, p_all)

