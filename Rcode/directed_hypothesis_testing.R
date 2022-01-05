#######################
## directed hypotheses
#######################

## source or run merging_features.R to generate for_directed_hypothesis_testing object
# source(merging_features.R)
## load libraries
library(car)
library(patchwork)

## function of SEM 
## thanks John https://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean
se <- function(x) sqrt(var(x)/length(x))

## 1. test average fiber intake ASA24
leveneTest(for_directed_hypothesis_testing$avg_fibe_tnfs, group = for_directed_hypothesis_testing$cluster)
mod <- aov(avg_fibe_tnfs ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
## transform
for_directed_hypothesis_testing$avg_fibe_tnfs_boxcox <- (bestNormalize::boxcox(for_directed_hypothesis_testing$avg_fibe_tnfs))$x.t
mod <- aov(avg_fibe_tnfs_boxcox ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$avg_fibe_tnfs_boxcox, group = for_directed_hypothesis_testing$cluster)
summary(mod)
TukeyHSD(x = mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(avg_fibe_tnfs),
                                                                       sd = sd(avg_fibe_tnfs),
                                                                       sem = se(avg_fibe_tnfs))

## 2. test total fiber intake FFQ
leveneTest(for_directed_hypothesis_testing$total_fiber, group = for_directed_hypothesis_testing$cluster)
mod <- aov(total_fiber ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
## transform
bestNormalize(for_directed_hypothesis_testing$total_fiber)
for_directed_hypothesis_testing$total_fiber_log <- (bestNormalize::log_x(for_directed_hypothesis_testing$total_fiber))$x.t
mod <- aov(total_fiber_log ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$total_fiber_log, group = for_directed_hypothesis_testing$cluster)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(total_fiber),
                                                                       sd = sd(total_fiber),
                                                                       sem = se(total_fiber))

## 2b. test total fiber intake FFQ per kcal
for_directed_hypothesis_testing$total_fiber_per_kcal <- ((for_directed_hypothesis_testing$total_fiber / for_directed_hypothesis_testing$dt_kcal) * 1000)
leveneTest(for_directed_hypothesis_testing$total_fiber_per_kcal, group = for_directed_hypothesis_testing$cluster)
mod <- aov(total_fiber_per_kcal ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
## transform
bestNormalize(for_directed_hypothesis_testing$total_fiber_per_kcal)
for_directed_hypothesis_testing$total_fiber_per_kcal_log <- (bestNormalize::log_x(for_directed_hypothesis_testing$total_fiber_per_kcal))$x.t
mod <- aov(total_fiber_per_kcal_log ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$total_fiber_per_kcal_log, group = for_directed_hypothesis_testing$cluster)
bartlett.test(for_directed_hypothesis_testing$total_fiber_per_kcal_log ~ for_directed_hypothesis_testing$cluster)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(total_fiber_per_kcal),
                                                                       sd = sd(total_fiber_per_kcal),
                                                                       sem = se(total_fiber_per_kcal))

## 3. test for fiber (asa) per 1000 kcal
leveneTest(for_directed_hypothesis_testing$perKcal_fiber_tnfs, group = for_directed_hypothesis_testing$cluster)
mod <- aov(perKcal_fiber_tnfs ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
## transform
bestNormalize(for_directed_hypothesis_testing$perKcal_fiber_tnfs)
for_directed_hypothesis_testing$perKcal_fiber_tnfs_log <- (bestNormalize::log_x(for_directed_hypothesis_testing$perKcal_fiber_tnfs))$x.t
mod <- aov(perKcal_fiber_tnfs_log ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$perKcal_fiber_tnfs_log, group = for_directed_hypothesis_testing$cluster)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(perKcal_fiber_tnfs),
                                                                       sd = sd(perKcal_fiber_tnfs),
                                                                       sem = se(perKcal_fiber_tnfs))

## 4. test for solublefiber (ffq)
leveneTest(for_directed_hypothesis_testing$dt_fiber_sol, group = for_directed_hypothesis_testing$cluster)
mod <- aov(dt_fiber_sol ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
## transform
bestNormalize(for_directed_hypothesis_testing$dt_fiber_sol)
for_directed_hypothesis_testing$dt_fiber_sol_log <- (bestNormalize::log_x(for_directed_hypothesis_testing$dt_fiber_sol))$x.t
mod <- aov(dt_fiber_sol_log ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$dt_fiber_sol_log, group = for_directed_hypothesis_testing$cluster)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_fiber_sol),
                                                                       sd = sd(dt_fiber_sol),
                                                                       sem = se(dt_fiber_sol))

## 5. test for solublefiber per kcal (ffq)
for_directed_hypothesis_testing$dt_fiber_sol_per_kcal <- ((for_directed_hypothesis_testing$dt_fiber_sol / for_directed_hypothesis_testing$dt_kcal) * 1000)
leveneTest(for_directed_hypothesis_testing$dt_fiber_sol_per_kcal, group = for_directed_hypothesis_testing$cluster)
mod <- aov(dt_fiber_sol_per_kcal ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
## transform
bestNormalize(for_directed_hypothesis_testing$dt_fiber_sol_per_kcal)
for_directed_hypothesis_testing$dt_fiber_sol_per_kcal_log <- (bestNormalize::log_x(for_directed_hypothesis_testing$dt_fiber_sol_per_kcal))$x.t
mod <- aov(dt_fiber_sol_per_kcal_log ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$dt_fiber_sol_per_kcal_log, group = for_directed_hypothesis_testing$cluster)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_fiber_sol_per_kcal),
                                                                       sd = sd(dt_fiber_sol_per_kcal),
                                                                       sem = se(dt_fiber_sol_per_kcal))

## 6. Habitual protein intake from animals, g (ffq)
leveneTest(for_directed_hypothesis_testing$dt_prot_animal, group = for_directed_hypothesis_testing$cluster)
mod <- aov(dt_prot_animal ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
## transform
bestNormalize(for_directed_hypothesis_testing$dt_prot_animal)
for_directed_hypothesis_testing$dt_prot_animal_boxcox <- (bestNormalize::boxcox(for_directed_hypothesis_testing$dt_prot_animal))$x.t
mod <- aov(dt_prot_animal_boxcox ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$dt_prot_animal_boxcox, group = for_directed_hypothesis_testing$cluster)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_prot_animal),
                                                                       sd = sd(dt_prot_animal),
                                                                       sem = se(dt_prot_animal))

## 7. Saturated fatty acid intake, g (ffq)
leveneTest(for_directed_hypothesis_testing$dt_sfat, group = for_directed_hypothesis_testing$cluster)
mod <- aov(dt_sfat ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
## transform
bestNormalize(for_directed_hypothesis_testing$dt_sfat)
for_directed_hypothesis_testing$dt_sfat_boxcox <- (bestNormalize::boxcox(for_directed_hypothesis_testing$dt_sfat))$x.t
mod <- aov(dt_sfat_boxcox ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$dt_sfat_boxcox, group = for_directed_hypothesis_testing$cluster)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_sfat),
                                                                       sd = sd(dt_sfat),
                                                                       sem = se(dt_sfat))

## 8. Total energy intake, kcal (ffq)
leveneTest(for_directed_hypothesis_testing$dt_kcal, group = for_directed_hypothesis_testing$cluster)
mod <- aov(kcal_total ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids) # close enough
## transform
bestNormalize(for_directed_hypothesis_testing$dt_kcal)
for_directed_hypothesis_testing$dt_kcal_boxcox <- (bestNormalize::boxcox(for_directed_hypothesis_testing$dt_kcal))$x.t
mod <- aov(dt_kcal_boxcox ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$dt_kcal_boxcox, group = for_directed_hypothesis_testing$cluster)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(dt_kcal),
                                                                       sd = sd(dt_kcal),
                                                                       sem = se(dt_kcal))

## 8. Beef/pork consumption excluding cured meat, g (ffq)
leveneTest(for_directed_hypothesis_testing$pf_meat, group = for_directed_hypothesis_testing$cluster)
mod <- aov(pf_meat ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
## transform
bestNormalize(for_directed_hypothesis_testing$pf_meat)
for_directed_hypothesis_testing$pf_meat_boxcox <- (bestNormalize::boxcox(for_directed_hypothesis_testing$pf_meat))$x.t
mod <- aov(pf_meat_boxcox ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$pf_meat_boxcox, group = for_directed_hypothesis_testing$cluster)
summary(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(pf_meat),
                                                                       sd = sd(pf_meat),
                                                                       sem = se(pf_meat))

## 8. Beef/pork consumption including cured meat, oz eq (ffq)
leveneTest(for_directed_hypothesis_testing$pf_mps_total, group = for_directed_hypothesis_testing$cluster)
mod <- aov(pf_mps_total ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
## transform
bestNormalize(for_directed_hypothesis_testing$pf_mps_total)
for_directed_hypothesis_testing$pf_mps_total_boxcox <- (bestNormalize::boxcox(for_directed_hypothesis_testing$pf_mps_total))$x.t
mod <- aov(pf_mps_total_boxcox ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$pf_mps_total_boxcox, group = for_directed_hypothesis_testing$cluster)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(pf_mps_total),
                                                                       sd = sd(pf_mps_total),
                                                                       sem = se(pf_mps_total))

## 9. Beef/pork consumption including cured meat, oz eq per 1000kcal (ffq)
for_directed_hypothesis_testing$pf_mps_total_per_kcal <- (for_directed_hypothesis_testing$pf_mps_total / for_directed_hypothesis_testing$dt_kcal) * 1000
leveneTest(for_directed_hypothesis_testing$pf_mps_total_per_kcal, group = for_directed_hypothesis_testing$cluster)
mod <- aov(pf_mps_total_per_kcal ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
## transform
bestNormalize(for_directed_hypothesis_testing$pf_mps_total_per_kcal)
for_directed_hypothesis_testing$pf_mps_total_per_kcal_boxcox <- (bestNormalize::boxcox(for_directed_hypothesis_testing$pf_mps_total_per_kcal))$x.t
mod <- aov(pf_mps_total_per_kcal_boxcox ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$pf_mps_total_per_kcal_boxcox, group = for_directed_hypothesis_testing$cluster)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(pf_mps_total_per_kcal),
                                                                       sd = sd(pf_mps_total_per_kcal),
                                                                       sem = se(pf_mps_total_per_kcal))

## 10. Beef/pork consumption including cured meat, oz eq (ffq)
leveneTest(for_directed_hypothesis_testing$pf_meat, group = for_directed_hypothesis_testing$cluster)
mod <- aov(pf_meat ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
## transform
bestNormalize(for_directed_hypothesis_testing$pf_meat)
for_directed_hypothesis_testing$pf_meat_boxcox <- (bestNormalize::boxcox(for_directed_hypothesis_testing$pf_meat))$x.t
mod <- aov(pf_meat_boxcox ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$pf_meat_boxcox, group = for_directed_hypothesis_testing$cluster)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(pf_meat),
                                                                       sd = sd(pf_meat),
                                                                       sem = se(pf_meat))

## 11. Beef/pork consumption excluding cured meat, oz eq per 1000kcal (ffq)
for_directed_hypothesis_testing$pf_meat_per_kcal <- (for_directed_hypothesis_testing$pf_meat / for_directed_hypothesis_testing$dt_kcal) * 1000
leveneTest(for_directed_hypothesis_testing$pf_meat_per_kcal, group = for_directed_hypothesis_testing$cluster)
mod <- aov(pf_meat_per_kcal ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
## transform
bestNormalize(for_directed_hypothesis_testing$pf_meat_per_kcal)
for_directed_hypothesis_testing$pf_meat_per_kcal_boxcox <- (bestNormalize::boxcox(for_directed_hypothesis_testing$pf_meat_per_kcal))$x.t
mod <- aov(pf_meat_per_kcal_boxcox ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
leveneTest(for_directed_hypothesis_testing$pf_meat_per_kcal_boxcox, group = for_directed_hypothesis_testing$cluster)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(pf_meat_per_kcal),
                                                                       sd = sd(pf_meat_per_kcal),
                                                                       sem = se(pf_meat_per_kcal))

## 12. HEI total score (asa)
leveneTest(for_directed_hypothesis_testing$hei_asa24_totalscore, group = for_directed_hypothesis_testing$cluster)
mod <- aov(hei_asa24_totalscore ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(hei_asa24_totalscore),
                                                                       sd = sd(hei_asa24_totalscore),
                                                                       sem = se(hei_asa24_totalscore))

## 13. Stool pH
leveneTest(for_directed_hypothesis_testing$fecal_ph, group = for_directed_hypothesis_testing$cluster)
mod <- aov(fecal_ph ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(fecal_ph),
                                                                       sd = sd(fecal_ph),
                                                                       sem = se(fecal_ph))

## 13. BMI
leveneTest(for_directed_hypothesis_testing$bmi_final, group = for_directed_hypothesis_testing$cluster)
mod <- aov(bmi_final ~ cluster, data = for_directed_hypothesis_testing)
resids <- residuals(mod)
shapiro.test(resids)
for_directed_hypothesis_testing$bmi_final_log <- (bestNormalize::log_x(for_directed_hypothesis_testing$bmi_final))$x.t
mod <- aov(bmi_final_log ~ cluster, data = for_directed_hypothesis_testing)
summary(mod)
TukeyHSD(mod)
## describe
for_directed_hypothesis_testing %>% group_by(., cluster) %>% summarise(., average = mean(bmi_final),
                                                                       sd = sd(bmi_final),
                                                                       sem = se(bmi_final))
