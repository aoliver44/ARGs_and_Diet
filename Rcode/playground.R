danielle_boxplots <- for_directed_hypothesis_testing %>% select(., subject_id, cluster, t_flavones, t_flavonoids, sex) %>% drop_na() #%>%
  dplyr::mutate(., `dheas:ur_cort` = (dheas_bd1 / ur_cort_ug_gcreat)) %>% dplyr::mutate(., `dheas:saliva` = (dheas_bd1 / as.numeric(`Fasted Salivary Cortisol nmol/L`)))
danielle_boxplots <- melt(id.vars = c("subject_id", "cluster", "sex"), data = danielle_boxplots)

danielle_boxplots <- danielle_boxplots %>% mutate(., value = ifelse(value == "soft", 0,
                                                                    ifelse(value == "normal", 1, 
                                                                           ifelse(value == "hard", 2, value))))

danielle_boxplots$cluster <- factor(x = danielle_boxplots$cluster, levels = c("low", "medium", "high"), ordered = T)
danielle_boxplots %>%
  ggplot() + aes(x = as.factor(cluster), y = as.numeric(value)) + 
  geom_boxplot(aes(fill = as.factor(cluster)),outlier.alpha = 0) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.2) +
  facet_wrap( ~ variable, scales = "free_y") +
  scale_fill_jama() + theme_bw(base_size = 20) +
  labs(x = "Cluster", y = "Value") + theme(legend.position = "none") +
  stat_compare_means(method = "kruskal") 

tmp <- for_directed_hypothesis_testing %>% 
  select(., cluster, sex) %>%
  group_by(., cluster, sex) %>%
  tally
tmp <- as.data.frame(tmp)
tmp$cluster <- factor(x = tmp$cluster, levels = c("low", "medium", "high"), ordered = T)

ggplot(tmp) + aes(x = as.factor(sex), weight = as.numeric(n)) +
  geom_bar() + facet_wrap( ~ as.factor(cluster)) + labs(y = "Number of individuals", title = "Stool consistancy faceted by AMR Cluster")
 

for_directed_hypothesis_testing %>%
  group_by(., cluster) %>%
  summarise(., mean_age = mean(per_kcal_fiber_tnfs), sd_age = sd(per_kcal_fiber_tnfs))

mpa_taxa_group_otu = mpa_taxa_group_otu[,colSums(mpa_taxa_group_otu) > 500]
tmp <- mpa_taxa_group_otu %>%
  column_to_rownames(var = "subject_id") %>% 
  corrr::correlate(method = "pearson") %>% 
  corrr::focus(plasma_lbp_bd1) 
  # cut off of 0.318 determined below...keeps all padjusted sig stuff in
  #filter(abs(plasma_lbd_bd1) > 0.08) %>%
  #mutate(term = factor(term, levels = term[order(plasma_lbd_bd1)])) %>%
ggplot(data = tmp, aes(x = term, y = tmp$plasma_lbp_bd1)) +
  geom_bar(stat = "identity") +
  ylab("Correlation with Enterobacteriaceae\n (Pearson)") +
  xlab("Genus") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()


pre_corr_clean %>%
  pivot_longer(GNP.deflator:Armed.Forces, names_to="x_var", values_to="x_val") %>% 
  pivot_longer(Population:Employed, names_to="y_var", values_to="y_val") %>% 
  nest(data=c(x_val, y_val)) %>%
  mutate(cor_test = map(data, ~cor.test(.x$x_val, .x$y_val)),
         tidied = map(cor_test, tidy)) %>% 
  unnest(tidied)


mtx = cor(pre_corr_clean, method = "spearman")
drop = caret::findCorrelation(mtx, cutoff = .85, exact = T)
drop = names(pre_corr_clean)[drop]

## correlating food diversity and AMR alpha
for_elizabeth <- merge(alpha_diversity_data, fiber_alpha, by = "subject_id")
for_elizabeth$cluster <- factor(for_elizabeth$cluster, levels = c("low", "medium", "high"), ordered = T)
cor(for_elizabeth$richness, for_elizabeth$pd_whole_tree_fiber, method = "pearson")

for_elizabeth1 <- merge(for_directed_hypothesis_testing, alpha_diversity_data, by = "subject_id")
for_elizabeth1 %>%
  dplyr::select(where(is.numeric)) %>%
  column_to_rownames(var = "subject_id") %>% 
  corrr::correlate(method = "spearman") %>%
  corrr::focus(shannon_div) %>% View()
# cut off of 0.318 determined below...keeps all padjusted sig stuff in
filter(abs(shannon_div) > 0.15) %>%
mutate(term = factor(term, levels = term[order(shannon_div)])) %>%
ggplot() + aes(x = term, y = shannon_div) +
  geom_bar(stat = "identity") +
  ylab("Correlation with shannon_div\n (Spearman)") +
  xlab("Genus") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()



## kmeans clustering for age, sex, bmi

#' Plots a chart showing the sum of squares within a group for each execution of the kmeans algorithm. 
#' In each execution the number of the initial groups increases by one up to the maximum number of centers passed as argument.
#'
#' @param data The dataframe to perform the kmeans 
#' @param nc The maximum number of initial centers
#'
wssplot <- function(data, nc=15, seed=123){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of groups",
       ylab="Sum of squares within a group")}

kmeans_input <- for_directed_hypothesis_testing %>% select(., age, sex, bmi_final)

wssplot(kmeans_input, nc = 10)

set.seed(123)
clustering <- kmeans(kmeans_input, centers = 6, nstart = 20)
clustering


## output file for lefse
for_lefse <- merge(abx_cluster, mpa_taxa_group_otu, by.x = "subject_id", by.y = "variable")


## correlation with pH
library(ggpubr)

fiber_ph <- for_directed_hypothesis_testing %>% select(., fecal_ph, dt_fibe, cluster)
fiber_ph$cluster <- factor(fiber_ph$cluster, levels = c("low", "medium", "high"), ordered = T)
ggscatter(data = fiber_ph, x = "dt_fibe", y = "fecal_ph", color = "cluster", 
          add = "reg.line", conf.int = TRUE, palette = pal_jama()(3), rug = T) 


