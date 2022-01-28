##################
## Supp Fig. Lefse
##################

library(tidyverse)
library(ggsci)
setwd("/home/datasets/new_datasets/")
## supp figure 2

lefse <- read.table(file = "lefse_family_raw.csv", sep = ",", header = T)[1:4]
lefse$Cluster <- factor(lefse$Cluster, levels = c("low", "medium", "high"), ordered = T)
lefse_plot <- lefse %>% ggplot() + aes(x = reorder(Family,-LDA), weight = LDA) +
  geom_bar(aes(fill = Cluster)) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 90),
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  facet_grid(~ Cluster, scales = "free_x") +
  scale_fill_jama()
