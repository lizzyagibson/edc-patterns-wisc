library(tidyverse)
library(GGally)
library(R.matlab)
library(MNdata) # my data
library(janitor)
library(NMF)
library(reshape2)
library(viridis)
library(wesanderson)

source(here::here("../BN2MF/functions/compare_functions.R"))
source(here::here("../BN2MF/functions/fig_set.R"))

# Read output

lower = readMat(here::here("./Data/mn2_WA_lower.mat"))[[1]] %>% as_tibble()
ewa = readMat(here::here("./Data/mn2_EWA.mat"))[[1]] %>% as_tibble()
upper = readMat(here::here("./Data/mn2_WA_upper.mat"))[[1]] %>% as_tibble()

lower = lower %>% 
  mutate(id = 1:nrow(.)) %>% 
  pivot_longer(V1:V2,
               names_to = "pattern",
               values_to = "lower")

ewa = ewa %>% 
  mutate(id = 1:nrow(.)) %>% 
  pivot_longer(V1:V2,
               names_to = "pattern",
               values_to = "ewa")

upper = upper %>% 
  mutate(id = 1:nrow(.)) %>% 
  pivot_longer(V1:V2,
               names_to = "pattern",
               values_to = "upper")

vci = full_join(lower, ewa) %>% full_join(., upper) %>% 
  mutate(pattern = ifelse(pattern == "V1", "Pattern 1", "Pattern 2"))

vci %>% group_by(pattern) %>% 
  summarise(max = max(ewa),
            min = min(ewa))
vci %>% 
  arrange((ewa))

pdf("./Figures/scores_w_vci.pdf")
vci %>% 
  mutate(id = as.factor(id)) %>% 
  #filter(id %in% sample(1:343, 10)) %>% 
  filter(id %in% c(78, 172, 164, 90, 319, 307, 102, 22, sample(1:343, 2))) %>% # Highest and lowest
  arrange(pattern, id) %>% 
  mutate(id = fct_reorder(id, ewa)) %>% 
  arrange(pattern, id) %>% 
  mutate(id2 = as.factor(rep(1:10, 2))) %>% 
  ggplot(aes(x = id2, y = ewa, group = pattern, fill = pattern)) +
  #geom_point(size = 4, position = position_dodge(.5)) +
  geom_col(width = 0.7, position = position_dodge(0.8), color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(.8), color = "black", width = 0.3) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        legend.direction = "vertical",
        legend.position = c(0.2, 0.87), # c(1,0) right bottom, c(1,1) right top.
        legend.background = element_rect(fill = "#ffffffaa", colour = NA)) +
  labs(x ="Study participants",
       y = "Individual scores") +
  scale_fill_manual(values = c("#0072B5", "#E18727"))
dev.off()

# ggsci_db$"nejm"$"default" <- c(
#   "TallPoppy" = "#BC3C29", "DeepCerulean" = "#0072B5",
#   "Zest" = "#E18727", "Eucalyptus" = "#20854E",
#   "WildBlueYonder" = "#7876B1", "Gothic" = "#6F99AD",
#   "Salomie" = "#FFDC91", "FrenchRose" = "#EE4C97"
# )

