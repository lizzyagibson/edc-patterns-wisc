library(tidyverse)
library(GGally)
library(R.matlab)
library(MNdata) # my data
library(janitor)
library(NMF)
library(reshape2)
library(viridis)
library(wesanderson)
library(ggridges)

source(here::here("../BN2MF/functions/compare_functions_2.R"))

theme_set(theme_bw(base_size = 20) + 
            theme(strip.background = element_rect(fill="white"),
                  legend.direction = "horizontal",
                  legend.title = element_blank()) # top, then right, bottom and left. 
)

# Read output
wa_dist = readMat(here::here("./Data/mn2_WA_dist.mat"))[[1]] 
dim(wa_dist)

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

#pdf("./Figures/scores_w_vci.pdf")
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
        axis.title = element_text(),
        legend.position = c(0.2, 0.87), # c(1,0) right bottom, c(1,1) right top.
        legend.background = element_rect(fill = "#ffffffaa", colour = NA)) +
  labs(x ="Study participants",
       y = "Individual scores") +
  scale_fill_manual(values = c("#0072B5", "#E18727"))
#dev.off()

# ggsci_db$"nejm"$"default" <- c(
#   "TallPoppy" = "#BC3C29", "DeepCerulean" = "#0072B5",
#   "Zest" = "#E18727", "Eucalyptus" = "#20854E",
#   "WildBlueYonder" = "#7876B1", "Gothic" = "#6F99AD",
#   "Salomie" = "#FFDC91", "FrenchRose" = "#EE4C97"
# )

# Plot distribution ####
dist_1 =  wa_dist[,1,] %>% 
  as_tibble() %>% 
  mutate(id = as_factor(1:nrow(.))) %>% 
  pivot_longer(V1:V1000) %>% 
  mutate(pattern = "Pattern 1")

dist_2 =  wa_dist[,2,] %>% 
  as_tibble() %>% 
  mutate(id = as_factor(1:nrow(.))) %>% 
  pivot_longer(V1:V1000) %>% 
  mutate(pattern = "Pattern 2")

vci_dist = bind_rows(dist_1, dist_2)

vci_dist %>% group_by(pattern) %>% 
  summarise(max = max(value),
            min = min(value))

pdf("./Figures/scores_w_vci_ridge.pdf")
vci_dist %>% 
  #filter(id %in% sample(1:343, 10)) %>% 
  filter(id %in% c(78, 172, 164, 90, 319, 307, 102, 22, sample(1:343, 2))) %>% # Highest and lowest
  arrange(pattern, id) %>% 
  mutate(id = fct_reorder(id, value)) %>% 
  arrange(pattern, id) %>% 
  mutate(id2 = as_factor(rep(rep(1:10, each=1000),2))) %>% 
  ggplot(aes(x = value, y = id2)) +
  geom_density_ridges(aes(fill=pattern), bandwidth=1, scale = 1, alpha = 0.75,
                      quantile_lines=TRUE,
                      quantile_fun=function(x,...) mean(x)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        # legend.position = "bottom",
        legend.direction = "vertical",
        legend.position = c(0.8, 0.15), # c(1,0) right bottom, c(1,1) right top.
        legend.background = element_rect(fill = "#ffffffaa", colour = NA)) +
  guides(fill = guide_legend(override.aes = list(linetype = 0))) +
  labs(y ="Study participants",
       x = "Individual scores") +
  scale_fill_manual(values = c("#0072B5", "#E18727"))
dev.off()
