library(tidyverse)
library(R.matlab)

ewa  <- readMat(here::here("./Data/mn2_EWA_un.mat"))[[1]] %>% as_tibble()
var_wa <- readMat(here::here("./Data/mn2_WA_var_un.mat"))[[1]] %>% as_tibble() %>% rename(var1 = V1, var2 = V2)

Walpha <- readMat(here::here("./Data/mn2_Walpha.mat"))[[1]] %>% as_tibble() %>% rename(Walpha1 = V1, Walpha2 = V2)
Wbeta <- readMat(here::here("./Data/mn2_Wbeta.mat"))[[1]] %>% as_tibble() %>% rename(Wbeta1 = V1, Wbeta2 = V2)

Aalpha <- readMat(here::here("./Data/mn2_Aalpha.mat"))[[1]] %>% as_tibble() %>% rename(Aalpha1 = V1, Aalpha2 = V2)
Abeta  <- readMat(here::here("./Data/mn2_Abeta.mat"))[[1]] %>% as_tibble() %>% rename(Abeta1 = V1, Abeta2 = V2)

ewa_final  <- readMat(here::here("./Data/mn2_EWA_sd1.mat"))[[1]] %>% as_tibble() %>% rename(ewa_final1 = V1, ewa_final2 = V2)
var_wa_final <- readMat(here::here("./Data/mn2_WA_var_sd1.mat"))[[1]] %>% as_tibble() %>% rename(var_final1 = V1, var_final2 = V2)

Aalpha
Abeta

# scaled normal
final_norm_dist = bind_cols(ewa_final, var_wa_final) %>% 
  mutate(FinalNormal_dist1 = map2(ewa_final1, var_final1, function(x,y) rnorm(1000, mean = x, sd = sqrt(y))),
         FinalNormal_dist2 = map2(ewa_final2, var_final2, function(x,y) rnorm(1000, mean = x, sd = sqrt(y)))) %>% 
  dplyr::select(-c(ewa_final1, ewa_final2, var_final1, var_final2))

# raw normal
norm_dist = bind_cols(ewa, var_wa) %>% 
  mutate(Normal_dist1 = map2(V1, var1, function(x,y) rnorm(1000, mean = x, sd = sqrt(y))),
         Normal_dist2 = map2(V2, var2, function(x,y) rnorm(1000, mean = x, sd = sqrt(y)))) %>% 
  dplyr::select(-c(V1,     V2,    var1,    var2))

# raw gamma
w_dist = bind_cols(Walpha, Wbeta) %>% 
  mutate(Wdist1 = map2(Walpha1, Wbeta1, function(x,y) rgamma(1000, shape = x, rate = y)),
         Wdist2 = map2(Walpha2, Wbeta2, function(x,y) rgamma(1000, shape = x, rate = y)))

a_dist = bind_cols(Aalpha, Abeta) %>% 
  mutate(Adist1 = map2(Aalpha1, Abeta1, function(x,y) rgamma(1000, shape = x, rate = y)),
         Adist2 = map2(Aalpha2, Abeta2, function(x,y) rgamma(1000, shape = x, rate = y)))  

adist1 = a_dist$Adist1[[1]]
adist2 = a_dist$Adist2[[1]]

wa_dist = w_dist %>% 
  mutate(Gamma_dist1 = map(Wdist1, function(x) x*adist1),
         Gamma_dist2 = map(Wdist2, function(x) x*adist2)) %>% 
  dplyr::select(-grep("(alpha|beta|Wdist)", names(.)))

# scaled gamma
w_dist_s = bind_cols(Walpha, Wbeta) %>% 
  mutate(Wdist1_s = map2(Walpha1, Wbeta1, function(x,y) rgamma(1000, shape = x, rate = (y/17.3043)*4.4074)),
         Wdist2_s = map2(Walpha2, Wbeta2, function(x,y) rgamma(1000, shape = x, rate = (y/17.2979)*3.3636)))

a_dist_s = bind_cols(Aalpha, Abeta) %>% 
  mutate(Adist1_s = map2(Aalpha1, Abeta1, function(x,y) rgamma(1000, shape = x, rate = (y/17.3043)*4.4074)),
         Adist2_s = map2(Aalpha2, Abeta2, function(x,y) rgamma(1000, shape = x, rate = (y/17.2979)*3.3636)))  

adist1_s = a_dist_s$Adist1_s[[1]]
adist2_s = a_dist_s$Adist2_s[[1]]

wa_dist_s = w_dist %>% 
  mutate(ScaledGamma_dist1 = map(Wdist1, function(x) x*adist1_s), # don't need scaled a // double scaling
         ScaledGamma_dist2 = map(Wdist2, function(x) x*adist2_s)) %>% 
  dplyr::select(-grep("(alpha|beta|Wdist)", names(.)))
  
both_dists = bind_cols(wa_dist, norm_dist) %>% 
  mutate(id = 1:nrow(.)) %>% 
  unnest(Gamma_dist1:Normal_dist2) %>% 
  pivot_longer(Gamma_dist1:Normal_dist2,
               names_to = c("dist", "pattern"),
               names_sep = "_")
both_dists

both_dists %>% 
  mutate(dist = ifelse(dist == "Gamma", "Gamma*Gamma", "Normal")) %>% 
  filter(id %in% 1:12 & pattern == "dist1") %>% 
  ggplot(aes(x = value)) +
  geom_histogram(aes(fill = dist), bins = 20,
                 position = "identity", alpha = 0.55) +
  facet_wrap(.~id, scales = "free") +
  theme_bw() + 
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  labs(x = "", y = "Count")

all_dists = bind_cols(wa_dist, norm_dist, final_norm_dist, wa_dist_s) %>% 
  mutate(id = 1:nrow(.)) %>% 
  unnest(Gamma_dist1:ScaledGamma_dist2) %>% 
  pivot_longer(Gamma_dist1:ScaledGamma_dist2,
               names_to = c("dist", "pattern"),
               names_sep = "_")
all_dists

all_dists %>% 
  filter(id %in% 1:4 & pattern == "dist1") %>% 
  #filter(dist != "ScaledGamma") %>% 
  ggplot(aes(x = value)) +
  geom_histogram(aes(fill = dist), bins = 20,
                 position = "identity", alpha = 0.55) +
  facet_wrap(.~id, scales = "free") +
  theme_bw() + 
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  labs(x = "", y = "Count")

# This looks good!