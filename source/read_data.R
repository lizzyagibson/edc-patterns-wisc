# Robbie: cannot source this without MNData. I could comment and skip but thought I'd check with you before proceeding
source("./source/packages.R")

# this is a local data package
library(MNdata) 

# Results from BNMF

## scores
ewa <- readMat(here::here("./Data/mn2_EWA_sd1.mat"))[[1]] %>% as_tibble() %>% rename(P1 = V1, P2 = V2)

## unstandardized scores (only for viz)
ewa_raw <- readMat(here::here("./Data/mn2_EWA.mat"))[[1]] %>% 
  as_tibble() %>% 
  rename(P1 = V1, P2 = V2)

## loadings
eh <- readMat(here::here("./Data/mn2_EH.mat"))[[1]]
colnames(eh) <- c("MEHHP", "MECPP", "MEOHP", "MEHP", "MCPP", "MIBP", "MBP", "MBZP", "MEP", "BP_3", "B_PB", "M_PB", "P_PB", "TCS", "DCP_24", "DCP_25", "BPA")

# Mothers and Newborns data
whole = bind_cols(ppp_cov, ewa) %>% left_join(., mn_outcome) %>% rename(HOME_SCORE = HOMETOT)

# Filter
mn = whole %>% 
  filter(P1 != max(P1)) %>%
  filter(P1 != max(P1)) %>% # removed two outliers
  mutate_at(vars(c(HOME_SCORE, M_AGE, M_IQ)), scale) # scaled continuous values

# Flip so female is reference
mn_flip = mn %>% mutate(SEX = fct_rev(SEX))

# Subset so that n matches regression results
mn_subset = mn %>% dplyr::select(SID, WISC, P1, SEX, M_IQ, ALCOHOL, M_EDU, M_AGE,mat_hard,
                                 MARITAL_STATUS, HOME_SCORE) %>% drop_na() %>% 
  mutate(model = "Main Model")

# whole with all covariates (for sensitivity analysis)
whole_cc = whole %>% dplyr::select(SID, WISC, P1, SEX, M_IQ, ALCOHOL, M_EDU, M_AGE, mat_hard,
                                   MARITAL_STATUS, HOME_SCORE) %>% drop_na() %>% 
  mutate_at(vars(c(HOME_SCORE, M_AGE, M_IQ)), ~scale(.)[,1]) %>% 
  mutate(SEX = as_factor(str_c(SEX, "s")))
