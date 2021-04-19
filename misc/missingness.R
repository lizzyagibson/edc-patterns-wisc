source("./packages.R")
source("./read_data.R")

model_dat = whole %>% 
  filter(!is.na(WISC)) %>% 
  dplyr::select(-c(SID, GEST, WISC_VC, WISC_PR, WISC_PS, WISC_WM, P2))

model_dat %>% drop_na()
model_dat %>% names()

model_dat %>% ff_glimpse()

mice::md.pattern(model_dat, rotate.names = T)
# 289 complete samples
# 14 missing only home score
# 6 missing only maternal iq
# 2 missing both

aggr(model_dat, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, 
     labels=names(model_dat), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

marginplot(model_dat[c(6,9)])

explanatory = c("ETH","M_EDU","MARITAL_STATUS","SMOKER_IN_HOME","SEX","HOME_SCORE",
                "ALCOHOL","M_AGE","mat_hard","WISC","P1")

model_m = model_dat %>% 
  mutate(home_miss = as_factor(ifelse(is.na(HOME_SCORE), 1, 0)),
         iq_miss = as_factor(ifelse(is.na(M_IQ), 1, 0)))

model_m %>% missing_compare(., dependent = "M_IQ", explanatory)

explanatory = c("ETH","M_EDU","MARITAL_STATUS","SMOKER_IN_HOME","SEX","M_IQ",
                "ALCOHOL","M_AGE","mat_hard","WISC","P1")

model_m %>% missing_compare(., dependent = "HOME_SCORE", explanatory) 



