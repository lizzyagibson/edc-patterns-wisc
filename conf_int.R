# Functions ####
tidy_ci = function(fit) {
  tidy(fit) %>% bind_cols(., as_tibble(confint(fit)))
}

add_ci4interaction <- function(fit, term1, term2) {
  # CONFIDENCE INTERVAL for pattern in females
  # Compute association and its uncertainty 
  
  # here we create tables of coefficients and covariance
  coef.mat <- summary(fit)$coefficients
  var.mat  <- vcov(fit)
  
  # the total term for the association is the 
  # sum of P2 in the reference sex plus the term for P2:female
  beta.Pf <- coef.mat[term1,1] + coef.mat[term2,1]
  
  # Compute variance in order to compute standard error
  # We must compute the variance for the total term 
  # Var(Beta1 + Beta2) = Var(Beta1) + Var(Beta2) + 2*CoVar(Beta1, Beta2) 
  var.Pf <- var.mat[term1, term1] + 
    var.mat[term2, term2] +
    2*var.mat[term1, term2]
  
  # this is st error NOT std
  ste.Pf  <-  sqrt(abs(var.Pf))
  
  # compute confidence intervals 
  lci.Pf <- beta.Pf - 1.96*ste.Pf
  uci.Pf <- beta.Pf + 1.96*ste.Pf
  
  # 2 calculate the test statistic: t = Est/SE
  # 3 calculate the P value2: P = exp(−0.717×z − 0.416×z2).
  test_stat = abs(beta.Pf/ste.Pf)
  pvalue = exp(-0.717*test_stat - 0.416*test_stat^2)
  pvalue
  
  tidy_ci(fit) %>% 
    bind_rows(., tibble(term = paste0(term1, " in females"), estimate = beta.Pf, std.error = ste.Pf,
                        statistic = test_stat, p.value = pvalue, `2.5 %` = lci.Pf, `97.5 %` = uci.Pf)) %>% 
    mutate(term = ifelse(term == term1, paste(term1, "in males"), term))
}

add_ci4int_bayes = function(fit) {
  int_sum = summary(fit, c("beta_int", "beta_sex", "beta_p"))$summary
  
  params = names(p1_non_const)[grep("(alpha|beta)", names(p1_non_const))]
  model_sum = summary(fit, params)$summary
  
  ext <- extract(fit)
  
  post_coef = cbind(ext$alpha, ext$beta_c, 
                    ext$beta_int, ext$beta_sex, ext$beta_p) %>% as_tibble()
  
  colnames(post_coef) = params
  int_coefs = post_coef %>% dplyr::select(beta_p, beta_int)
  
  coef1 = int_coefs$beta_p
  coef2 = int_coefs$beta_int
  
  covmat = cov(int_coefs)
  beta_f = coef1 + coef2
  
  var_f = covmat[1, 1] + 
    covmat[2, 2] +
    2*covmat[1, 2]
  
  # this is std dev not std error
  sd_f  <- sqrt(abs(var_f))
  
  # take min sample size as conservative est.
  sampsize = min(int_sum[1,9], int_sum[3,9])
  
  # this is std error
  se_f  <- sd_f/sqrt(sampsize)
  
  lci.f <- mean(beta_f) - 1.96*se_f
  uci.f <- mean(beta_f) + 1.96*se_f
  
  # 2 calculate the test statistic: t = Est/SE
  # 3 calculate the P value2: P = exp(−0.717×z − 0.416×z2).
  # test_stat = abs(mean(beta_f)/se_f)
  # pvalue = exp(-0.717*test_stat - 0.416*test_stat^2)
  # NO BAYESIAN PVALUE
  
  q = quantile(beta_f, probs = c(0.025, .25, .50, .75, .975)) %>% as.data.frame() %>% 
    rownames_to_column() %>% 
    pivot_wider(names_from = "rowname",
                values_from = ".")
  
  new_var = tibble(beta = "p in females", mean = round(mean(beta_f), 4), 
                   se_mean = round(se_f,4), sd = round(sd_f,4)) %>% 
    bind_cols(., q)
  
  # M_AGE + M_EDU + MARITAL_STATUS + HOME_SCORE + M_IQ + ALCOHOL + mat_hard, 
  as.data.frame(model_sum)[,-c(9:10)] %>% 
    rownames_to_column(var="beta") %>% 
    mutate(beta = c("alpha", "age", "edu", "marital", "home", "iq", "alcohol", "mat_hard",
                    "sex*p", "sex", "p in males")) %>% 
    bind_rows(.,new_var) %>% 
    mutate_at(vars(2:9), round, 2)
}

# Data ####
ewa <- readMat(here::here("./Data/mn2_EWA_sd1.mat"))[[1]] %>% as_tibble() %>% rename(P1 = V1, P2 = V2)
whole = bind_cols(ppp_cov, ewa) %>% left_join(., mn_outcome) %>% rename(HOME_SCORE = HOMETOT)

mn = whole %>% filter(P1 < mean(whole$P1) + 5*sd(whole$P1)) %>% 
  mutate_at(vars(c(HOME_SCORE, M_AGE, M_IQ)), scale)

# so female is reference
mn_flip = mn %>% mutate(SEX = fct_rev(SEX))

# Subset so that n matches regression results
mn_subset = mn %>% dplyr::select(SID, WISC, P = P1, SEX, M_IQ, ALCOHOL, M_EDU, M_AGE,mat_hard,
                                 MARITAL_STATUS, HOME_SCORE) %>% drop_na() %>% 
  mutate(model = "Main Model")

# OLS ####
load(file = "reg_pat1_fit.rda")
reg_sum = add_ci4interaction(fit_p1a, "P1", "P1:SEXFemale")
reg_sum

fit_flip <- lm(WISC ~ P1 + SEX*P1 + SEX + M_IQ + ALCOHOL + M_EDU  + M_AGE +
                 MARITAL_STATUS + HOME_SCORE + mat_hard
               , data = mn_flip)
summary(fit_flip)

# Get slope and se from regression model
f_slope = reg_sum[12,2][[1]]
m_slope = reg_sum[2,2][[1]]

fint = 102.1678 + (0.8503) - 
  (5.4022*0.2517483 ) - # alcohol
  (1.3726*0.3531469) - # edu 
  (0.1713*0.3251748) - # married
  (2.6254*0.4020979 ) # material hardship

mint = 102.1678 - 
  (5.4022*0.2517483 ) - # alcohol
  (1.3726*0.3531469) - # edu 
  (0.1713*0.3251748) - # married
  (2.6254*0.4020979 ) # material hardship

# Get predicted values
# this is conditional on everything else!
pred = as.data.frame(predict(fit_p1a, se.fit = TRUE, type = c("terms")))
head(pred)[,c(1,2,10,11,12,20)]

# Get predicted values for females
pred_flip = as.data.frame(predict(fit_flip, se.fit = TRUE, type = c("terms")))
head(pred_flip)[,c(1,2,10,11,12,20)] # se.fit.P1 = 0.5440310

# Combine predicted outcomes with data
main = 
  mn_subset %>%
  bind_cols(., pred) %>% 
  bind_cols(., flip.se = (pred_flip$se.fit.P1 + pred_flip$se.fit.P1.SEX)) %>% # just get se for female
  dplyr::select(SID, WISC, grep("(SEX|P)", names(.)), flip.se) %>% 
   # -c(M_IQ:model, df, residual.scale, fit.M_IQ:fit.mat_hard, se.fit.M_IQ:se.fit.mat_hard)) %>% 
  mutate(# predIQm = fit.P1 + mean(WISC) + fit.SEX + fit.P1.SEX,
         # predIQf = fit.P1 + mean(WISC) + fit.SEX + fit.P1.SEX, 
         # fit.sex changes intercept, interaction changes slope
         predIQ = fit.P1 + mean(WISC) + fit.SEX + fit.P1.SEX, # ifelse(SEX == "Female", predIQf, predIQm),
         se.male = se.fit.P1 + se.fit.P1.SEX,
         se.female = flip.se,
         se = ifelse(SEX == "Female", se.female, se.male),
         lci = predIQ - 1.96*se,
         uci = predIQ + 1.96*se)
main

ols_plot = main %>% 
  ggplot(aes(x = P)) +
  geom_point(aes(y = WISC), color = "lightgray", size = 0.5) + 
  geom_abline(intercept = mean(main$WISC), slope = 0, 
              color = "darkgray", linetype = "dashed") + 
  # geom_abline(slope = f_slope, intercept = fint,color="darkgray") +
  # geom_abline(slope = m_slope, intercept = mint,color="darkgray") +
  geom_line(aes(y = predIQ)) +
  #geom_line(aes(y = lci), linetype = "dotted") +
  #geom_line(aes(y = uci), linetype = "dotted") +
  geom_ribbon(aes(ymin = lci, ymax = uci),
              alpha = 0.25, linetype = "dotted", size = 0.25) +
  facet_grid(.~SEX, scales = "free_x") +
  ylim(45, 135) +
  theme_minimal(base_size = 15) + 
  theme(legend.title = element_blank(),
        legend.position = "none") +
  scale_color_lancet() + scale_fill_lancet() +
  labs(y = "WISC full scale IQ", x = "PHT pattern concentration") 

# BAYESIAN ####
load("./Stan/fits/pattern1_noninf_fit.rda")
bayes_sum = add_ci4int_bayes(p1_non_const)
bayes_sum

summary(p1_non_const)$summary[1:5,]
# shinystan::launch_shinystan(p1_non_const)

# Extract fit
ext_p1 <- extract(p1_non_const)

y_f = ext_p1$y_f
y_m = ext_p1$y_m

f_median = apply(y_f, 2, median)
f_upper  = apply(y_f, 2, quantile, .975)
f_lower  = apply(y_f, 2, quantile, .025)

m_median = apply(y_m, 2, median)
m_upper  = apply(y_m, 2, quantile, .975)
m_lower  = apply(y_m, 2, quantile, .025)

# hist(y_f[,2])
# hist(y_m)
# mean(f_upper-f_lower)

y_all = tibble(y_median, y_lower, y_upper)

yall = tibble(f_median,f_upper ,f_lower ,m_median,m_upper ,m_lower)

pred_bayes = bind_cols(mn_subset, yall) %>% 
  mutate(y_median = ifelse(SEX == "Male", m_median, f_median),
         y_lower  = ifelse(SEX == "Male", m_lower, f_lower),
         y_upper  = ifelse(SEX == "Male", m_upper, f_upper))

ols_plot +
  #geom_point(data = pred_bayes, aes(y = y_median, color = SEX), alpha=0.25, size = 0.5) +
  #geom_point(data = pred_bayes, aes(y  = y_lower, color = SEX), alpha=0.25, size = 0.5) +
  #geom_point(data = pred_bayes, aes(y = y_upper, color = SEX), alpha=0.25, size = 0.5) +
  geom_smooth(data = pred_bayes, aes(y = y_median, color = SEX), size = 0.5,
              method = "lm", fullrange = TRUE, se=F) +
  geom_smooth(data = pred_bayes, aes(y = y_lower, color = SEX),
              se=F, linetype = "dotted", size = 0.5,
              fullrange = TRUE) +
  geom_smooth(data = pred_bayes, aes(y = y_upper, color = SEX),
              se=F, linetype = "dotted", size = 0.5) +
  labs(y = "WISC full scale IQ", x = "PHT pattern concentration")

