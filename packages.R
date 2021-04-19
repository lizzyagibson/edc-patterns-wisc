options(scipen = 999)
library(haven)
library(tidyverse)
library(RColorBrewer)
library(mgcv)
library(janitor)
library(gt)
library(reshape2)
library(broom)
library(tableone)
library(xtable)
library(GGally)
library(gtsummary)
library(huxtable)
library(rcompanion)
library(R.matlab)
library(MNdata) # Local package
library(ggsci)
library(patchwork)
library(rstan)
library(bayesplot)
library(shinystan)
library(rstanarm)
library(tictoc)
library(recipes)
library(ggridges)
library(bayestestR)
library(olsrr)
library(mvoutlier)
library(outliers)
library(EnvStats)
library(finalfit) 
library(mice)
library(VIM)

theme_set(theme_bw(base_size = 20) + 
            theme(strip.background = element_rect(fill="white"),
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.title = element_text(size = 12),
                  legend.position = "bottom", 
                  legend.text = element_text(size = 7)))

options(
  ggplot2.discrete.color = pal_nejm()(8),
  ggplot2.discrete.fill = pal_nejm()(8))

# Get lower triangle of the correlation matrix
get_lower_tri<-function(x){
  x[upper.tri(x)] <- NA
  
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      if (j == i) {x[i,j] <- NA}
    }
  }
  
  return(x)
}

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
  
  params = names(fit)[grep("(alpha|beta)", names(fit))]
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
