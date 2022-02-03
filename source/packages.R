options(scipen = 999)

# packages used
list.of.packages <- c( "haven", "tidyverse", "RColorBrewer", "mgcv", "janitor", "gt", "reshape2", "broom", 
                       "tableone", "xtable", "GGally", "gtsummary", "huxtable", "rcompanion", "R.matlab", 
                       "ggsci", "patchwork", "rstan", "bayesplot", "shinystan",  "tictoc", "recipes", 
                       #"rstanarm",
                       "ggridges", "bayestestR", "olsrr", "mvoutlier", "outliers", "EnvStats", "finalfit",  
                       "mice", "VIM")

# if these aren't installed, install them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load all packages
lapply(list.of.packages, library, character.only = TRUE)

# set ggplot theme
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

# quick tidy output with confidence intervals
tidy_ci = function(fit) {
  tidy(fit) %>% bind_cols(., as_tibble(confint(fit)))
}

# add CI for interaction term to tidy output
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

# add distribution for interaction term to Bayesian output
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

# add CI for interaction term to tidy output
add_ci4interaction_both <- function(fit, term1, interaction1, term2, interaction2) {
  # CONFIDENCE INTERVAL for pattern in females
  # Compute association and its uncertainty 
  
  # here we create tables of coefficients and covariance
  coef.mat <- summary(fit)$coefficients
  var.mat  <- vcov(fit)
  
  # the total term for the association is the 
  # sum of P1 in the reference sex plus the term for P1:female
  beta.Pf_1 <- coef.mat[term1,1] + coef.mat[interaction1,1]
  
  # need this for both interactions
  # sum of P2 in the reference sex plus the term for P2:female
  beta.Pf_2 <- coef.mat[term2,1] + coef.mat[interaction2,1]
  
  # Compute variance in order to compute standard error
  # We must compute the variance for the total term s
  # Var(Beta1 + Beta2) = Var(Beta1) + Var(Beta2) + 2*CoVar(Beta1, Beta2) 
  var.Pf_1 <- var.mat[term1, term1] + 
    var.mat[interaction1, interaction1] +
    2*var.mat[term1, interaction1]
  
  var.Pf_2 <- var.mat[term2, term2] + 
    var.mat[interaction2, interaction2] +
    2*var.mat[term2, interaction2]
  
  # this is st error NOT std
  ste.Pf_1  <-  sqrt(abs(var.Pf_1))
  ste.Pf_2  <-  sqrt(abs(var.Pf_2))
  
  # compute confidence intervals 
  lci.Pf_1 <- beta.Pf_1 - 1.96*ste.Pf_1
  uci.Pf_1 <- beta.Pf_1 + 1.96*ste.Pf_1
  
  lci.Pf_2 <- beta.Pf_2 - 1.96*ste.Pf_2
  uci.Pf_2 <- beta.Pf_2 + 1.96*ste.Pf_2
  
  # 2 calculate the test statistic: t = Est/SE
  # 3 calculate the P value2: P = exp(−0.717×z − 0.416×z2).
  test_stat_1 = abs(beta.Pf_1/ste.Pf_1)
  pvalue_1 = exp(-0.717*test_stat_1 - 0.416*test_stat_1^2)
  
  test_stat_2 = abs(beta.Pf_2/ste.Pf_2)
  pvalue_2 = exp(-0.717*test_stat_2 - 0.416*test_stat_2^2)
  
  tidy_ci(fit) %>% 
    bind_rows(., tibble(term = paste0(term1, " in females"), estimate = beta.Pf_1, std.error = ste.Pf_1,
                        statistic = test_stat_1, p.value = pvalue_1, `2.5 %` = lci.Pf_1, `97.5 %` = uci.Pf_1)) %>% 
    bind_rows(., tibble(term = paste0(term2, " in females"), estimate = beta.Pf_2, std.error = ste.Pf_2,
                        statistic = test_stat_2, p.value = pvalue_2, `2.5 %` = lci.Pf_2, `97.5 %` = uci.Pf_2)) %>% 
    mutate(term = ifelse(term == term1, paste(term1, "in males"), term),
           term = ifelse(term == term2, paste(term2, "in males"), term))
}

# add distribution for interaction term to Bayesian output
# FOR MODEL WITH BOTH PATTERNS
add_ci4int_bayes_both = function(fit) {
  # need interaction sum for each pattern
  int_sum_1 = summary(fit, c("beta_int_1", "beta_sex", "beta_p_1"))$summary
  int_sum_2 = summary(fit, c("beta_int_2", "beta_sex", "beta_p_2"))$summary
  
  params = names(fit)[grep("(alpha|beta)", names(fit))]
  model_sum = summary(fit, params)$summary
  
  ext <- extract(fit)
  
  post_coef = cbind(ext$alpha, ext$beta_c, 
                    ext$beta_int_1, ext$beta_int_2,ext$beta_sex, 
                    ext$beta_p_1, ext$beta_p_2) %>% as_tibble()
  
  colnames(post_coef) = params
  # interaction coefficients for each pattern
  int_coefs_1 = post_coef %>% dplyr::select(beta_p_1, beta_int_1)
  int_coefs_2 = post_coef %>% dplyr::select(beta_p_2, beta_int_2)
  
  # all steps for both patterns
  coef1_p1 = int_coefs_1$beta_p_1
  coef2_p1 = int_coefs_1$beta_int_1
  
  coef1_p2 = int_coefs_2$beta_p_2
  coef2_p2 = int_coefs_2$beta_int_2
  
  # get covariance and beta in females for both patterns
  covmat_1 = cov(int_coefs_1)
  beta_f_1 = coef1_p1 + coef2_p1
  
  covmat_2 = cov(int_coefs_2)
  beta_f_2 = coef1_p2 + coef2_p2

  # get variance in females for both patterns
  var_f_1 = covmat_1[1, 1] + 
    covmat_1[2, 2] +
    2*covmat_1[1, 2]
  
  var_f_2 = covmat_2[1, 1] + 
    covmat_2[2, 2] +
    2*covmat_2[1, 2]

  # this is std dev not std error
  # get sd for each pattern in females
  sd_f_1  <- sqrt(abs(var_f_1))
  sd_f_2  <- sqrt(abs(var_f_2))
  
  # take min sample size as conservative est.
  # sample size for each patter
  sampsize_1 = min(int_sum_1[1,9], int_sum_1[3,9])
  sampsize_2 = min(int_sum_2[1,9], int_sum_2[3,9])
  
  # this is std error
  # std err for both patterns in females
  se_f_1  <- sd_f_1/sqrt(sampsize_1)
  se_f_2  <- sd_f_2/sqrt(sampsize_2)
  
  # get confidence intervals for each
  lci.f_1 <- mean(beta_f_1) - 1.96*se_f_1
  uci.f_1 <- mean(beta_f_1) + 1.96*se_f_1
  
  lci.f_2 <- mean(beta_f_2) - 1.96*se_f_2
  uci.f_2 <- mean(beta_f_2) + 1.96*se_f_2
  
  # 2 calculate the test statistic: t = Est/SE
  # 3 calculate the P value2: P = exp(−0.717×z − 0.416×z2).
  # test_stat = abs(mean(beta_f)/se_f)
  # pvalue = exp(-0.717*test_stat - 0.416*test_stat^2)
  # NO BAYESIAN PVALUE
  
  # get quanties for both patterns in females
  q_1 = quantile(beta_f_1, probs = c(0.025, .25, .50, .75, .975)) %>% as.data.frame() %>% 
    rownames_to_column() %>% 
    pivot_wider(names_from = "rowname",
                values_from = ".")
  
  q_2 = quantile(beta_f_2, probs = c(0.025, .25, .50, .75, .975)) %>% as.data.frame() %>% 
    rownames_to_column() %>% 
    pivot_wider(names_from = "rowname",
                values_from = ".")
  
  var_1 = tibble(beta = "Pattern 1 in females", mean = round(mean(beta_f_1), 4), 
                   se_mean = round(se_f_1,4), sd = round(sd_f_1,4)) %>% 
    bind_cols(., q_1)
  
  var_2 = tibble(beta = "Pattern 2 in females", mean = round(mean(beta_f_2), 4), 
                   se_mean = round(se_f_2,4), sd = round(sd_f_2,4)) %>% 
    bind_cols(., q_2)
  
  # M_AGE + M_EDU + MARITAL_STATUS + HOME_SCORE + M_IQ + ALCOHOL + mat_hard, 
  as.data.frame(model_sum)[,-c(9:10)] %>% 
    rownames_to_column(var="beta") %>% 
    mutate(beta = c("alpha", "age", "edu", "marital", "home", "iq", "alcohol", "mat_hard",
                    "Sex*Pattern 1", "Sex*Pattern 2", "Sex", "Pattern 1 in males", "Pattern 2 in males")) %>%
    bind_rows(.,var_1, var_2) %>% 
    mutate_at(vars(2:9), round, 2)
}
