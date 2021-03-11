# Load packages ####
library(MNdata)
library(tidyverse)
library(R.matlab)
library(rstan)
library(bayesplot)
library(shinystan)
library(rstanarm)
library(tictoc)
library(recipes)
options(mc.cores = parallel::detectCores())

# Data ####
predictors = inner_join(mn_pht[,1:10], mn_phenol[,1:9], by = "SID") %>% drop_na()
iq = mn_outcome %>% dplyr::select(SID, WISC = WSC_CSFS_84)
mn_dat = predictors %>% left_join(., mn_demo, by = "SID") %>% left_join(., iq, by = "SID")
# 32 without WISC

summary(mn_outcome$WSC_CSFS_84)
summary(mn_outcome$WSC_SSFS_84)
summary(mn_demo$M_IQ)

e_wa <- readMat(here::here("./Data/mn2_EWA_un.mat"))[[1]] %>% as_tibble() %>% rename(P2 = V1, P1 = V2)
var_wa <- readMat(here::here("./Data/mn2_WA_var.mat"))[[1]] %>% as_tibble() %>% rename(varP2 = V1, varP1 = V2)

summary(e_wa)
summary(var_wa)

bayes_int = bind_cols(mn_dat, e_wa, var_wa) %>% 
        filter(P2 < mean(P2) + 5*sd(P2)) %>% # removes two females, 1 male
        dplyr::select(-(2:18), -SMOKER_IN_HOME , -M_AGE) %>% drop_na()
bayes_int # 286 = n

# For stan model ####
N = nrow(bayes_int)
C = ncol(bayes_int) - 6 - 1
K = 2
x = model.matrix(WISC ~ ETH + M_EDU + MARITAL_STATUS + HOME_SCORE + M_IQ + ALCOHOL, data = bayes_int)[,-1]
sex = model.matrix(WISC ~ SEX, data = bayes_int)[,-1]
y = bayes_int$WISC
ewa = bayes_int %>% dplyr::select(P2:P1)
sd_ewa = bayes_int %>% dplyr::select(varP2:varP1) %>% mutate_all(~sqrt(.))

# Run stan model ####
data_int = list(N = N,
                  C = C,
                  K = K,
                  x = x,
                  y = y,
                  sex = sex,
                  ewa = ewa,
                  sd_ewa = sd_ewa)

tic()
fit_int <- stan(
  file = "bn2mf_score_reg_int.stan",  # Stan program
    data = data_int,           # named list of data
    chains = 4,                   # number of Markov chains
    warmup = 1000,                # number of warmup iterations per chain
    iter = 5000,                  # total number of iterations per chain
    cores = 2                     # number of cores
  )         
toc()

# Check model fit ####

params_int = tibble(names = names(fit_int)) %>% 
  filter(grepl("(alpha|sigma|beta)", names)) %>% 
  as.matrix()

# shinystan::launch_shinystan(fit_int)
traceplot(fit_int, pars = params_int, inc_warmup = TRUE)

# check posterior
lp_female = log_posterior(fit_int)
#head(lp_female)
#dim(lp_female)
np_female <- nuts_params(fit_int)

color_scheme_set("brightblue")
mcmc_nuts_divergence(np_female, lp_female)

# scatterplot
# divergences will be colored in the plot (by default in red).
# this indicate that something is wrong with the model and the results should not be trusted
color_scheme_set("darkgray")
mcmc_scatter(
  as.matrix(fit_int),
  pars = params_int[7:8],
  np = nuts_params(fit_int),
  np_style = scatter_style_np(div_color = "red", div_alpha = 0.8)
)

# Extract fit
ext_fit <- extract(fit_int, inc_warmup = TRUE, permuted = FALSE)
#str(ext_fit)

# same as traceplots above
# color_scheme_set("mix-blue-pink")
# mcmc_trace(ext_fit, pars = params_int, n_warmup = 1000,
#                 facet_args = list(nrow = 2, labeller = label_parsed))

pred_all = ext_fit[1001:5000,,]
str(pred_all)

y_pred = pred_all[,,grepl("y_tilde", dimnames(pred_all)$parameters)]
#dim(y_pred)
#dimnames(y_pred)

# posterior predict = A draws by nrow(newdata) matrix of simulations from the posterior predictive 
# distribution. Each row of the matrix is a vector of predictions generated using a single draw of 
# the model parameters from the posterior distribution. 
#dim(y_pred)

y_pred_df = tibble()
for (i in 1:4) {
  stack = as_tibble(y_pred[,i,])
  y_pred_df = bind_rows(y_pred_df, stack)
}
#dim(y_pred_df)  

samp = sample(1:16000, 75) # (5000-1000 warmup) * 4
#dim(y_pred_df[samp,])
post_samp = y_pred_df[samp,]
#dim(post_samp)

# error
color_scheme_set("brightblue")
ppc_error_hist(y = as.vector(y),
               yrep = as.matrix(post_samp)[1:6,])

# log posterior
ggplot(lp_female, aes(x = Iteration, y = Value, color = as.factor(Chain))) +
  geom_line(alpha = 0.5) +
  geom_smooth(method = "loess") + 
  labs(y = "Log Posterior Value",
       color = "Chains") +
  theme_bw() +
  theme(legend.position = "bottom")

# 1 is good
# if all chains are at equilibrium, these will be the same and 𝑅̂  will be one. 
# If the chains have not converged to a common distribution, the 𝑅̂  statistic 
# will be greater than one 
rhats <- rhat(fit_int)
mcmc_rhat(rhats)

# The larger the ratio of 𝑛𝑒𝑓𝑓 to 𝑁 the better
# a useful heuristic is to worry about any 𝑛𝑒𝑓𝑓/𝑁 less than 0.1.
ratios_cp <- neff_ratio(fit_int)
#names(ratios_cp)[1:11]
mcmc_neff(ratios_cp[names(ratios_cp) %in% params_int], size = 2)

# check autocorrelation of draws
# autocorrelation for each Markov chain separately up to a user-specified number of lags. 
# Positive autocorrelation is bad (it means the chain tends to stay in the same area between 
# iterations) and you want it to drop quickly to zero with increasing lag. Negative autocorrelation 
# is possible and it is useful as it indicates fast convergence of sample mean towards true mean.
mcmc_acf_bar(fit_int, pars = params_int, lags = 10)

# Summarize model ####

# RMSE
sqrt(mean((apply(y_pred, 3, median) - y)^2))

# n_eff is a crude measure of effective sample size, 
# Rhat is the potential scale reduction factor on split chains (at convergence Rhat=1).
print(fit_int, c("lp__"))
print(fit_int, params_int)

# plot coefficients
plot(fit_int, pars = params_int)

# Plot pattern beta coefficient with uncertainty
color_scheme_set("brightblue")
mcmc_areas(fit_int,
           pars = params_int[9:12],
           prob = 0.95) + ggtitle("Posterior distributions",
                                 "with medians and 95% intervals")

# plot posterior prediction
ppc_dens_overlay(y = as.vector(y),
                 yrep = as.matrix(post_samp))

# posterior pred by grouping variable
as.matrix(post_samp) %>%
  ppc_stat_grouped(y = as.vector(bayes_int$WISC),
                   group = as.vector(bayes_int$M_EDU),
                   stat = "median")

# prediction
color_scheme_set("brightblue")
ppc_intervals(
  y = as.vector(y),
  yrep = as.matrix(post_samp),
  x = as.vector(as.matrix(ewa[,1])),
  prob = 0.5
)

# Plot interaction
y_median = apply(y_pred_df, 2, median)
y_lower = apply(y_pred_df, 2, quantile, 0.025)
y_upper = apply(y_pred_df, 2, quantile, 0.975)

int_viz = bind_cols(ewa, y_median = y_median, upper = y_upper, lower = y_lower, Sex = sex, WISC = y)
print(fit_int, params_int[9:12,])

int_viz %>% 
  mutate(Sex = as_factor(Sex)) %>% 
  ggplot(aes(x = P2, color = Sex, fill = Sex)) +
    geom_point(aes(y = WISC)) +
    geom_smooth(aes(y = y_median,
                    ymin = lower,
                    ymax = upper), method = "lm")
