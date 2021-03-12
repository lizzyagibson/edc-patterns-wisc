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
e_wa <- readMat(here::here("./Data/mn2_EWA_un.mat"))[[1]] %>% as_tibble() %>% rename(P2 = V1, P1 = V2)
var_wa <- readMat(here::here("./Data/mn2_WA_var.mat"))[[1]] %>% as_tibble() %>% rename(varP2 = V1, varP1 = V2)

summary(e_wa)
summary(var_wa)

bayes = bind_cols(mn_ppp, e_wa, var_wa) %>% 
        filter(SEX == "Female") %>% 
        filter(P2 < mean(P2) + 5*sd(P2)) %>% # removes two females
        drop_na()
bayes # 133 = n

# For stan model ####
N = nrow(bayes)
C = ncol(bayes) - 6
K = 2
ewa = bayes %>% dplyr::select(P2:P1)
sd_ewa = bayes %>% dplyr::select(varP2:varP1) %>% mutate_all(~sqrt(.))
x = model.matrix(WISC ~ M_EDU + MARITAL_STATUS + HOME_SCORE + M_IQ + ALCOHOL + M_AGE, data = bayes)[,-1]
y = bayes$WISC

# Run stan model ####
female_data = list(N = N,
                  C = C,
                  K = K,
                  x = x,
                  y = y,
                  ewa = ewa,
                  sd_ewa = sd_ewa)

tic()
female_fit <- stan(
  file = "bn2mf_score_reg.stan",  # Stan program
    data = female_data,           # named list of data
    chains = 4,                   # number of Markov chains
    warmup = 1000,                # number of warmup iterations per chain
    iter = 5000,                  # total number of iterations per chain
    cores = 2                     # number of cores
  )         
toc()

params = tibble(names = names(female_fit)) %>% 
  filter(grepl("(alpha|sigma|beta)", names)) %>% 
  as.matrix()

# Extract fit
ext_fit <- extract(female_fit, inc_warmup = TRUE, permuted = FALSE)
#str(ext_fit)

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

# Summarize model ####

# RMSE
sqrt(mean((apply(y_pred, 3, median) - y)^2))

# n_eff is a crude measure of effective sample size, 
# Rhat is the potential scale reduction factor on split chains (at convergence Rhat=1).
print(female_fit, c("lp__"))
print(female_fit, params)

# plot coefficients
plot(female_fit, pars = params)

# Plot pattern beta coefficient with uncertainty
color_scheme_set("brightblue")
mcmc_areas(female_fit,
           pars = params[(length(params)-1):length(params)],
           prob = 0.95) + ggtitle("Posterior distributions",
                                 "with medians and 95% intervals")

# plot posterior prediction
ppc_dens_overlay(y = as.vector(y),
                 yrep = as.matrix(post_samp))

# posterior pred by grouping variable
as.matrix(post_samp) %>%
  ppc_stat_grouped(y = as.vector(bayes$WISC),
                   group = as.vector(bayes$M_EDU),
                   stat = "median")

# prediction
color_scheme_set("brightblue")
ppc_intervals(
  y = as.vector(y),
  yrep = as.matrix(post_samp),
  x = as.vector(as.matrix(ewa[,1])),
  prob = 0.5
)

# Check model fit ####

# shinystan::launch_shinystan(female_fit)
traceplot(female_fit, pars = params, inc_warmup = TRUE)

# same as traceplots above
# color_scheme_set("mix-blue-pink")
# mcmc_trace(ext_fit, pars = params, n_warmup = 1000,
#                 facet_args = list(nrow = 2, labeller = label_parsed))

# check posterior
lp_female = log_posterior(female_fit)
#head(lp_female)
#dim(lp_female)
np_female <- nuts_params(female_fit)

color_scheme_set("brightblue")
mcmc_nuts_divergence(np_female, lp_female)

# scatterplot
# divergences will be colored in the plot (by default in red).
# this indicate that something is wrong with the model and the results should not be trusted
color_scheme_set("darkgray")
mcmc_scatter(
  as.matrix(female_fit),
  pars = params[7:8],
  np = nuts_params(female_fit),
  np_style = scatter_style_np(div_color = "red", div_alpha = 0.8)
)

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
rhats <- rhat(female_fit)
mcmc_rhat(rhats)

# The larger the ratio of 𝑛𝑒𝑓𝑓 to 𝑁 the better
# a useful heuristic is to worry about any 𝑛𝑒𝑓𝑓/𝑁 less than 0.1.
ratios_cp <- neff_ratio(female_fit)
#names(ratios_cp)[1:11]
mcmc_neff(ratios_cp[names(ratios_cp) %in% params], size = 2)

# check autocorrelation of draws
# autocorrelation for each Markov chain separately up to a user-specified number of lags. 
# Positive autocorrelation is bad (it means the chain tends to stay in the same area between 
# iterations) and you want it to drop quickly to zero with increasing lag. Negative autocorrelation 
# is possible and it is useful as it indicates fast convergence of sample mean towards true mean.
mcmc_acf_bar(female_fit, pars = params, lags = 10)

