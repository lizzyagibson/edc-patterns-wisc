library(MNdata)
library(tidyverse)
library(R.matlab)
library(rstan)
library(bayesplot)
library(shinystan)
library(rstanarm)
options(mc.cores = parallel::detectCores())

predictors = inner_join(mn_pht[,1:10], mn_phenol[,1:9], by = "SID") %>% drop_na()
iq = mn_outcome %>% dplyr::select(SID, WISC = WSC_CSFS_84)
mn_dat = predictors %>% left_join(., mn_demo, by = "SID") %>% left_join(., iq, by = "SID")
# 32 without WISC

summary(mn_outcome$WSC_CSFS_84)
summary(mn_outcome$WSC_SSFS_84)

e_wa <- readMat(here::here("./Data/mn2_EWA_un.mat"))[[1]] %>% as_tibble() %>% rename(P2 = V1, P1 = V2)
var_wa <- readMat(here::here("./Data/mn2_WA_var.mat"))[[1]] %>% as_tibble() %>% rename(varP2 = V1, varP1 = V2)

bayes = bind_cols(mn_dat, e_wa, var_wa) %>% 
        filter(SEX == "Female") %>% 
        filter(P2 < mean(P2) + 5*sd(P2)) %>% # removes two females
        dplyr::select(-(2:18), -SMOKER_IN_HOME, -M_AGE, -SEX) %>% drop_na()
bayes

# For stan
N = nrow(bayes)
C = ncol(bayes) - 6
K = 2
x = bayes %>% dplyr::select(2:7)
y = bayes$WISC
ewa = bayes[,9:10]
sd_ewa = bayes[,11:12] %>% mutate_all(~sqrt(.))

female_data = list(N = N,
                  C = C,
                  K = 2,
                  x = x,
                  y = y,
                  ewa = ewa,
                  sd_ewa = sd_ewa)

female_fit <- stan(
  file = "bn2mf_score_reg.stan",  # Stan program
    data = female_data,           # named list of data
    chains = 4,                   # number of Markov chains
    warmup = 1000,                # number of warmup iterations per chain
    iter = 5000,                  # total number of iterations per chain
    cores = 2                    # number of cores (could use one per chain)
  )         

betas = c("beta_c[1]","beta_c[2]",
          "beta_c[3]","beta_c[4]",
          "beta_c[5]","beta_c[6]",
          "beta_p[1]","beta_p[2]")

print(female_fit, "alpha")
print(female_fit, betas)
# summary(female_fit)$summary[1:10,]
plot(female_fit, pars = betas)

# Extract fit
ext_fit <- extract(female_fit, inc_warmup = TRUE, permuted = FALSE)
str(ext_fit)
#dimnames(ext_fit)

# RMSE
y_pred = ext_fit[,,275:406]
#dimnames(y_pred)
sqrt(mean((apply(y_pred, 3, median) - y)^2))

traceplot(female_fit, pars = betas, inc_warmup = TRUE)

# shinystan::launch_shinystan(female_fit)

# Plot beta coeffs with uncertainty
mcmc_areas(female_fit,
           pars = betas[7:8],
           prob = 0.8) + ggtitle("Posterior distributions",
                                 "with medians and 80% intervals")

# posterior predict = A draws by nrow(newdata) matrix of simulations from the posterior predictive 
# distribution. Each row of the matrix is a vector of predictions generated using a single draw of 
# the model parameters from the posterior distribution. 
dim(y_pred)

y_pred_df = tibble()
for (i in 1:4) {
  stack = as_tibble(y_pred[,i,])
  y_pred_df = bind_rows(y_pred_df, stack)
  }
dim(y_pred_df)  

rownames(y_pred_df)

samp = sample(1:20000, 75)
dim(y_pred_df[samp,])
post_samp = y_pred_df[samp,]

# plot posterior prediction
color_scheme_set("brightblue")
ppc_dens_overlay(y = as.vector(y),
                 yrep = as.matrix(post_samp))

# posterior pred by grouping variable
as.matrix(post_samp) %>%
  ppc_stat_grouped(y = as.vector(bayes$WISC),
           group = as.vector(bayes$M_EDU),
                   stat = "median")

# same as traceplots above
# color_scheme_set("mix-blue-pink")
# mcmc_trace(ext_fit, pars = betas, n_warmup = 1000,
#                 facet_args = list(nrow = 2, labeller = label_parsed))

# scatterplot
color_scheme_set("darkgray")
mcmc_scatter(
  as.matrix(female_fit),
  pars = betas[7:8],
  np = nuts_params(female_fit),
  np_style = scatter_style_np(div_color = "red", div_alpha = 0.8)
)

# prediction
color_scheme_set("brightblue")
ppc_intervals(
  y = as.vector(y),
  yrep = as.matrix(post_samp),
  x = as.vector(as.matrix(ewa[,1])),
  prob = 0.5
)
