
# OLS ####
load(file = "reg_pat1_fit.rda")

fit_flip <- lm(WISC ~ P1 + SEX*P1 + SEX + M_IQ + ALCOHOL + M_EDU  + M_AGE +
                 MARITAL_STATUS + HOME_SCORE + mat_hard
               , data = mn_flip)
summary(fit_flip)

# Get predicted values
# this is conditional on everything else!
pred = as.data.frame(predict(fit_p1a, se.fit = TRUE, type = c("terms")))

# Get predicted values for females
pred_flip = as.data.frame(predict(fit_flip, se.fit = TRUE, type = c("terms")))

# Combine predicted outcomes with data
main = 
  mn_subset %>%
  bind_cols(., pred) %>% 
  bind_cols(., flip.se = (pred_flip$se.fit.P1 + pred_flip$se.fit.P1.SEX)) %>% # just get se for female
  dplyr::select(SID, WISC, grep("(SEX|P)", names(.)), flip.se) %>% 
   # -c(M_IQ:model, df, residual.scale, fit.M_IQ:fit.mat_hard, se.fit.M_IQ:se.fit.mat_hard)) %>% 
  mutate(predIQ = fit.P1 + mean(WISC) + fit.SEX + fit.P1.SEX, # ifelse(SEX == "Female", predIQf, predIQm),
         se.male = se.fit.P1 + se.fit.P1.SEX,
         se.female = flip.se,
         se = ifelse(SEX == "Female", se.female, se.male),
         lci = predIQ - 1.96*se,
         uci = predIQ + 1.96*se) %>% 
  dplyr::select(-c(fit.SEX, fit.P1, fit.P1.SEX, se.female, se.male, flip.se,
                   se.fit.SEX, se.fit.P1.SEX, se.fit.P1)) %>% 
  mutate(model = "Main")

main
save(main, file = "./sense_plot.rda")

main %>% 
  ggplot(aes(x = P)) +
  geom_point(aes(y = WISC), color = "lightgray", size = 0.5) + 
  geom_abline(intercept = mean(main$WISC), slope = 0, 
              color = "darkgray", linetype = "dashed") + 
  geom_line(aes(y = predIQ)) +
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

