pal = brewer.pal(11, "RdBu")[c(2,4,11,10)]

# DATA ####
ewa <- readMat(here::here("./Data/mn2_EWA_sd1.mat"))[[1]] %>% as_tibble() %>% rename(P1 = V1, P2 = V2)
whole = bind_cols(ppp_cov, ewa) %>% left_join(., mn_outcome) %>% rename(HOME_SCORE = HOMETOT)

whole_cc = whole %>% dplyr::select(SID, WISC, P1, SEX, M_IQ, ALCOHOL, M_EDU, M_AGE, mat_hard,
                                 MARITAL_STATUS, HOME_SCORE) %>% drop_na() %>% 
  mutate_at(vars(c(HOME_SCORE, M_AGE, M_IQ)), ~scale(.)[,1])

mn = whole_cc %>% filter(P1 < mean(P1) + 5.5*sd(P1))
mn_flip = mn %>% mutate(SEX = fct_rev(SEX))

summary(mn$P1)
summary(whole$P1)

# Sensitivity model ####
fit_sense <- gam(WISC ~ s(P1, by = SEX) + SEX + M_IQ + ALCOHOL + M_EDU + 
                   MARITAL_STATUS + HOME_SCORE + M_AGE + mat_hard, 
                 data = whole_cc)
summary(fit_sense)

pred = as.data.frame(predict.gam(fit_sense, se.fit = T, type="terms"))
as_tibble(pred)

const = tibble(P1 = seq(min(whole_cc$P1), max(whole_cc$P1), length=289)) %>% 
  mutate(SEX = "Female",
         M_IQ = 0, 
         ALCOHOL = as_factor("No"),
         M_EDU = "≥ High school degree or equivalent", 
         M_AGE = 0, 
         mat_hard = "No",
         MARITAL_STATUS = "Never Married", 
         HOME_SCORE = 0)

pred_f = as.data.frame(predict.gam(fit_sense, newdata = const, se.fit = T, type="terms"))

# This gets smooth CI for females
sense_f = as_tibble(pred_f) %>% 
  bind_cols(., const) %>% 
  dplyr::select(grep("(SEX|P1)", names(.))) %>%
  mutate(predIQ = mean(whole_cc$WISC) + fit.SEX + fit.s.P1..SEXMale + fit.s.P1..SEXFemale,
         se= se.fit.SEX + se.fit.s.P1..SEXMale + se.fit.s.P1..SEXFemale,
         lci = predIQ - 1.96*se,
         uci = predIQ + 1.96*se) %>% 
  mutate(model = "Sensitivity",
         SEX = str_c(SEX, "s")) %>% 
  dplyr::select(-c(fit.SEX, fit.s.P1..SEXMale, fit.s.P1..SEXFemale,
                   se.fit.SEX, se.fit.s.P1..SEXMale, se.fit.s.P1..SEXFemale))

sense_all = whole_cc %>%
  bind_cols(., pred) %>% 
  mutate(extreme = ifelse(SID %in% c(1209, 1229), "Yes", "No")) %>% 
  dplyr::select(WISC, grep("(SEX|P1)", names(.)), extreme) %>% 
  mutate(predIQ = mean(WISC) + fit.SEX + fit.s.P1..SEXMale + fit.s.P1..SEXFemale,
         se= se.fit.SEX + se.fit.s.P1..SEXMale + se.fit.s.P1..SEXFemale,
         lci = predIQ - 1.96*se,
         uci = predIQ + 1.96*se,
         model = "Sensitivity") %>% 
  dplyr::select(-c(fit.SEX, fit.s.P1..SEXMale, fit.s.P1..SEXFemale,
                   se.fit.SEX, se.fit.s.P1..SEXMale, se.fit.s.P1..SEXFemale))

# Just whole data
sense_all %>% 
  ggplot(aes(x = P1, group = interaction(SEX, model),fill = interaction(SEX, model))) +
  geom_point(aes(y = WISC), alpha=0.25, color="grey") +
  geom_point(aes(y = WISC), color="black",
             data = subset(sense_all, extreme =="Yes")) +
  geom_ribbon(aes(ymin = lci,
                  ymax = uci), 
              alpha = 0.5, data = subset(sense_all, SEX == "Male")) +
  geom_line(aes(y = predIQ),
            data = subset(sense_all, SEX == "Male")) + 
  geom_ribbon(aes(ymin = lci,
                  ymax = uci), 
              alpha = 0.5, data = sense_f) +
  geom_line(aes(y = predIQ), data = sense_f) +
  facet_grid(.~SEX, scales = "free_x") +
  labs(y = "WISC Full Scale IQ", x = "Pattern 1 concentration") +
  theme(legend.title = element_blank())

# Main ####
main_flip <- lm(WISC ~ P1 + SEX*P1 + SEX + M_IQ + ALCOHOL + M_EDU  + M_AGE +
                 MARITAL_STATUS + HOME_SCORE + mat_hard
               , data = mn_flip)

main_fit <- lm(WISC ~ P1 + SEX*P1 + SEX + M_IQ + ALCOHOL + M_EDU  + M_AGE +
                MARITAL_STATUS + HOME_SCORE + mat_hard
              , data = mn)
summary(main_fit)
# Get predicted values
# this is conditional on everything else!
pred_main = as.data.frame(predict(main_fit, se.fit = TRUE, type = c("terms")))

# Get predicted values for females
pred_flip = as.data.frame(predict(main_flip, se.fit = TRUE, type = c("terms")))

# Combine predicted outcomes with data
main = 
  mn %>%
  bind_cols(., pred_main) %>% 
  bind_cols(., flip.se = (pred_flip$se.fit.P1 + pred_flip$se.fit.P1.SEX)) %>% # just get se for female
  dplyr::select(SID, WISC, grep("(SEX|P)", names(.)), flip.se) %>% 
  mutate(predIQ = fit.P1 + mean(WISC) + fit.SEX + fit.P1.SEX, # ifelse(SEX == "Female", predIQf, predIQm),
         se.male = se.fit.P1 + se.fit.P1.SEX,
         se.female = flip.se,
         se = ifelse(SEX == "Female", se.female, se.male),
         lci = predIQ - 1.96*se,
         uci = predIQ + 1.96*se) %>% 
  dplyr::select(-c(fit.SEX, fit.P1, fit.P1.SEX, se.female, se.male, flip.se,
                   se.fit.SEX, se.fit.P1.SEX, se.fit.P1)) %>% 
  mutate(model = "Main")

main %>% 
  ggplot(aes(x = P1)) +
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

# Final plot ####
sense_plot = bind_rows(sense_all, main) %>% 
  mutate(SEX = ifelse(SEX == "Male", "Males", "Females"))

#pdf("./Figures/sense_plot.pdf")
sense_plot %>% 
  ggplot(aes(x = P1, group = model, fill = model)) +
  geom_point(aes(y = WISC), alpha=0.25, color="grey") +
  geom_point(aes(y = WISC), color="black",
             data = subset(sense_plot, extreme =="Yes")) +
  geom_line(aes(y = predIQ, color = model), data = sense_f) +
  geom_ribbon(aes(ymin = lci,
                  ymax = uci),
              alpha = 0.5, data = sense_f) +
  geom_ribbon(aes(ymin = lci,
                  ymax = uci), alpha = 0.5,
              data = filter(sense_plot, SEX == "Males" | model == "Main")) +
  geom_line(aes(y = predIQ, color = model),
            data = subset(sense_plot, SEX == "Males" | model == "Main")) +
  facet_grid(.~SEX, scales = "free_x") +
  labs(y = "WISC Full Scale IQ", x = "Pattern 1 concentration") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15))
#dev.off()
