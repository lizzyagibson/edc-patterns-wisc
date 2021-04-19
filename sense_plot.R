source("./packages.R")
source("./read_data.R")

pal = brewer.pal(11, "RdBu")[c(2,4,11,10)]

# DATA ####
summary(mn$P1)
summary(whole_cc$P1)

# Sensitivity model ####
fit_sense <- gam(WISC ~ s(P1, by = SEX) + SEX + M_IQ + ALCOHOL + M_EDU + 
                   MARITAL_STATUS + HOME_SCORE + M_AGE + mat_hard, 
                 data = whole_cc)
summary(fit_sense)

pred = as.data.frame(predict.gam(fit_sense, se.fit = T, type="terms"))

const = tibble(P1 = seq(min(whole_cc$P1), max(whole_cc$P1), length=289)) %>% 
  mutate(SEX = "Females",
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
  mutate(predIQ = mean(whole_cc$WISC) + fit.SEX + fit.s.P1..SEXMales + fit.s.P1..SEXFemales,
         se= se.fit.SEX + se.fit.s.P1..SEXMales + se.fit.s.P1..SEXFemales,
         lci = predIQ - 1.96*se,
         uci = predIQ + 1.96*se) %>% 
  mutate(model = "Sensitivity") %>% 
  dplyr::select(-grep("fit", names(.)))

sense_all = whole_cc %>%
  bind_cols(., pred) %>% 
  mutate(extreme = ifelse(SID %in% c(1209, 1229), "Yes", "No")) %>% 
  dplyr::select(WISC, grep("(SEX|P1)", names(.)), extreme) %>% 
  mutate(predIQ = mean(WISC) + fit.SEX + fit.s.P1..SEXMales + fit.s.P1..SEXFemales,
         se= se.fit.SEX + se.fit.s.P1..SEXMales + se.fit.s.P1..SEXFemales,
         lci = predIQ - 1.96*se,
         uci = predIQ + 1.96*se,
         model = "Sensitivity") %>% 
  dplyr::select(-grep("fit", names(.)))

# Just whole data
sense_all %>% 
  ggplot(aes(x = P1, group = interaction(SEX, model),fill = interaction(SEX, model))) +
  geom_point(aes(y = WISC), alpha=0.25, color="grey") +
  geom_point(aes(y = WISC), color="black",
             data = subset(sense_all, extreme =="Yes")) +
  geom_ribbon(aes(ymin = lci,
                  ymax = uci), 
              alpha = 0.5, data = subset(sense_all, SEX == "Males")) +
  geom_line(aes(y = predIQ),
            data = subset(sense_all, SEX == "Males")) + 
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
summary(main_flip)

# Get predicted values
# this is conditional on everything else!
pred_main = as.data.frame(predict(main_fit, se.fit = TRUE, type = c("terms")))

# Get predicted values for females
pred_flip = as.data.frame(predict(main_flip, se.fit = TRUE, type = c("terms")))
as_tibble(pred_flip)

# Combine predicted outcomes with data
main = 
  mn_subset %>%
  bind_cols(., pred_main) %>% 
  bind_cols(., flip.se = (pred_flip$se.fit.P1 + pred_flip$se.fit.P1.SEX)) %>% # just get se for female
  dplyr::select(SID, WISC, grep("(SEX|P)", names(.)), flip.se) %>% 
  mutate(SEX = str_c(SEX, "s"),
         predIQ = fit.P1 + mean(WISC) + fit.SEX + fit.P1.SEX, # ifelse(SEX == "Female", predIQf, predIQm),
         se.male = se.fit.P1 + se.fit.P1.SEX,
         se.female = flip.se,
         se = ifelse(SEX == "Females", se.female, se.male),
         lci = predIQ - 1.96*se,
         uci = predIQ + 1.96*se) %>% 
  dplyr::select(-grep("fit", names(.))) %>% 
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
sense_plot = bind_rows(sense_all, main)

# Figure S1
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
        legend.text = element_text(size = 15)) + scale_fill_nejm() + scale_color_nejm()
#dev.off()

# Missingness
par(mfrow=c(1,1))

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