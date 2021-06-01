source("./packages.R")
source("./read_data.R")

library(gridExtra)
library(gtable)
library(grid)
library(cowplot)

names = colnames(mn_ppp[,-1]) %>% 
  str_replace(., "_", "-") %>% as.character() %>% str_pad(., 6)

dat = tibble(x = 1:17, y = 1:17, lbs = names)

#pdf("./Figures/ppp_corr_prez.pdf")
ggcorr(mn_ppp[,-1], method = c("complete.obs", "spearman"),
       label = TRUE, label_size = 2.75, label_round = 2, 
       size = 0, color = "grey50", layout.exp = 2) +
  geom_text(data=dat, aes(x = x, y = y, label=lbs),
            size = 5, nudge_x = 1) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) + 
  theme(legend.text = element_text(size = 12),
        panel.border = element_blank()) +
  labs(fill = "")
#dev.off()

applied <- eh %>% 
  as_tibble() %>% 
  mutate(Pattern = 1:2) %>% 
  gather(key = Chemicals, value = Loadings, -Pattern) %>%
  mutate(Class = as_factor(case_when(grepl("PB", Chemicals) ~ "Parabens",
                                     grepl("^M[A-Z]", Chemicals) ~ "Phthalates",
                                     TRUE ~ "Phenols")),
         Chemicals = str_replace(Chemicals, "_", "-")) %>% 
  mutate(Pattern = paste0("Pattern ", Pattern),
         Chemicals = fct_inorder(Chemicals))

## Figure 1b
min(applied$Loadings)

# pdf("./Figures/prez_eh_loadings_mock.pdf", height = 2.5)
applied %>%
  mutate(Pattern = fct_rev(Pattern)) %>% 
  mutate(Chemicals = fct_rev(Chemicals)) %>% 
  ggplot(aes(x = Pattern, y = Chemicals)) + 
  geom_tile(aes(fill = Loadings), color = "black", size=.5) + 
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(x = "", y = "", title = "") + 
  theme_test(base_size = 20) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        panel.border = element_blank(),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "#ffffffaa", colour = NA),
        legend.position = c(0.4, -2.5), # c(1,0) right bottom, c(1,1) right top.
        #legend.key.size = unit(0.85, "lines"),
        plot.margin = margin(0,0,1,0, "cm")
        )
# dev.off()

# WORD CLOUD ####
library(wordcloud)
library(RColorBrewer)

matblue = "#1976D2"
orange = "#d27519"

pattern1 = applied %>% filter(Pattern == "Pattern 1") %>% 
  mutate(color = ifelse(Class == "Phthalates", matblue, orange)) %>% 
  arrange(desc(Loadings))

text1 = as.character(pattern1$Chemicals)
sqrttimes1 = round((pattern1$Loadings*100))

#pdf("./Figures/cloud1.pdf")
wordcloud(words = text1, freq = sqrttimes1, min.freq = 2,ordered.colors=T,random.order=F,
          rot.per = .3, colors = pattern1$color)
#dev.off()

pattern2 = applied %>% filter(Pattern == "Pattern 2") %>% 
  mutate(color = ifelse(Class == "Phthalates", matblue, orange)) %>% 
  arrange(desc(Loadings))

text2 = as.character(pattern2$Chemicals)
sqrttimes2 = round((pattern2$Loadings*100))

#pdf("./Figures/cloud2.pdf")
wordcloud(words = text2, freq = sqrttimes2, min.freq = 4,ordered.colors=T,random.order=F,
          rot.per = 0.3, colors = pattern2$color)
#dev.off()


# Effect size figure ####
betas = rbind(c(-3.4, -6.5, -0.4,  -3.4, -6.7, -0.4), 
              c(-0.4, -2.4, 1.7, -0.4, -2.5, 1.7), 
              c(2.0, -0.3, 4.2, 2.1, -0.3, 4.5), 
              c(-0.8, -2.7, 1.1, -0.8, -2.8, 1.1))

colnames(betas) = c("trad_beta", "trad_lower", "trad_upper", "bayes_beta", "bayes_lower", "bayes_upper")

models = bind_cols(betas, patt = c("p1f", "p1m", "p2f", "p2m"))

modplot = models %>% pivot_longer(1:6,
                        names_to = c("model", "point"),
                        names_sep = "_") %>% 
  pivot_wider(values_from = value,
              names_from = point) %>% 
  mutate(pattern = str_sub(patt, 1, 2),
         sex = str_sub(patt, 3)) %>% 
  mutate(pattern = ifelse(pattern == "p1", "Pattern 1", "Pattern 2"),
         model = ifelse(model == "bayes", "Bayesian", "Traditional"),
         sex = ifelse(sex == "f", "Female", "Male"))

model = modplot %>% 
  ggplot(aes(x = pattern, 
             y = beta, ymin = lower, ymax = upper, 
             color = model)) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_errorbar(position = position_dodge(width = .5),
                width = 0.25, size = 1.5) +
  geom_point(position = position_dodge(width = .5), size = 4) +
  theme(axis.text.x = element_text(angle=0, hjust=0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.direction = "vertical") +
  scale_color_manual(values = c("#E18727", "#1976D2"))

modellegend = as_ggplot(get_legend(model))

sex = modplot %>% 
  ggplot(aes(x = pattern, 
             y = beta, ymin = lower, ymax = upper, 
             shape = sex)) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_errorbar(position = position_dodge(width = .5),
                width = 0.25, size = 1.5) +
  geom_point(position = position_dodge(width = .5), size = 4) +
  theme(axis.text.x = element_text(angle=0, hjust=0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.direction = "vertical") +
  scale_color_manual(values = c("#E18727", "#1976D2"))

sexlegend = as_ggplot(get_legend(sex))

# "TallPoppy" = "#BC3C29", "DeepCerulean" = "#0072B5",
# "Zest" = "#E18727", "Eucalyptus" = "#20854E",
# "WildBlueYonder" = "#7876B1", "Gothic" = "#6F99AD",
# "Salomie" = "#FFDC91", "FrenchRose" = "#EE4C97"

nolegend = modplot %>% 
  ggplot(aes(x = pattern, 
             y = beta, ymin = lower, ymax = upper, 
             shape = sex, color = model)) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_errorbar(position = position_dodge(width = .5),
                width = 0.25, size = 1.5) +
  geom_point(position = position_dodge(width = .5), size = 4) +
  labs(y = "Beta coefficient with \n95% confidence intervals",
       x = "BN2MF-identified patterns") +
  theme(axis.text.x = element_text(angle=0, hjust=0.5),
        legend.position = "none") +
  scale_color_manual(values = c("#E18727", "#1976D2"))
nolegend

ggdraw(nolegend) +
  draw_plot(sexlegend, x=.4, y=-0.3) +
  draw_plot(modellegend, x=.165, y=-0.3)

#pdf("./Figures/modelplot.pdf", height=5)
ggdraw(nolegend) +
  draw_plot(sexlegend, x=.375, y=-0.24) +
  draw_plot(modellegend, x=.185, y=-0.24)
#dev.off()
