source("./packages.R")
source("./read_data.R")

names = colnames(mn_ppp[,-1]) %>% 
  str_replace(., "_", "-") %>% as.character() %>% str_pad(., 6)

dat = tibble(x = 1:17, y = 1:17, lbs = names)

#pdf("./Figures/ppp_corr_prez.pdf")
ggcorr(mn_ppp[,-1], method = c("complete.obs", "spearman"),
       label = TRUE, label_size = 2.75, label_round = 2, 
       size = 0, color = "grey50", layout.exp = 2) +
  geom_text(data=dat, aes(x = x, y = y, label=lbs),
            size = 5, nudge_x = 1) +
  scale_fill_distiller(palette = "RdYlBu", limits = c(-1,1)) + 
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

# pdf("./Figures/prez_eh_loadings.pdf", height = 2.5)
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

pdf("./Figures/cloud1.pdf")
wordcloud(words = text1, freq = sqrttimes1, min.freq = 2,ordered.colors=T,random.order=F,
          rot.per = .3, colors = pattern1$color)
dev.off()

pattern2 = applied %>% filter(Pattern == "Pattern 2") %>% 
  mutate(color = ifelse(Class == "Phthalates", matblue, orange)) %>% 
  arrange(desc(Loadings))

text2 = as.character(pattern2$Chemicals)
sqrttimes2 = round((pattern2$Loadings*100))

pdf("./Figures/cloud2.pdf")
wordcloud(words = text2, freq = sqrttimes2, min.freq = 4,ordered.colors=T,random.order=F,
          rot.per = 0.3, colors = pattern2$color)
dev.off()

