
# 01b_Fresnel Comparison ####

Fresnel2021 <- 
  "Cleaned Files" %>% 
  list.files(full.names = T) %>% 
  map(read.csv)

Fresnel2020 <- 
  "Cleaned Files_2020" %>% 
  list.files(full.names = T) %>% 
  map(read.csv)

Fresnel2021 %>% map_dbl(nrow) %>% cbind(
Fresnel2020 %>% map_dbl(nrow))

FullFresnel <- 
1:length(Fresnel2021) %>% 
  map(~full_join(Fresnel2021[[.x]], Fresnel2020[[.x]], 
                 by = "Sp", 
                 suffix = c(".2021", ".2020")))

names(FullFresnel) <- c("Bats_In", "Bats_Out", "Mammals_In", "Mammals_Out")

FullFresnel$Bats_In$Betacov.2021 %>% sum(na.rm = T)
FullFresnel$Bats_In$Betacov.2020 %>% sum(na.rm = T)

dir_create("Full Comparison")

names(FullFresnel) %>% map(~write.csv(FullFresnel[[.x]], 
                                      row.names = F, 
                                      file = paste0("Full Comparison/", ., ".csv")))


# Correlation figures ####

library(tidyverse); library(cowplot); library(colorspace); library(patchwork); library(ggpubr); library(ggregplot)
library(ggtext)

dir_create("Figures")

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")

AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])

AlberColours[length(AlberColours)+1:2] <- 
  
  RColorBrewer::brewer.pal(11, AlberPalettes[[4]])[c(2,10)]

AlberColours <- c(AlberColours, Pink = "#FD6396", Blue = "#3C78D8")

Relabel <- c(glue::glue("Network.{1:4}"),
             glue::glue("Hybrid.{1}"), 
             glue::glue("Trait.{1:3}"))

names(Relabel) <- c("R.Po2", "R.Po3",
                    "R.Dal1","R.Far1",
                    "R.Stock1",
                    "R.Gut1", "R.Car3",
                    "R.Alb")

Relabel[intersect(names(Relabel), names(BatModels_IS))] ->
  
  Relabel

# Figure 1: Observed v Predicted panels ####

Fresnel2021[[1]] %>% 
  #gather("Key", "Value", -c(Sp, Betacov, Rank, PropRank, InSample)) %>%
  gather("Key", "Value", starts_with("R.")) %>%
  mutate_at("Key", ~.x %>% recode(!!!Relabel) %>% factor(levels = Relabel)) %>%
  ggplot(aes(Value, Betacov)) + 
  geom_point(alpha = 0.3, colour = AlberColours[[3]]) + 
  geom_smooth(method = glm, 
              method.args = list(family = "binomial"),
              fill = NA, colour = "black") +
  #geom_smooth(method = lm, fill = NA, colour = "black") +
  #geom_smooth(fill = NA, colour = "black") +
  facet_wrap(~Key, nrow = 2) +
  stat_cor(label.y = 1.2,
           aes(label = ..rr.label..)) +
  #coord_fixed() +
  scale_x_continuous(breaks = c(0, 0.5, 1.3)) +
  scale_y_continuous(breaks = c(0:5/5), limits = c(-0, 1.25)) +
  labs(x = "Proportional rank") -> 
  
  SingleCorrelations_2021

Fresnel2021[[1]] %>% 
  mutate(Key = "Multi-model ensemble") %>%
  ggplot(aes(PropRank, Betacov)) + 
  #ggtitle("Model assemblage") +
  geom_point(alpha = 0.6, colour = AlberColours[[3]], 
             position = position_jitter(h = 0.05)) + 
  #geom_smooth(method = lm, fill = NA, colour = "black") +
  geom_smooth(method = glm, 
              method.args = list(family = "binomial"),
              fill = NA, colour = "black") +
  # geom_smooth(fill = NA, colour = "black") +
  #coord_fixed() + 
  lims(x = c(0, 1)) +
  scale_y_continuous(breaks = c(0:5/5), 
                     limits = c(-0.1, 1.25)) +
  stat_cor(label.y = 1.2, method = "spearman",
           aes(label = ..rr.label..)) +
  facet_wrap(~Key) +
  labs(x = "Proportional rank") ->
  
  OverallCorrelations_2021

Fresnel2020[[1]] %>% 
  #gather("Key", "Value", -c(Sp, Betacov, Rank, PropRank, InSample)) %>%
  gather("Key", "Value", starts_with("R.")) %>%
  mutate_at("Key", ~.x %>% recode(!!!Relabel) %>% factor(levels = Relabel)) %>%
  ggplot(aes(Value, Betacov)) + 
  geom_point(alpha = 0.3, colour = AlberColours[[3]]) + 
  geom_smooth(method = glm, 
              method.args = list(family = "binomial"),
              fill = NA, colour = "black") +
  #geom_smooth(method = lm, fill = NA, colour = "black") +
  #geom_smooth(fill = NA, colour = "black") +
  facet_wrap(~Key, nrow = 2) +
  stat_cor(label.y = 1.2,
           aes(label = ..rr.label..)) +
  #coord_fixed() +
  scale_x_continuous(breaks = c(0, 0.5, 1.3)) +
  scale_y_continuous(breaks = c(0:5/5), limits = c(-0, 1.25)) +
  labs(x = "Proportional rank") -> 
  
  SingleCorrelations_2020

Fresnel2020[[1]] %>% 
  mutate(Key = "Multi-model ensemble") %>%
  ggplot(aes(PropRank, Betacov)) + 
  #ggtitle("Model assemblage") +
  geom_point(alpha = 0.6, colour = AlberColours[[3]], 
             position = position_jitter(h = 0.05)) + 
  #geom_smooth(method = lm, fill = NA, colour = "black") +
  geom_smooth(method = glm, 
              method.args = list(family = "binomial"),
              fill = NA, colour = "black") +
  # geom_smooth(fill = NA, colour = "black") +
  #coord_fixed() + 
  lims(x = c(0, 1)) +
  scale_y_continuous(breaks = c(0:5/5), 
                     limits = c(-0.1, 1.25)) +
  stat_cor(label.y = 1.2, method = "spearman",
           aes(label = ..rr.label..)) +
  facet_wrap(~Key) +
  labs(x = "Proportional rank") ->
  
  OverallCorrelations_2020

(SingleCorrelations_2020|OverallCorrelations_2020)/
  (SingleCorrelations_2021|OverallCorrelations_2021) + 
  plot_layout(widths = c(1.35,1)) +
  ggsave("Figures/Obs_Pred_CorrelationsHorizontal.jpeg", 
         units = "mm", width = 325, height = 150)
