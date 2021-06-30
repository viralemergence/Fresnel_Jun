
library(PresenceAbsence)
library(tidyverse)

os <- read_csv("Cleaned Files/BatModels_OS.csv")
is <- read_csv("Cleaned Files/BatModels_iS.csv")

sum(os$Betacov)
sum(is$Betacov)
BatModels <- bind_rows(is, os)

rnames <- c("R.Alb", "R.Car3", "R.Dal1", "R.Far1", "R.Gut1", "R.Po2", "R.Po3", "R.Stock1")

# This isn't sustainable but it's a one time thing.
w1 <- read_rds("C:/Users/cjcar/Documents/Github/Fresnel/EnsembleWeights.rds")

ensemble2 <- matrixStats::rowWeightedMeans(as.matrix(BatModels[,rnames]), 
                                           w = w1, 
                                           na.rm = TRUE)

BatModels %>%
  select(Sp, Betacov) %>%
  mutate(Ensemble.2 = ensemble2) -> df

df %>%
  mutate(Negative = 1 - Ensemble.2) %>%
  select(Sp, Betacov, Negative, Ensemble.2) -> threshdf

tvalues <- optimal.thresholds(threshdf[,c("Sp","Betacov","Negative")],
                              threshold = 10001,
                              opt.methods = 10,
                              req.sens = 0.9,
                              na.rm = TRUE)

threshdf %>% mutate(EnsembleBinary = (Ensemble.2 < (1 - tvalues$Negative))) %>%
  select(-Negative) %>%
  write_csv("WeightedEnsemble2021.csv")
