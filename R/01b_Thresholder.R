

####### RUN THIS TO SKIP IMPORT STEP

# # Importing predictions and creating ranked predictions ####
# 
# library(tidyverse); library(fs); library(magrittr); library(ggregplot)
# library(conflicted); library(dplyr)
# 
# conflict_prefer("map", "purrr")
# 
# c("select", "filter", "intersect", "summarise", "mutate", "rename", "arrange") %>%
#   map(~conflict_prefer(.x, "dplyr"))
# 
# rm(list = ls())
# 
# here::here() %>% setwd()

# Thresholding Model Predictions ####

library(PresenceAbsence); library(tidyverse)

conflict_prefer("extract", "magrittr")

BatModels <- "Cleaned Files" %>% list.files(full.names = T) %>% extract(1:2) %>% 
  map(read.csv) %>% bind_rows()

RNames <- c('R.Alb','R.Car3','R.Dal1','R.Far1','R.Gut1','R.Po2','R.Po3','R.Stock1')

RNames %>% str_replace_all("^R.", "P.") -> PNames

thresh <- c('T.Alb','T.Car','T.Dal','T.Far','T.Gut','T.Po1','T.Po2','T.Stock1')

negatory <- function(x) {1-x}

BatModels %>% mutate(n = 1:nrow(BatModels)) %>%
  mutate_at(RNames, negatory) %>%
  mutate_at('PropRank', negatory) -> BatModels2

####### GET AUC'S

for (i in 1:8) {
  
  print(c(PNames)[i])
  print(auc(data.frame(BatModels2[,c('n','Betacov',PNames[i])]), na.rm = TRUE))
  
}

auc(data.frame(BatModels2[,c('n','Betacov',"PropRank")]), na.rm = TRUE)

####### 

tvalues <- sapply(c(RNames,PNames,'PropRank'), function(x){
  o <- optimal.thresholds(data.frame(BatModels2[,c('n','Betacov',x)]),
                          threshold = 10001,
                          opt.methods = 10,
                          req.sens = 0.9,
                          na.rm = TRUE)
  return(o[,2])
})

for (name in RNames) {BatModels2[,name] <- as.vector(BatModels2[,name] > tvalues[name])}
for (name in PNames) {BatModels2[,name] <- as.vector(BatModels2[,name] > tvalues[name])}

colSums(BatModels2[BatModels2$Betacov==0,RNames], na.rm = TRUE)
colSums(BatModels2[BatModels2$Betacov==0,PNames], na.rm = TRUE)

# TOTAL RANK

BatModels2[,'PropRank'] <- as.vector(BatModels2[,'PropRank'] > tvalues['PropRank'])

table(BatModels2[BatModels2$Betacov==0,'PropRank'])

# Rhinolophus

# BatModels2[grep('Rhinolophus',BatModels2$Sp),] %>% View()

BatModels2[grep('Rhinolophus',BatModels2$Sp),] %>%
  filter(Betacov==0) %>% select(PropRank) %>% table()

# Clean it up to write out

BatModels2 %>% select(Sp, Betacov, P.Alb, P.Car3, P.Dal1, P.Far1, P.Gut1, P.Po2, P.Po3, P.Stock1, PropRank) %>%
  rename(Trait.1 = P.Gut1,
         Trait.2 = P.Car3,
         Trait.3 = P.Alb,
         Network.1 = P.Po2,
         Network.2 = P.Po3,
         Network.3 = P.Dal1,
         Network.4 = P.Far1,
         Hybrid.1 = P.Stock1,
         Ensemble = PropRank) -> BatModels2 

BatModels2 %>% write_csv("BinaryPredictions.csv")

# WHO WON! WHO'S NEXT!

verify <- c('Acerodon celebensis',
            'Artibeus jamaicensis',
            'Carollia sowelli',
            'Chaerephon pumilus', 
            'Desmodus rotundus',
            'Epomops buettikoferi',
            'Emballonura alecto',
            'Glauconycteris variegata', 
            'Hipposideros cervinus',
            'Hipposideros fuliginosus', 
            'Hipposideros gigas',
            'Hipposideros larvatus', 
            'Hipposideros lekaguli',
            'Hipposideros pomona', 
            'Hypsugo pulveratus',
            'Macroglossus minimus',
            'Megaerops ecaudatus',
            'Megaerops kusnotoi',
            'Miniopterus magnater',
            'Myonycteris angolensis',
            'Myonycteris torquata', 
            'Myotis horsfieldii',
            'Myotis pequinius', 
            'Myotis punicus',
            'Nanonycteris veldkampii',
            'Neoromicia somalicus',
            'Neoromicia zuluensis',
            'Nycteris gambiensis',
            'Nycteris macrotis',
            'Nycteris thebaica',
            'Pipistrellus coromandra',
            'Pipistrellus deserti',
            'Pipistrellus tenuis',
            'Plecotus auritus',
            'Pteronotus personatus',
            'Pteropus conspicillatus',
            'Pteropus lylei', 
            'Rhinolophus acuminatus',
            'Rhinolophus malayanus',
            'Rhinolophus rufus',
            'Rhinolophus shameli',
            'Rhinolophus stheno',
            'Rousettus madagascariensis',
            'Scotophilus heathii', 
            'Scotophilus kuhlii',
            'Tadarida teniotis',
            'Vespertilio murinus'
)
# BatModels2 %>% filter(Sp %in% verify) %>% View()

#################################

key1 <- c(`1` = "Reported", `0` = "Unreported")

key2 <- c(`0` = "Unlikely", `1` = "Suspected",
          `2` = "False -", `3` = "True +")

BatModels2 %>% mutate_at(vars(contains(".")), ~(as.numeric(.))) %>%
  mutate(Ensemble = as.numeric(Ensemble)) %>%
  mutate(Betacov2 = Betacov*2) %>%
  mutate_at(vars(contains(".")), ~(. + Betacov2)) %>% 
  mutate(Ensemble = Ensemble + Betacov2) %>% 
  mutate_at(vars(contains(".")), ~(recode(.,!!!key2))) %>% 
  mutate(Ensemble =  recode(Ensemble,!!!key2)) %>% 
  mutate(Betacov = recode(Betacov, !!!key1)) %>% 
  mutate(`New data` = Betacov) %>% 
  rename(`Training data` = Betacov) -> BatWeb

BatWeb$`New data`[which(BatWeb$Sp %in% verify)] <- 'New data'

BatWeb %>% select(c(Sp,
                    `Training data`, 
                    `New data`,
                    Ensemble,
                    Trait.1,
                    Trait.2,
                    Trait.3,
                    Hybrid.1,
                    Network.1,
                    Network.2,
                    Network.3,
                    Network.4)) %>% rename(Hybrid = Hybrid.1) -> BatWeb

BatWeb %>% mutate(Source = '') -> BatWeb

BatWeb$Source[BatWeb$Sp %in% c('Hipposideros pomona',
                               'Scotophilus kuhlii', 
                               'Myotis pequinius', 
                               'Myotis horsfieldii')] <- 'https://www.biorxiv.org/content/10.1101/2020.05.31.116061v1'

BatWeb$Source[BatWeb$Sp %in% c('Pteropus lylei')] <- 'https://virologyj.biomedcentral.com/articles/10.1186/s12985-018-0950-6'

BatWeb$Source[BatWeb$Sp %in% c('Hipposideros larvatus',
                               'Scotophilus heathii',
                               'Hipposideros lekaguli')] <- 'https://virologyj.biomedcentral.com/articles/10.1186/s12985-015-0289-1'

BatWeb$Source[BatWeb$Sp %in% c('Desmodus rotundus')] <- 'https://www.scielo.br/scielo.php?script=sci_arttext&pid=S1413-86702008000600003'

BatWeb$Source[BatWeb$Sp %in% c('Macroglossus minimus')] <- 'http://philjournalsci.dost.gov.ph/96-next-issue/vol-149-no-1-march-2020/1160-first-molecular-evidence-for-bat-betacoronavirus-in-mindanao'

BatWeb$Source[BatWeb$Sp %in% c('Hipposideros gigas')] <- 'https://www.nature.com/articles/s41598-020-64159-1.pdf'

BatWeb$Source[BatWeb$Sp %in% c('Plecotus auritus',
                               'Tadarida teniotis')] <- 'https://link.springer.com/content/pdf/10.1007/s11262-018-1614-8.pdf'

BatWeb$Source[BatWeb$Sp %in% c('Artibeus jamaicensis',
                               'Carollia sowelli')] <- 'https://onlinelibrary.wiley.com/doi/epdf/10.1111/tbed.13751'

BatWeb$Source[BatWeb$Sp %in% c('Myonycteris angolensis',
                               'Nycteris macrotis',
                               'Nanonycteris veldkampii')] <- 'https://www.mdpi.com/1999-4915/12/8/855'

BatWeb$Source[BatWeb$Sp %in% c('Pipistrellus deserti')] <- 'https://wwwnc.cdc.gov/eid/article/22/1/15-1397_article'

BatWeb$Source[BatWeb$Sp %in% c('Pipistrellus coromandra')] <- 'https://www.sciencedirect.com/science/article/abs/pii/S1567134816305135'

BatWeb$Source[BatWeb$Sp %in% c('Hipposideros lekaguli')] <- 'https://www.sciencedirect.com/science/article/abs/pii/S1567134816305135'

# BatWeb$Source[BatWeb$Sp %in% c('Megaerops kusnotei')] <- 'https://link.springer.com/content/pdf/10.1007/s12250-016-3727-3.pdf'

# BatWeb$Source[BatWeb$Sp %in% c('Hipposideros lekaguli')] <- 'https://www.sciencedirect.com/science/article/abs/pii/S1567134816305135'

BatWeb$Source[BatWeb$Sp %in% c('Megaerops kusnotei')] <- 'https://link.springer.com/content/pdf/10.1007/s12250-016-3727-3.pdf'

BatWeb$Source[BatWeb$Sp %in% c('Hypsugo pulveratus')] <- 'GenBank: MN312842, MN312848, MN312849, MN312852, MN312853, MN312854'

BatWeb$Source[BatWeb$Sp %in% c('Myotis punicus')] <- 'GenBank: MN823619'

BatWeb$Source[BatWeb$Sp %in% c('Rhinolophus shameli')] <- 'https://www.biorxiv.org/content/10.1101/2021.01.26.428212v1'

BatWeb$Source[BatWeb$Sp %in% c('Rhinolophus acuminatus')] <- 'https://www.nature.com/articles/s41467-021-21240-1'

BatWeb$Source[BatWeb$Sp %in% c('Nycteris gambiensis',
                               'Pteronotus personatus')] <- 'https://journals.plos.org/plospathogens/article/authors?id=10.1371/journal.ppat.1008758'

BatWeb$Source[BatWeb$Sp %in% c('Rhinolophus stheno',
                               'Rhinolophus malayanus')] <- 'https://www.biorxiv.org/content/10.1101/2021.03.08.434390v1?rss=1'

BatWeb$Source[BatWeb$Sp %in% c('Vespertilio murinus')] <- 'https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0252534'

BatWeb$Source[BatWeb$Sp %in% c('Myonycteris torquata', 
                               'Hipposideros fuliginosus', 
                               'Hipposideros cervinus',
                               'Chaerephon pumilus', 
                               'Glauconycteris variegata', 
                               'Neoromicia somalicus',
                               'Megaerops ecaudatus',
                               'Epomops buettikoferi',
                               'Acerodon celebensis',
                               'Pteropus conspicillatus')] <- 'USAID PREDICT data - PREDICT_PCR_Tests.csv - June 28, 2021'

BatWeb$Source[BatWeb$Sp %in% c('Emballonura alecto')] <- 'https://link.springer.com/article/10.1007/s00705-012-1410-z'

BatWeb$Source[BatWeb$Sp %in% c('Miniopterus magnater')] <- 'https://www.microbiologyresearch.org/content/journal/jgv/10.1099/vir.0.82203-0#tab2'

BatWeb$Source[BatWeb$Sp %in% c('Neoromicia zuluensis')] <- 'https://wwwnc.cdc.gov/eid/article/19/10/13-0946_article'

BatWeb$Source[BatWeb$Sp %in% c('Nycteris thebaica',
                               'Rousettus madagascariensis')] <- 'https://www.nature.com/articles/s41598-020-63799-7'

BatWeb$Source[BatWeb$Sp %in% c('Pipistrellus tenuis')] <- 'https://journals.asm.org/doi/10.1128/JVI.00116-18?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed'

BatWeb$Source[BatWeb$Sp %in% c('Rhinolophus rufus')] <- 'https://link.springer.com/article/10.1007/s00705-012-1410-z'


# Compile and print out 

# BatWeb %>% write_csv("BinaryWebsite.csv")

BatWeb %>% write_csv("BinaryWebsiteNew.csv")

# Writing out 2020 predictions ####

conflict_prefer("extract", "magrittr")

BatModels <- "Cleaned Files_2020" %>% list.files(full.names = T) %>% extract(1:2) %>% 
  map(read.csv) %>% bind_rows()

RNames <- c('R.Alb','R.Car3','R.Dal1','R.Far1','R.Gut1','R.Po2','R.Po3','R.Stock1')

RNames %>% str_replace_all("^R.", "P.") -> PNames

thresh <- c('T.Alb','T.Car','T.Dal','T.Far','T.Gut','T.Po1','T.Po2','T.Stock1')

negatory <- function(x) {1-x}

BatModels %>% mutate(n = 1:nrow(BatModels)) %>%
  mutate_at(RNames, negatory) %>%
  mutate_at('PropRank', negatory) -> BatModels2

####### GET AUC'S

for (i in 1:8) {
  
  print(c(PNames)[i])
  print(auc(data.frame(BatModels2[,c('n','Betacov',PNames[i])]), na.rm = TRUE))
  
}

auc(data.frame(BatModels2[,c('n','Betacov',"PropRank")]), na.rm = TRUE)

####### 

tvalues <- sapply(c(RNames,PNames,'PropRank'), function(x){
  o <- optimal.thresholds(data.frame(BatModels2[,c('n','Betacov',x)]),
                          threshold = 10001,
                          opt.methods = 10,
                          req.sens = 0.9,
                          na.rm = TRUE)
  return(o[,2])
})

for (name in RNames) {BatModels2[,name] <- as.vector(BatModels2[,name] > tvalues[name])}
for (name in PNames) {BatModels2[,name] <- as.vector(BatModels2[,name] > tvalues[name])}

colSums(BatModels2[BatModels2$Betacov==0,RNames], na.rm = TRUE)
colSums(BatModels2[BatModels2$Betacov==0,PNames], na.rm = TRUE)

# TOTAL RANK

BatModels2[,'PropRank'] <- as.vector(BatModels2[,'PropRank'] > tvalues['PropRank'])

table(BatModels2[BatModels2$Betacov==0,'PropRank'])

# Rhinolophus

# BatModels2[grep('Rhinolophus',BatModels2$Sp),] %>% View()

BatModels2[grep('Rhinolophus',BatModels2$Sp),] %>%
  filter(Betacov==0) %>% select(PropRank) %>% table()

# Clean it up to write out

BatModels2 %>% select(Sp, Betacov, P.Alb, P.Car3, P.Dal1, P.Far1, P.Gut1, P.Po2, P.Po3, P.Stock1, PropRank) %>%
  rename(Trait.1 = P.Gut1,
         Trait.2 = P.Car3,
         Trait.3 = P.Alb,
         Network.1 = P.Po2,
         Network.2 = P.Po3,
         Network.3 = P.Dal1,
         Network.4 = P.Far1,
         Hybrid.1 = P.Stock1,
         Ensemble = PropRank) -> BatModels2 

BatModels2 %>% write_csv("BinaryPredictions_2020.csv")
