
assoc <- read_csv('~/Github/virionette/03_interaction_data/virionette.csv')

BatModels %>% as_tibble %>% mutate(InAssoc = (Sp %in% assoc$host_species)) -> BatModels

library(tidyverse)

BatModels <- read_csv("~/Github/Fresnel_Jun/WeightedEnsemble2021.csv")

#
library(fasterize)
library(rgdal)
library(raster)
library(sf)

iucn <- st_read(dsn = "C:/Users/cjcar/Dropbox/CurrentIUCN", layer='MAMMALS')

r <- disaggregate(getData("worldclim",var="alt",res=2.5)*0,2) # Make a blank raster

######################################################## TOP 50 DE NOVO BETACOV-

BatModels %>% filter(EnsembleBinary==TRUE,
                     Betacov==0) %>% 
              dplyr::pull(Sp) -> 
              unsampled.0

#unsampled.0 <- unsampled.0[1:50] # Top 50 predictions

### Some manual fixes

unsampled.0[unsampled.0=='Myotis oxygnathus'] <- 'Myotis blythii'
unsampled.0[unsampled.0=='Pteropus leucopterus'] <- 'Desmalopex leucopterus'
unsampled.0[unsampled.0=='Chaerephon leucogaster'] <- 'Chaerephon pumilus'
unsampled.0[unsampled.0=='Harpiocephalus mordax'] <- 'Harpiocephalus harpia'
unsampled.0[unsampled.0=='Paracoelops megalotis'] <- 'Hipposideros pomona'
unsampled.0[unsampled.0=='Megaderma lyra'] <- 'Lyroderma lyra'
unsampled.0[unsampled.0=='Falsistrellus affinis'] <- 'Hypsugo affinis'
unsampled.0[unsampled.0=='Rhinolophus imaizumii'] <- 'Rhinolophus perditus'
unsampled.0[unsampled.0=='Eptesicus nasutus'] <- 'Rhyneptesicus nasutus'
unsampled.0[unsampled.0=='Arielulus aureocollaris'] <- 'Thainycteris aureocollaris'
unsampled.0[unsampled.0=='Scotoecus hindei'] <- 'Scotoecus hirundo'
unsampled.0[unsampled.0=='Rousettus bidens'] <- 'Boneia bidens'
unsampled.0[unsampled.0=='Murina silvatica'] <- 'Murina ussuriensis'
unsampled.0[unsampled.0=='Hypsugo bodenheimeri'] <- 'Hypsugo ariel'
unsampled.0[unsampled.0=='Scotoecus albigula'] <- 'Scotoecus hirundo'
unsampled.0[unsampled.0=='Pipistrellus subflavus'] <- 'Perimyotis subflavus'
unsampled.0[unsampled.0=='Artibeus cinereus'] <- 'Dermanura cinerea'
unsampled.0[unsampled.0=='Murina grisea'] <- 'Harpiola grisea'
unsampled.0[unsampled.0=='Neoromicia brunneus'] <- 'Neoromicia brunnea'
unsampled.0[unsampled.0=='Hypsugo imbricatus'] <- 'Pipistrellus imbricatus'
unsampled.0[unsampled.0=='Eptesicus matroka'] <- 'Neoromicia matroka'
unsampled.0[unsampled.0=='Hypsugo crassulus'] <- 'Pipistrellus crassulus'
unsampled.0[unsampled.0=='Hypsugo anchietae'] <- 'Pipistrellus anchietae'
unsampled.0[unsampled.0=='Nyctophilus timoriensis'] <- 'Nyctophilus corbeni'


###

# Pipistrellus tenuis is unmapped
# Pteropus brunneus and pilosus are extinct

iucn.sub <- iucn[iucn$binomial %in% unsampled.0,] # Subset IUCN maps to the right ones

length(unique(iucn.sub$binomial)) # How many of the species made it in?

unsampled.0[!(unsampled.0 %in% iucn$binomial)] # Which names are bonked?

unsampled50.map <- fasterize(iucn.sub, r, fun="sum")

#outline <- rasterToPolygons(r, dissolve=TRUE)

library(maps)
library(rasterVis)
library(RColorBrewer)

par(mfrow=c(3,1))

unsampled50.map <- sum(unsampled50.map, r, na.rm=TRUE)
unsampled50.map <- sum(unsampled50.map, r)

mycolors <- colorRampPalette(rev(brewer.pal(10,"Spectral")))(21)
mycolors[1] <- "#C0C0C0"

rasterVis::levelplot(unsampled50.map,  
                     col.regions = mycolors,
                     #at = seq(0, 20, 1),
                     alpha = 0.5, 
                     scales=list(alternating=FALSE),
                     par.strip.text=list(cex=0),
                     xlab = NULL, ylab = NULL,
                     margin = FALSE, 
                     maxpixels = 5e6,
                     colorkey=list(space="right")) #+


