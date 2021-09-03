## phylofactor betacov ensemble V2*
## danbeck@ou.edu

## clean environment & plots

library(tidyverse); library(fs); library(magrittr); library(ggregplot)
library(conflicted); library(dplyr)

conflict_prefer("map", "purrr")

c("select", "filter", "intersect", "summarise", "mutate", "rename", "arrange") %>%
  map(~conflict_prefer(.x, "dplyr"))

rm(list = ls())

setwd("~/Desktop/Fresnel_Jun")

## source helper files for gpf

source('R/02a_Phylofactor helper.R')

#detach(package:dplyr)

## packages
library(ape)
library(phylofactor)
library(tidyverse)
library(data.table)
library(ggtree)
library(plyr)
library(caper)
library(patchwork)
library(ggpubr)

## load final ensemble
batin=read.csv("Cleaned Files/BatModels_IS.csv")
batout=read.csv("Cleaned Files/BatModels_OS.csv")

## combine
batnew=rbind.data.frame(batin[intersect(names(batin),names(batout))],
                        batout[intersect(names(batin),names(batout))])
batnew=batnew[c("Sp","Betacov","Rank","PropRank","InSample")]

## load weighted final ensemble
batnew2=read.csv("WeightedEnsemble2021.csv")
names(batnew2)=c("Sp","Betacov","PropRank","Ensemble")

## bind InSample
batnew2=merge(batnew2,batnew[c("Sp","InSample")],by="Sp")

## fix InSample
batnew2$InSample=ifelse(batnew2$Betacov==1 & batnew2$InSample==0,1,batnew2$InSample)

## fix dataset
bats=batnew2
rm(batnew,batnew2)

## load mammal supertree
setwd("~/Desktop/virionette")
tree = readRDS(paste0(getwd(), 
                      "/04_predictors/Full Supertree.rds"))

## extract node 6580 (Chiroptera)
btree = extract.clade(tree, 6580)
length(btree$tip.label)

## ladderize
btree=ladderize(btree)

## make tree names
bats$treenames = gsub(' ', '_', bats$Sp)

## load in taxonomy
setwd("~/Desktop/becker-betacov")
taxonomy = read.csv('mammal taxonomy.csv', 
                    header = T)
taxonomy$X = NULL
taxonomy$Sp = NULL

## trim to genus
taxonomy = taxonomy[!duplicated(taxonomy$hGenus), ]

## get genus
bats$hGenus = sapply(strsplit(bats$treenames, '_'), function(x) x[1])

## merge
bats = merge(bats, taxonomy, by = 'hGenus', all.x = T)

## fix bats
bats$hFamily = as.character(bats$hFamily)
bats$hFamily2 = revalue(bats$hGenus, 
                        c('Aproteles' = 'Pteropodidae', 
                          'Paracoelops' = 'Hipposideridae'))
bats$hFamily = ifelse(is.na(bats$hFamily), bats$hFamily2, bats$hFamily)
bats$hFamily2 = NULL

## make taxonomy 
bats$taxonomy = with(bats, paste(hOrder, hFamily, hGenus, treenames, sep = '; '))

## comparative frame
bats = bats[match(btree$tip.label, bats$treenames), ]
bats=bats[!is.na(bats$Sp),]
bats = comparative.data(phy = btree, data = bats, names.col = treenames, vcv = T, na.omit = F, warn.dropped = T)
bats$data$treenames = rownames(bats$data)
bats$data$Species = rownames(bats$data)

## phylfactor of ranks
temp=bats
set.seed(1)
temp_pf = gpf(Data = temp$data, 
              tree = temp$phy, 
              frmla.phylo = PropRank~phylo, 
              family = gaussian, algorithm = 'phylo', 
              nfactors = 9)

## results
temp_results = pfsum(temp_pf)

## split data from results
temp_data = temp_results$set
temp_results = temp_results$results

## visualize
gg = ggtree(temp$phy, 
            size = 0.15, 
            layout = 'circular') + 
  theme(legend.position = "none")

## add clades
for(i in 1:nrow(temp_results)){
  
  gg = gg+
    geom_hilight(node = temp_results$node[i], 
                 alpha = 0.5, 
                 fill = ifelse(temp_results$clade<
                                 temp_results$other, 
                               pcols[2], pcols[1])[i])+
    geom_cladelabel(node = temp_results$node[i], 
                    label = temp_results$factor[i], 
                    offset = 25, 
                    offset.text = 20)
}

## visualize
final = gg+
  
  ## title
  #ggtitle('bat rank, in-sample')+
  #theme(plot.title = element_text(hjust = 0.5))+
  
  ## add predictions
  geom_segment(data = segfun(temp, 25), 
               aes(x = x, y = y, xend = xend, yend = yend, colour = Betacov), size = 0.2)+
  scale_colour_manual(values = c('grey70', 'black'))#+

## export fig
setwd("~/Desktop/Fresnel_Jun")
png("Figures/Figure 5_phylo.png",width=4,height=4,units="in",res=300)
final
dev.off()

## export pdf
setwd("~/Desktop/Fresnel_Jun")
pdf("Figures/Figure 5_phylo.pdf",width=4,height=4)
final
dev.off()

## export table
write.csv(temp_results, 'final ensemble phylofactor results.csv')