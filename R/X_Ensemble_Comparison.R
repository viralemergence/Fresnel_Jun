## ensemble comparison
## danbeck@ou.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages for Figure 4
library(plyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(ggpubr)

## load 2020 ranks
setwd("~/Desktop/Fresnel/Cleaned Files")
batin=read.csv("BatModels_IS.csv")
batout=read.csv("BatModels_OS.csv")

## restrict BatModels_OS to out-of-sample
batout=batout[batout$InSample==0,]

## combine
batold=rbind.data.frame(batin[intersect(names(batin),names(batout))],
                        batout[intersect(names(batin),names(batout))])
batold=batold[c("Sp","Betacov","PropRank","InSample")]

## merge with ensemble threshold
setwd("~/Desktop/Fresnel")
batpred=read.csv("BinaryPredictions.csv")

## get model names
mnames=names(batpred)
mnames=mnames[!mnames%in%c("Sp","Ensemble","Betacov")]

## suspects host per model
shost=rep(NA,length(mnames))
names(shost)=mnames
for(i in mnames){
  
  ## subset batpred
  set=batpred[c("Betacov",i)]
  
  ## set model
  set$model=set[i]
  
  ## ifelse
  set$suspect=ifelse(set$model==T & set$Betacov==0,1,0)
  
  ## return
  ns=as.numeric(table(set$suspect)["1"])
  shost[i]=ns
  
}
range(shost)

## trim to main columns
batpred=batpred[c("Sp","Ensemble")]

## merge and clean
batold=merge(batold,batpred,by="Sp")
rm(batpred)

## make category
batold$category=ifelse(batold$Betacov==0 & batold$Ensemble==F,"unlikely",
                       ifelse(batold$Betacov==0 & batold$Ensemble==T,"suspected",
                              ifelse(batold$Ensemble==F & batold$Betacov==1,"false -","true +")))

## category 2
batold$category2=ifelse(batold$category%in%c("false -","true +"),"known host",batold$category)
table(batold$category2)

## tabulate
table(batold$category,batold$InSample)

## save names
nm=names(batold)

## fix names
names(batold)[2:ncol(batold)]=paste(names(batold)[2:ncol(batold)],"_V1",sep="")

## load new ranks
setwd("~/Desktop/Fresnel_Jun/Cleaned Files")
batin=read.csv("BatModels_IS.csv")
batout=read.csv("BatModels_OS.csv")

## combine
batnew=rbind.data.frame(batin[intersect(names(batin),names(batout))],
                        batout[intersect(names(batin),names(batout))])
batnew=batnew[c("Sp","Betacov","Rank","PropRank","InSample")]
table(batnew$Betacov)

## merge with ensemble threshold
setwd("~/Desktop/Fresnel_Jun")
batpred=read.csv("BinaryPredictions.csv")
batpred=batpred[c("Sp","Ensemble")]

## merge and clean
batnew=merge(batnew,batpred,by="Sp")
rm(batpred)

## make category
batnew$category=ifelse(batnew$Betacov==0 & batnew$Ensemble==F,"unlikely",
                       ifelse(batnew$Betacov==0 & batnew$Ensemble==T,"suspected",
                              ifelse(batnew$Ensemble==F & batnew$Betacov==1,"false -","true +")))

## category 2
batnew$category2=ifelse(batnew$category%in%c("false -","true +"),"known host",batnew$category)

## trim columns
batnew=batnew[nm]

## fix names
names(batnew)[2:ncol(batnew)]=paste(names(batnew)[2:ncol(batnew)],"_V2",sep="")

## load updated ensemble
setwd("~/Desktop/Fresnel_Jun")
batnew2=read.csv("WeightedEnsemble2021.csv")
names(batnew2)=c("Sp","Betacov","PropRank","Ensemble")

## bind InSample
batnew2=merge(batnew2,batnew[c("Sp","InSample_V2")],by="Sp")
batnew2$InSample=batnew2$InSample_V2
batnew2$InSample_V2=NULL

## fix InSample
batnew2$InSample=ifelse(batnew2$Betacov==1 & batnew2$InSample==0,1,batnew2$InSample)

## make category
batnew2$category=ifelse(batnew2$Betacov==0 & batnew2$Ensemble==F,"unlikely",
                        ifelse(batnew2$Betacov==0 & batnew2$Ensemble==T,"suspected",
                               ifelse(batnew2$Ensemble==F & batnew2$Betacov==1,"false -","true +")))

## category 2
batnew2$category2=ifelse(batnew2$category%in%c("false -","true +"),"known host",batnew2$category)

## trim columns
batnew2=batnew2[nm]

## fix names
names(batnew2)[2:ncol(batnew2)]=paste(names(batnew2)[2:ncol(batnew2)],"_V3",sep="")

## combine
bats=merge(batold,batnew,by="Sp")
bats=merge(bats,batnew2,by="Sp")

## tally known
table(bats$Betacov_V1)
table(bats$Betacov_V2)
table(bats$Betacov_V3)

## tally suspect
table(bats$category2_V1)
table(bats$category2_V2)
table(bats$category2_V3)

## designate hosts lost V1-V2
bats$lost1=ifelse(bats$category2_V1=="suspected" & bats$category2_V2=="unlikely","yes","no")

## designate hosts gained V1-V2
bats$new1=ifelse(bats$category2_V1=="unlikely" & bats$category2_V2=="suspected","yes","no")

## tabulate
table(bats$lost1)
table(bats$new1)

## correlate propranks
cor(bats$PropRank_V1,bats$PropRank_V2,method="spearman")

## designate hosts lost and gained in weighted ensemble
bats$lost2=ifelse(bats$category2_V1=="suspected" & bats$category2_V3=="unlikely","yes","no")
bats$new2=ifelse(bats$category2_V1=="unlikely" & bats$category2_V3=="suspected","yes","no")

## tabulate
table(bats$lost2)
table(bats$new2)

## correlate predictions
cor(bats$PropRank_V1,bats$PropRank_V3,method="spearman")

## designate new category for V1 vs V2
bats$revised=ifelse(bats$Betacov_V1!=bats$Betacov_V2,"new host",
                    ifelse(bats$category2_V1=="suspected" & bats$category2_V2=="suspected","retained suspect",
                           ifelse(bats$category2_V1=="suspected" & bats$category2_V2=="unlikely","lost",
                                  ifelse(bats$category2_V1=="unlikely" & bats$category2_V2=="suspected",
                                         "gained","unlikely host"))))

## order revised
bats$revised=factor(bats$revised,
                    levels=c("unlikely host",
                             "retained suspect",
                             "new host",
                             "lost",
                             "gained"))
table(bats$revised)

## repeat for V1 vs V3
bats$revised2=ifelse(bats$Betacov_V1!=bats$Betacov_V3,"new host",
                     ifelse(bats$category2_V1=="suspected" & bats$category2_V3=="suspected","retained suspect",
                            ifelse(bats$category2_V1=="suspected" & bats$category2_V3=="unlikely","lost",
                                   ifelse(bats$category2_V1=="unlikely" & bats$category2_V3=="suspected",
                                          "gained","unlikely host"))))
bats$revised2=factor(bats$revised2,levels=levels(bats$revised))

## matching colors
cs=c("grey80",
     brewer.pal(n=4,name="Pastel1"))

## sizes
sz=c(2,2,3,5,5)

## alpha
as=c(0.25,0.25,0.75,0.75,0.75)

## function to take dataset and subset to top 10 (lowest proprank)
subfun=function(V,type){
  
  ## get proprank and category for comparison ensemble
  set=bats
  
  ## set PropRank and category
  set$PropRank=set[paste("PropRank_V",V,sep="")][,1]
  set$category=set[paste("category_V",V,sep="")][,1]
  
  ## select revised
  if(V==2){
    
    set$change=set$revised
  }else{
    set$change=set$revised2
  }
  
  ## clean revised
  set$revised=NULL
  set$revised2=NULL
  
  ## trim to type
  set2=set[set$change==type,]
  
  ## rank by proprank
  set2=set2[order(set2$PropRank),]
  
  ## top 10
  set2=set2[1:10,]
  
  ## return
  return(set2)
}

## 4A
set.seed(1)
f4A=ggplot(bats,aes(PropRank_V1,PropRank_V2))+
  
  ## add 1:1 line
  geom_abline(intercept=0,slope=1,linetype=2,alpha=0.5,size=0.5)+
  
  ## add ensemble 1 threshold
  geom_vline(xintercept=1-0.5286,linetype=2,size=0.25,alpha=0.5)+
  
  ## add ensemble 2 threshold
  geom_hline(yintercept=1-0.5684,linetype=2,size=0.25,alpha=0.5)+
  
  ## add smooth
  geom_line(stat="smooth",method="lm",
            formula=y~x,
            size=1,
            alpha=0.5)+
  
  ## add point
  geom_point(aes(colour=revised,
                 size=revised,
                 alpha=revised))+
  
  ## colour
  scale_colour_manual(values=cs)+
  
  ## size
  scale_size_manual(values=sz)+
  
  ## alpha
  scale_alpha_manual(values=as)+
  
  ## guides
  guides(size=F,alpha=F,
         colour=guide_legend(title="ensemble status",
                             override.aes=list(alpha=1,
                                               size=sz)))+
  
  ## add labels for gain
  geom_text_repel(data=subfun(V=2,type="gained"),
                  aes(label=Sp),
                  direction="y",
                  hjust="left",
                  nudge_x=0.15,
                  force=2,
                  size=3,
                  nudge_y=-0.02,
                  segment.colour=cs[5],
                  segment.alpha=0.5,
                  fontface="italic")+
  
  ## add labels for loss
  geom_text_repel(data=subfun(V=2,type="lost"),
                  aes(label=Sp),
                  direction="y",
                  hjust="left",
                  nudge_x=-0.37,
                  force=2,
                  size=3,
                  nudge_y=0.33,
                  segment.colour=cs[4],
                  segment.alpha=0.5,
                  fontface="italic")+
  
  ## theme
  theme_bw()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12),
        strip.text=element_text(size=10),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12),
        legend.position="top")+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="original proportional rank",
       y="revised proportional rank")

## f4C
set.seed(1)
f4C=ggplot(bats,aes(PropRank_V1,PropRank_V3))+
  
  ## add 1:1 line
  geom_abline(intercept=0,slope=1,linetype=2,alpha=0.5,size=0.5)+
  
  ## add ensemble 1 threshold
  geom_vline(xintercept=1-0.5286,linetype=2,size=0.25,alpha=0.5)+
  
  ## add ensemble 2* threshold
  geom_hline(yintercept=1-0.4886,linetype=2,size=0.25,alpha=0.5)+
  
  ## add smooth
  geom_line(stat="smooth",method="lm",
            formula=y~x,
            size=1,
            alpha=0.5)+
  
  ## add point
  geom_point(aes(colour=revised2,
                 size=revised2,
                 alpha=revised2))+
  
  ## colour
  scale_colour_manual(values=cs)+
  
  ## size
  scale_size_manual(values=sz)+
  
  ## alpha
  scale_alpha_manual(values=as)+
  
  ## guides
  guides(size=F,alpha=F,
         colour=guide_legend(title="ensemble status",
                             override.aes=list(alpha=1,
                                               size=sz)))+
  
  ## add labels for gain
  geom_text_repel(data=subfun(V=3,type="gained"),
                  aes(label=Sp),
                  direction="y",
                  hjust="left",
                  nudge_x=0.25,
                  force=2,
                  size=3,
                  nudge_y=0.035,
                  segment.colour=cs[5],
                  segment.alpha=0.5,
                  fontface="italic")+
  
  ## add labels for loss
  geom_text_repel(data=subfun(V=3,type="lost"),
                  aes(label=Sp),
                  direction="y",
                  hjust="right",
                  nudge_x=-0.1,
                  force=2,
                  size=3,
                  nudge_y=0.4,
                  segment.colour=cs[4],
                  segment.alpha=0.5,
                  fontface="italic")+
  
  ## theme
  theme_bw()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12),
        strip.text=element_text(size=10),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12),
        legend.position="top")+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="original proportional rank",
       y="weighted revised proportional rank")

## number of top hosts
th=10 ## per in and out of sample

## plot top suspect hosts from V1
set1=bats[bats$category_V1=="suspected",]
set1=set1[order(set1$PropRank_V1),]

## extract in
in1=set1[set1$InSample_V1==1,][1:th,]

## extract out
out1=set1[set1$InSample_V1==0,][1:th,]

## combine
set1=rbind.data.frame(in1,out1)

## assign y
set1$idy=rep(th:1,2)

## fix in/out
set1$InSample_V1=ifelse(set1$InSample_V1==1,1.5,0)

## add astericks
set1$Sp2=ifelse(set1$revised=="new host",paste(set1$Sp,"*",sep=""),set1$Sp)

## plot
f4B=ggplot(set1,aes(InSample_V1,idy))+
  geom_text(aes(label=Sp2),fontface="italic")+
  scale_x_reverse(breaks=sort(unique(set1$InSample_V1)),
                  limits=rev(c(-0.5,2.25)),
                  labels=rev(c("in-sample","out-of-sample")))+
  ylim(c(0.5,10.5))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        panel.border=element_blank(),
        axis.text.x=element_text(size=12,colour="black"))

## plot top suspect hosts from V3
set2=bats[bats$category_V3=="suspected",]
set2=set2[order(set2$PropRank_V3),]

## extract in
in2=set2[set2$InSample_V3==1,][1:th,]

## extract out
out2=set2[set2$InSample_V3==0,][1:th,]

## combine
set2=rbind.data.frame(in2,out2)

## assign y
set2$idy=rep(th:1,2)

## fix in/out
set2$InSample_V3=ifelse(set2$InSample_V3==1,1.5,0)

## plot
f4D=ggplot(set2,aes(InSample_V3,idy))+
  geom_text(aes(label=Sp),fontface="italic")+
  scale_x_reverse(breaks=sort(unique(set2$InSample_V3)),
                  limits=rev(c(-0.5,2.25)),
                  labels=rev(c("in-sample","out-of-sample")))+
  ylim(c(0.5,10.5))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        panel.border=element_blank(),
        axis.text.x=element_text(size=12,colour="black"))

## plot
setwd("~/Desktop/Fresnel_Jun/Figures")
png("Figure 4.png",width=10,height=8,units="in",res=300)
ggarrange(f4A,f4C,f4B,f4D,common.legend = T,
          labels=c("(A)","(B)","(C)","(D)"),
          font.label=list(size=11.5,face="plain"),
          heights=c(1,0.65),
          label.x=-0.01,
          la0el.y=0.98)
dev.off()

## as PDF
setwd("~/Desktop/Fresnel_Jun/Figures")
pdf("Figure 4.pdf",width=10,height=8)
ggarrange(f4A,f4C,f4B,f4D,common.legend = T,
          labels=c("(A)","(B)","(C)","(D)"),
          font.label=list(size=11.5,face="plain"),
          heights=c(1,0.65),
          label.x=-0.01,
          label.y=0.98)
dev.off()