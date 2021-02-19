ls()
rm(list=ls())
ls()

getwd()
setwd("C:/Users/korell/R/Meta.analysis/results")
dir()

library(lme4)
library(car)
library(boot)
library(dplyr)

life_hist <- read.csv("Life_history.csv" ,dec=".", sep=";",h=T)
life_hist$study_site<-paste(life_hist$study, "-", life_hist$site)
Div.data_local <- read.csv("Div.data_local.csv" ,dec=".", sep=";",h=T)
Div.data_local$study_site<-paste(Div.data_local$study, "-", Div.data_local$site)

agr<-aggregate(Div.data_local, by=list(Div.data_local$study_site,Div.data_local$study),FUN=mean)#get mean plot scale species richness/ evenness
life_hist1<-merge(life_hist, agr[,-(3:23)], by.x="study_site",by.y="Group.1")
life_hist2<-distinct(life_hist1,site,.keep_all = TRUE)#only keep each datapoint once
life_hist3<-subset(life_hist2,!is.na(polycarp), !is.na(monocarp))# remove datasets that do nort probide species identity 

mod1<-glmer(monocarp ~ MAP_rescale + (1|study), family = binomial, data=life_hist3)
summary(mod1)
Anova(mod1)

mod2<-glmer(monocarp ~ PET_rescale + (1|study), family = binomial, data=life_hist3)
summary(mod2)
Anova(mod2)

mod4<-lmer(LRR_S~gr_mono+ (1|study), REML=T,data=life_hist3, na.action="na.fail")
summary(mod4)
Anova(mod4)

mod5<-lmer(LRR_SPie~gr_mono+ (1|study), REML=T,data=life_hist3, na.action="na.fail")
summary(mod5)
Anova(mod5)

mod6<-lmer(LRR_gammaSn~gr_mono+ (1|study), REML=T,data=life_hist3, na.action="na.fail")
summary(mod6)
Anova(mod6)

mod7<-lmer(LRR_gammaSPie~gr_mono+ (1|study), REML=T,data=life_hist3, na.action="na.fail")
summary(mod7)
Anova(mod7)

