ls()
rm(list=ls())
ls()

getwd()
setwd("C:/Users/korell/R/Meta.analysis/results")
dir()

library(lme4)
library(MASS)
library(car)
library(ggplot2)
library(MuMIn)
library(arm)
library(boot)

Div.data_local <- read.csv("Precip_data_final.csv" ,dec=".", sep=";",h=T)
Div.data_gamma <- read.csv("Precip_data_gamma_final.csv" ,dec=".", sep=";",h=T)
Cum.abs_local <- read.csv("ES.cum.abs_final.csv" ,dec=".", sep=";",h=T)
Cum.abs_gamma<- read.csv("ES.cum.abs.gamma_final.csv" ,dec=".", sep=";",h=T)


levels(as.factor(Div.data_local$ID.x))
Div.data_local$gr<-paste(Div.data_local$ID.x,Div.data_local$site)
unique(Div.data_local$gr)

site_mean<-aggregate(Div.data_local, by=list(Div.data_local$gr,Div.data_local$study,Div.data_local$delta.P1_rescale, Div.data_local$treatment.level), FUN=mean)
site_mean1<-subset(site_mean, select = c(Group.1, Group.2, Group.3, size.local, size.gamma, N.replicates, treat.duration.y, MAP_rescale, PET_rescale, LRR_S, LRR_SPie, LRR_betaSn, LRR_betaSPie))
colnames(site_mean1)<-c("study", "site", "delta.P.year", "size.local","size.gamma", "N.replicates", "duration", "MAP", "PET", "LRR_S", "LRR_SPie", "LRR_betaSn", "LRR_betaSPie")

site_var<-aggregate(Div.data_local, by=list(Div.data_local$gr,Div.data_local$study, Div.data_local$delta.P.year, Div.data_local$treatment.level), FUN=var)
site_var1<-subset(site_var, select = c(LRR_S, LRR_SPie, LRR_betaSn, LRR_betaSPie))
colnames(site_var1)<-c("var_LRR_S", "var_LRR_SPie","var_LRR_betaS", "var_LRR_betaSPie")

site_div<-cbind(site_mean1, site_var1)

##Correlation between plot size and effect sizes
# local 
cor.test(site_div$LRR_S,site_div$size.local)
cor.test(site_div$var_LRR_S,site_div$size.local)
cor.test(site_div$LRR_SPie,site_div$size.local)
cor.test(site_div$var_LRR_SPie,site_div$size.local)

# turnover
cor.test(site_div$LRR_betaS,site_div$size.local)
cor.test(site_div$var_LRR_betaS,site_div$size.local)
cor.test(site_div$LRR_betaSPie,site_div$size.local)
cor.test(site_div$var_LRR_betaSPie,site_div$size.local)


# site
cor.test(Div.data_gamma$LRR_gammaSn,Div.data_gamma$size.gamma)
cor.test(Div.data_gamma$LRR_gammaSPie,Div.data_gamma$size.gamma)


## confound between plot size, replictaes & MAP
cor.test(site_div$size.local,site_div$MAP)
cor.test(site_div$size.gamma,site_div$MAP)
cor.test(site_div$N.replicates,site_div$MAP)

## confound between plot size, replictaes & PET
cor.test(site_div$size.local,site_div$PET)
cor.test(site_div$size.gamma,site_div$PET)
cor.test(site_div$N.replicates,site_div$PET)

## confound between plot size, replictaes & delta P
cor.test(site_div$size.local,site_div$delta.P.year)
cor.test(site_div$size.gamma,site_div$delta.P.year)
cor.test(site_div$N.replicates,site_div$delta.P.year)

## confound between plot size, replictaes & duration
cor.test(site_div$size.local,site_div$duration)
cor.test(site_div$size.gamma,site_div$duration)
cor.test(site_div$N.replicates,site_div$duration)


