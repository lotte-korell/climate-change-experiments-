
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



#####################
####### alpha #######
#####################

###############################
########## cum abs ######
###############################

###model selection 

m1.cum.abs<-lmer(LRR_cum.abs~MAP_rescale*delta.P1_rescale+
                   duration_rescale+
                   (1|study:site:SiteBlock), REML=F,data=Cum.abs_local, na.action="na.fail")

d1.cum.abs<-dredge(m1.cum.abs)#, m.lim=c(0,3))
print(d1.cum.abs)

#simples model delta P, delta AIC 2.51

d1.cum.absAIC<-model.avg(d1.cum.abs, subset = delta < 2)
summary(d1.cum.absAIC)

##best models < delta AIC 2

m1.cum.abs_best<-lmer(LRR_cum.abs~delta.P1_rescale*MAP_rescale+
                   (1|study:site:SiteBlock), REML=T,data=Cum.abs_local, na.action="na.fail")
summary(m1.cum.abs_best)
anova(m1.cum.abs_best)
Anova(m1.cum.abs_best)
r.squaredGLMM(m1.cum.abs_best)


m1.cum.abs_best1<-lmer(LRR_cum.abs~delta.P1_rescale+
                   (1|study:site:SiteBlock), REML=T,data=Cum.abs_local,na.action="na.fail")
summary(m1.cum.abs_best1)
anova(m1.cum.abs_best1)
Anova(m1.cum.abs_best1)
r.squaredGLMM(m1.cum.abs_best1)


###############################
####### Species richness ######
###############################


####### Precip ######

####model selection

m1.beta<-lmer(LRR_S~MAP_rescale*delta.P1_rescale+
              duration_rescale+
                (1|study:site:SiteBlock),REML=F,data=Div.data_local, na.action="na.fail")

##########---- multimodel inference----###########

d1.beta<-dredge(m1.beta)
print(d1.beta)

d1.AIC<-model.avg(d1.beta, subset = delta < 2)
summary(d1.AIC)

#best models

m1.best1<-lmer(LRR_S~delta.P1_rescale*MAP_rescale+
                 (1|study:site:SiteBlock), REML=T,data=Div.data_local, na.action="na.fail")
anova(m1.best1)
Anova(m1.best1)
r.squaredGLMM(m1.best1)
qqnorm(resid(m1.best1))
qqline(resid(m1.best1))
hist(resid(m1.best1))


m1.best2<-lmer(LRR_S~delta.P1_rescale*MAP_rescale+duration_rescale+
                 (1|study:site:SiteBlock),REML=T,data=Div.data_local, na.action="na.fail")
anova(m1.best2)
Anova(m1.best2)
r.squaredGLMM(m1.best2)
qqnorm(resid(m1.best2))
qqline(resid(m1.best2))
hist(resid(m1.best2))

m1.best3<-lmer(LRR_S~delta.P1_rescale+
                 (1|study:site:SiteBlock),REML=T,data=Div.data_local, na.action="na.fail")
anova(m1.best3)
Anova(m1.best3)
r.squaredGLMM(m1.best3)
qqnorm(resid(m1.best3))
qqline(resid(m1.best3))
hist(resid(m1.best3))



###############################
#######     ENS PIE     #######
##############################


m2.beta<-lmer(LRR_SPie~MAP_rescale*delta.P1_rescale+
                duration_rescale+
                (1|study:site:SiteBlock),REML=F,data=Div.data_local, na.action="na.fail")

##########---- multimodel inference----###########

d2.beta<-dredge(m2.beta)
print(d2.beta)

d2.AIC<-model.avg(d2.beta, subset = delta < 2)
summary(d2.AIC)

####best models 

m2.best<-lmer(LRR_SPie~delta.P1_rescale+
                (1|study:site:SiteBlock),REML=T,data=Div.data_local, na.action="na.fail")
anova(m2.best)
Anova(m2.best)
r.squaredGLMM(m2.best) 
qqnorm(resid(m2.best))
qqline(resid(m2.best))
hist(resid(m2.best))


m2.best2<-lmer(LRR_SPie~MAP_rescale+delta.P1_rescale+
                 (1|study:site:SiteBlock),REML=T,data=Div.data_local, na.action="na.fail")
anova(m2.best2)
Anova(m2.best2)
r.squaredGLMM(m2.best2)
qqnorm(resid(m2.best2))
qqline(resid(m2.best2))
hist(resid(m2.best2))


m2.best3<-lmer(LRR_SPie~delta.P1_rescale+duration_rescale+
                 (1|study:site:SiteBlock),REML=T,data=Div.data_local,na.action="na.fail")
anova(m2.best3)
Anova(m2.best3)
r.squaredGLMM(m2.best3)
qqnorm(resid(m2.best3))
qqline(resid(m2.best3))
hist(resid(m2.best3))



########################
######turnover #########
########################


################################
##### Species richness  ########
################################

m3a.beta<-lmer(LRR_betaSn~MAP_rescale*delta.P1_rescale+
                 duration_rescale+
                (1|study:site:SiteBlock),REML=F,data=Div.data_local, na.action="na.fail")


##########---- multimodel inference----###########

d3a.beta<-dredge(m3a.beta)
print(d3a.beta)

d3a.AIC<-model.avg(d3a.beta, subset = delta < 2)
summary(d3a.AIC)

m3a.best1<-lmer(LRR_betaSn~delta.P1_rescale+
                 (1|study:site:SiteBlock),REML=T,data=Div.data_local, na.action="na.fail")
anova(m3a.best1)
Anova(m3a.best1)
r.squaredGLMM(m3a.best1) 
qqnorm(resid(m3a.best1))
qqline(resid(m3a.best1))
hist(resid(m3a.best1))


m3a.best2<-lmer(LRR_betaSn~duration_rescale+
                 (1|study:site:SiteBlock),REML=T,data=Div.data_local, na.action="na.fail")
anova(m3a.best2)
Anova(m3a.best2)
r.squaredGLMM(m3a.best2) 
qqnorm(resid(m3a.best2))
qqline(resid(m3a.best2))
hist(resid(m3a.best2))



######################
####### ENS PIE#######
######################

m4.beta<-lmer(LRR_betaSPie~MAP_rescale*delta.P1_rescale+
                duration_rescale+
                (1|study:site:SiteBlock),REML=F,data=Div.data_local, na.action="na.fail")

##########model selection

d4.beta<-dredge(m4.beta)
print(d4.beta)

d4.AIC<-model.avg(d4.beta, subset = delta < 2)
summary(d4.AIC)


m4.best<-lmer(LRR_betaSPie~duration_rescale+*delta.P1_rescale
                +(1|study:site:SiteBlock),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best)
Anova(m4.best)
r.squaredGLMM(m4.best)
qqnorm(resid(m4.best))
qqline(resid(m4.best))
hist(resid(m4.best))


m4.best1<-lmer(LRR_betaSPie~delta.P1_rescale
                 +(1|study:site:SiteBlock),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best1)
Anova(m4.best1)
r.squaredGLMM(m4.best1) 
qqnorm(resid(m4.best1))
hist(resid(m4.best1))


m4.best3<-lmer(LRR_betaSPie~duration_rescale+delta.P1_rescale*MAP_rescale+
                 (1|study:site:SiteBlock),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best3)
Anova(m4.best3)
r.squaredGLMM(m4.best3) 
qqnorm(resid(m4.best3))
hist(resid(m4.best3))


m4.best4<-lmer(LRR_betaSPie~duration_rescale+delta.P1_rescale+MAP_rescale
                 +(1|study:site:SiteBlock),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best4)
Anova(m4.best4)
r.squaredGLMM(m4.best4) 
qqnorm(resid(m4.best4))
hist(resid(m4.best4))


m4.best5<-lmer(LRR_betaSPie~delta.P1_rescale+MAP_rescale+
                 (1|study:site:SiteBlock),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best5)
Anova(m4.best5)
r.squaredGLMM(m4.best5) 
qqnorm(resid(m4.best5))
hist(resid(m4.best5))

m4.best6<-lmer(LRR_betaSPie~delta.P1_rescale+MAP_rescale+
                 (1|study:site:SiteBlock),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best6)
Anova(m4.best6)
r.squaredGLMM(m4.best6) 
qqnorm(resid(m4.best6))
hist(resid(m4.best6))



####################
######gamma#########
####################

####################
##### cum abs ######
####################

m7.cum.abs<-lmer(LRR_cum.abs_gamma~MAP_rescale*delta.P1_rescale+
                 duration_rescale+
                   (1|study:site),REML=F,data=Cum.abs_gamma, na.action="na.fail")

d7.beta<-dredge(m7.cum.abs)
print(d7.beta)

d7.AIC<-model.avg(d7.beta, subset = delta < 2)
summary(d7.AIC)

###best model

m7.cum.abs1<-lmer(LRR_cum.abs_gamma~delta.P1_rescale+
                    (1|study:site), REML=T,data=Cum.abs_gamma, na.action="na.fail")
anova(m7.cum.abs1)
Anova(m7.cum.abs1)
r.squaredGLMM(m7.cum.abs1)


################################
##### Species richness  ########
################################


m7.beta<-lmer(LRR_gammaSn~MAP_rescale*delta.P1_rescale+
                duration_rescale + (1|study:site),  REML=F,data=Div.data_gamma, na.action="na.fail")


##########---- multimodel inference----###########

d7.beta<-dredge(m7.beta)
print(d7.beta)

d7.AIC<-model.avg(d7.beta, subset = delta < 2)
summary(d7.AIC)

##best models 

m7.best<-lmer(LRR_gammaSn~MAP_rescale*delta.P1_rescale
              + (1|study:site),  REML=F,data=Div.data_gamma, na.action="na.fail")

anova(m7.best)
Anova(m7.best)
r.squaredGLMM(m7.best) 
qqnorm(resid(m7.best))
hist(resid(m7.best))

m7.best1<-lmer(LRR_gammaSn~MAP_rescale*delta.P1_rescale+ duration_rescale
               + (1|study:site),  REML=F,data=Div.data_gamma, na.action="na.fail")

anova(m7.best1)
Anova(m7.best1)
r.squaredGLMM(m7.best1) 
qqnorm(resid(m7.best1))
hist(resid(m7.best1))

m7.best2<-lmer(LRR_gammaSn~delta.P1_rescale
               + (1|study:site),  REML=F,data=Div.data_gamma, na.action="na.fail")

anova(m7.best2)
Anova(m7.best2)
r.squaredGLMM(m7.best2) 
qqnorm(resid(m7.best2))
hist(resid(m7.best2))



####################
##### ENS PIE ######
####################


m6.beta<-lmer(LRR_gammaSPie~MAP_rescale*delta.P1_rescale+
                duration_rescale +
                (1|study:site),REML=F, data=Div.data_gamma, na.action="na.fail")
qqnorm(resid(m6.beta))
hist(resid(m6.beta))
r.squaredGLMM(m6.beta)

d6.beta<-dredge(m6.beta)
print(d6.beta)

d6.AIC<-model.avg(d6.beta, subset = delta < 2)
summary(d6.AIC)

###best models

m6.best<-lmer(LRR_gammaSPie~duration_rescale+
                (1|study:site),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m6.best)
Anova(m6.best)
r.squaredGLMM(m6.best)
qqnorm(resid(m6.best))
qqline(resid(m6.best))
hist(resid(m6.best))

m6.best2<-lmer(LRR_gammaSPie~delta.P1_rescale+duration_rescale+
                 (1|study:site),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m6.best2)
Anova(m6.best2)
r.squaredGLMM(m6.best2) 
qqnorm(resid(m6.best2))
hist(resid(m6.best2))

m6.best3<-lmer(LRR_gammaSPie~MAP_rescale+duration_rescale+
                 (1|study:site),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m6.best3)
Anova(m6.best3)
r.squaredGLMM(m6.best3) 
qqnorm(resid(m6.best3))
hist(resid(m6.best3))

m6.best4<-lmer(LRR_gammaSPie~delta.P1_rescale+MAP_rescale+duration_rescale+
                 (1|study:site),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m6.best4)
Anova(m6.best4)
r.squaredGLMM(m6.best4) 
qqnorm(resid(m6.best4))
hist(resid(m6.best4))

m6.best5<-lmer(LRR_gammaSPie~MAP_rescale+
                 (1|study:site),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m6.best5)
Anova(m6.best5)
r.squaredGLMM(m6.best5) 
qqnorm(resid(m6.best5))
hist(resid(m6.best5))


