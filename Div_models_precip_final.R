
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

##convert duration into factor to include it as random effect in the model
Div.data_local$treat.duration.m<-as.factor(Div.data_local$treat.duration.m)
Div.data_gamma$treat.duration.m<-as.factor(Div.data_gamma$treat.duration.m)
Cum.abs_local$treat.duration.m<-as.factor(Cum.abs_local$treat.duration.m)
Cum.abs_gamma$treat.duration.m<-as.factor(Cum.abs_gamma$treat.duration.m)

########################
######grand mean########
########################

##total abundance at local, turnover and gamma scale 


grand_mean.cum.abs_P <- lmer(LRR_cum.abs ~ 1 + (1|study:site:SiteBlock:treat.duration.m), data=Cum.abs_local) 
summary(grand_mean.cum.abs_P)
r.squaredGLMM(grand_mean.cum.abs_P)

grand_mean.gamma.cum.abs_P <- lmer(LRR_cum.abs_gamma ~ 1 + (1|study:site:treat.duration.m), data=Cum.abs_local) 
summary(grand_mean.gamma.cum.abs_P)
r.squaredGLMM(grand_mean.gamma.cum.abs_P)

##species richness at local, turnover and gamma scale 

grand_mean.S_P <- lmer(LRR_S ~ 1 + (1|study:site:SiteBlock:treat.duration.m), data=Div.data_local) 
summary(grand_mean.S_P)
r.squaredGLMM(grand_mean.S_P)

grand_mean.betaS_P <- lmer(LRR_betaS ~ 1 + (1|study:site:SiteBlock:treat.duration.m),  data=Div.data_local) 
summary(grand_mean.betaS_P)
r.squaredGLMM(grand_mean.betaS_P)


grand_mean.gammaS_P <- lmer(LRR_gammaS ~ 1 + (1|study:site:treat.duration.m), data=Div.data_gamma) 
summary(grand_mean.gammaS_P)
r.squaredGLMM(grand_mean.gammaS_P)

##evenness at local, turnover and gamma scale 

grand_mean.SPie_P <- lmer(LRR_SPie ~ 1 + (1|study:site:SiteBlock:treat.duration.m),  data=Div.data_local)  
summary(grand_mean.SPie_P)
r.squaredGLMM(grand_mean.SPie_P)

grand_mean.betaSPie_P <- lmer(LRR_betaSPie ~ 1 + (1|study:site:SiteBlock:treat.duration.m),  data=Div.data_local) 
summary(grand_mean.betaSPie_P)
r.squaredGLMM(grand_mean.betaSPie_P)

grand_mean.gammaSPie_P <- lmer(LRR_gammaSPie ~ 1 + (1|study:site:treat.duration.m), data=Div.data_gamma) 
summary(grand_mean.gammaSPie_P)
r.squaredGLMM(grand_mean.gammaSPie_P)

#####################
####### alpha #######
#####################

###############################
########## cum abs ######
###############################

###model selection 

m1.cum.abs<-lmer(LRR_cum.abs~MAT_rescale*delta.P1_rescale+MAP_rescale*delta.P1_rescale+
                   (1|study:site:SiteBlock:treat.duration.m), REML=F,data=Cum.abs_local, na.action="na.fail")

d1.cum.abs<-dredge(m1.cum.abs)#, m.lim=c(0,3))
print(d1.cum.abs)

d1.cum.absAIC<-model.avg(d1.cum.abs, subset = delta < 2)
summary(d1.cum.absAIC)

##best models < delta AIC 2

m1.cum.abs<-lmer(LRR_cum.abs~delta.P1_rescale+
                   (1|study:site:SiteBlock:treat.duration.m), REML=T,data=Cum.abs_local, na.action="na.fail")
summary(m1.cum.abs)
anova(m1.cum.abs)
Anova(m1.cum.abs)
r.squaredGLMM(m1.cum.abs)


m2.cum.abs<-lmer(LRR_cum.abs~MAP_rescale+delta.P1_rescale+
                   (1|study:site:SiteBlock:treat.duration.m), REML=T,data=Cum.abs_local,na.action="na.fail")
summary(m2.cum.abs)
anova(m2.cum.abs)
Anova(m2.cum.abs)
r.squaredGLMM(m2.cum.abs)

m3.cum.abs<-lmer(LRR_cum.abs~MAT_rescale+delta.P1_rescale+
                   (1|study:site:SiteBlock:treat.duration.m), REML=T,data=Cum.abs_local,na.action="na.fail")
summary(m3.cum.abs)
anova(m3.cum.abs)
Anova(m3.cum.abs)
r.squaredGLMM(m3.cum.abs)

###bootstrapping


####best models for estimates

m1.cum<-lmer(LRR_cum.abs~delta.P1_rescale+
               (1|study:site:SiteBlock:treat.duration.m)-1,REML=T, data=Cum.abs_local,na.action="na.fail")

m1.cum.sum<-summary(m1.cum)
m1.cum.est<-m1.cum.sum$coefficients

set.seed(1234)
m1.cum_boot<-bootMer(m1.cum, FUN = fixef, nsim = 1000)
m1.cum_boot
m1.cum_ci<-boot.ci(m1.cum_boot, index =1, type = "perc")
m1.cum_ci


m2.cum<-lmer(LRR_cum.abs~delta.P1_rescale+MAP_rescale+
               (1|study:site:SiteBlock:treat.duration.m)-1,REML=T, data=Cum.abs_local, na.action="na.fail")

m2.cum.sum<-summary(m2.cum)
m2.cum.est<-m2.cum.sum$coefficients

###bootstrapping
#set.seed(1234)
#m2.cum_boot<-bootMer(m2.cum, FUN = fixef, nsim = 1000)
#m2.cum_ci1<-boot.ci(m2.cum_boot, index =1, type = "perc")# delta P
#m2.cum_ci2<-boot.ci(m2.cum_boot, index =2, type = "perc")# MAP

m3.cum<-lmer(LRR_cum.abs~delta.P1_rescale+MAP_rescale+
               (1|study:site:SiteBlock:treat.duration.m)-1,REML=T, data=Cum.abs_local, na.action="na.fail")

m3.cum.sum<-summary(m3.cum)
m3.cum.est<-m3.cum.sum$coefficients

m1.beta<-lmer(LRR_S~MAT_rescale*delta.P1_rescale+MAP_rescale*delta.P1_rescale+
                (1|study:site:SiteBlock:treat.duration.m), REML=F,data=Div.data_local,  na.action="na.fail")
plot(ranef(m1.beta))
d1.beta<-dredge(m1.beta)#, m.lim=c(0,3))
print(d1.beta)

###############################
####### Species richness ######
###############################


####### Precip ######

####multimodel average

m1.beta<-lmer(LRR_S~MAT_rescale*delta.P1_rescale+MAP_rescale*delta.P1_rescale+
                (1|study:site:SiteBlock:treat.duration.m),REML=F,data=Div.data_local, na.action="na.fail")

##########---- multimodel inference----###########

d1.beta<-dredge(m1.beta)
print(d1.beta)

d1.AIC<-model.avg(d1.beta, subset = delta < 2)
summary(d1.AIC)



m1.best1<-lmer(LRR_S~delta.P1_rescale*MAP_rescale+
                 (1|study:site:SiteBlock:treat.duration.m), REML=T,data=Div.data_local, na.action="na.fail")
anova(m1.best1)
Anova(m1.best1)
r.squaredGLMM(m1.best1)
qqnorm(resid(m1.best1))
qqline(resid(m1.best1))
hist(resid(m1.best1))


m1.best2<-lmer(LRR_S~delta.P1_rescale+
                 (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m1.best2)
Anova(m1.best2)
r.squaredGLMM(m1.best2)
qqnorm(resid(m1.best2))
qqline(resid(m1.best2))
hist(resid(m1.best2))



####for estimates

m1.mod<-lmer(LRR_S~MAP_rescale*delta.P1_rescale+
               (1|study:site:SiteBlock:treat.duration.m)-1,REML=T,data=Div.data_local,  na.action="na.fail")

m1.sum<-summary(m1.mod)
m1.est<-m1.sum$coefficients
m1.est

##bootstrapping effect size and CI
set.seed(1234)
m1.est_boot<-bootMer(m1.mod, FUN = fixef, nsim = 1000)
m1.est_boot#t1 = MAP, t2 = delta P, t3 = elta P * MAP
m1.est_ci1<-boot.ci(m1.est_boot, index =1, type = "perc")
m1.est_ci1#MAP
m1.est_ci2<-boot.ci(m1.est_boot, index =2, type = "perc")
m1.est_ci2#deltaP
m1.est_ci3<-boot.ci(m1.est_boot, index =3, type = "perc")
m1.est_ci3#MAP*delta P



m1.mod2<-lmer(LRR_S~delta.P1_rescale+
                (1|study:site:SiteBlock:treat.duration.m)-1,REML=T,data=Div.data_local, na.action="na.fail")

m1.sum2<-summary(m1.mod2)
m1.est2<-m1.sum2$coefficients
m1.est2

##bootstrapping effect size and CI
set.seed(1234)
m1.est2_boot<-bootMer(m1.mod2, FUN = fixef, nsim = 1000)
m1.est2_boot#effect size delta P
m1.est2_ci1<-boot.ci(m1.est2_boot, index =1, type = "perc")#deltaP
m1.est2_ci1#CI delta P


###############################
#######     ENS PIE     #######
##############################


m2.beta<-lmer(LRR_SPie~MAT_rescale*delta.P1_rescale+MAP_rescale*delta.P1_rescale+
                (1|study:site:SiteBlock:treat.duration.m),REML=F,data=Div.data_local, na.action="na.fail")

##########---- multimodel inference----###########

d2.beta<-dredge(m2.beta)
print(d2.beta)

d2.AIC<-model.avg(d2.beta, subset = delta < 2)
summary(d2.AIC)

####best models 

m2.best<-lmer(LRR_SPie~delta.P1_rescale+
                (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m2.best)
Anova(m2.best)
r.squaredGLMM(m2.best) 
qqnorm(resid(m2.best))
qqline(resid(m2.best))
hist(resid(m2.best))


m2.best2<-lmer(LRR_SPie~MAP_rescale*delta.P1_rescale+
                 (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m2.best2)
Anova(m2.best2)
r.squaredGLMM(m2.best2)
qqnorm(resid(m2.best2))
qqline(resid(m2.best2))
hist(resid(m2.best2))


m2.best3<-lmer(LRR_SPie~MAP_rescale*delta.P1_rescale+MAT_rescale+
                 (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local,na.action="na.fail")
anova(m2.best3)
Anova(m2.best3)
r.squaredGLMM(m2.best3)
qqnorm(resid(m2.best3))
qqline(resid(m2.best3))
hist(resid(m2.best3))


m2.best4<-lmer(LRR_SPie~MAP_rescale+delta.P1_rescale+
                 (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local,na.action="na.fail")
anova(m2.best4)
Anova(m2.best4)
r.squaredGLMM(m2.best4)
qqnorm(resid(m2.best4))
qqline(resid(m2.best4))
hist(resid(m2.best4))


m2.best5<-lmer(LRR_SPie~MAT_rescale+delta.P1_rescale+
                 (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local,na.action="na.fail")
anova(m2.best5)
Anova(m2.best5)
r.squaredGLMM(m2.best5)
qqnorm(resid(m2.best5))
qqline(resid(m2.best5))
hist(resid(m2.best5))


###estimates 

m2.mod<-lmer(LRR_SPie~delta.P1_rescale+
               (1|study:site:SiteBlock:treat.duration.m)-1,REML=T, data=Div.data_local,na.action="na.fail")

m2.sum<-summary(m2.mod)
m2.est<-m2.sum$coefficients

###bootstrapping
set.seed(1234)
m2.est_boot<-bootMer(m2.mod, FUN = fixef, nsim = 1000)
m2.est_boot #effect size delta P
m2.est_ci1<-boot.ci(m2.est_boot, index =1, type = "perc")#deltaP
m2.est_ci1 #CI delta P



m2.mod2<-lmer(LRR_SPie~delta.P1_rescale*MAP_rescale+
                (1|study:site:SiteBlock:treat.duration.m)-1,REML=T,data=Div.data_local, na.action="na.fail")

m2.sum2<-summary(m2.mod2)
m2.est2<-m2.sum2$coefficients

###bootstrapping
set.seed(1234)
m2.est2_boot<-bootMer(m2.mod2, FUN = fixef, nsim = 1000)
m2.est2_boot#t1 = MAP, t2 = delta P, t3 = elta P * MAP
m2.est2_ci1<-boot.ci(m2.est2_boot, index =1, type = "perc")#MAP
m2.est2_ci2#CI delta P
m2.est2_ci2<-boot.ci(m2.est2_boot, index =2, type = "perc")#deltaP
m2.est2_ci2#CI MAP
m2.est2_ci3<-boot.ci(m2.est2_boot, index =3, type = "perc")#deltaP*MAP
m2.est2_ci3#CI delta P*MAP



m2.mod3<-lmer(LRR_SPie~delta.P1_rescale*MAP_rescale+MAT_rescale+
                (1|study:site:SiteBlock:treat.duration.m)-1,REML=T,data=Div.data_local, na.action="na.fail")

m2.sum3<-summary(m2.mod3)
m2.est3<-m2.sum3$coefficients


###bootstapping
set.seed(1234)
m2.est3_boot<-bootMer(m2.mod3, FUN = fixef, nsim = 1000)
m2.est3_boot#t1 = delta P, t2 = MAP, t3 = MAT, t4 = delta P * MAP
m2.est3_ci1<-boot.ci(m2.est3_boot, index =1, type = "perc")
m2.est3_ci1#delta P
m2.est3_ci2<-boot.ci(m2.est3_boot, index =2, type = "perc")
m2.est3_ci2#MAP
m2.est3_ci3<-boot.ci(m2.est3_boot, index =3, type = "perc")
m2.est3_ci3#MAT
m2.est3_ci4<-boot.ci(m2.est3_boot, index =4, type = "perc")
m2.est3_ci4#deltaP*MAP

###Hier noch mehr Modelle einfügen






###################
######beta#########
##################

##############
##### S ######
#############


m3.beta<-lmer(LRR_betaS~MAT_rescale*delta.P1_rescale+MAP_rescale*delta.P1_rescale+
                (1|study:site:SiteBlock:treat.duration.m),REML=F,data=Div.data_local, na.action="na.fail")


##########---- multimodel inference----###########

d3.beta<-dredge(m3.beta)
print(d3.beta)

d3.AIC<-model.avg(d3.beta, subset = delta < 2)
summary(d3.AIC)


###best  models


m3.best1<-lmer(LRR_betaS~delta.P1_rescale+
                 (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m3.best1)
Anova(m3.best1)
r.squaredGLMM(m3.best1) 
qqnorm(resid(m3.best1))
qqline(resid(m3.best1))
hist(resid(m3.best1))


m3.best3<-lmer(LRR_betaS~delta.P1_rescale+MAP_rescale+
                 +(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m3.best3)
Anova(m3.best3)
r.squaredGLMM(m3.best3) 
qqnorm(resid(m3.best3))
qqline(resid(m3.best3))
hist(resid(m3.best3))


###estimates 


m3.mod1<-lmer(LRR_betaS~delta.P1_rescale+
                (1|study:site:SiteBlock:treat.duration.m)-1,REML=T,data=Div.data_local, na.action="na.fail")

m3.sum1<-summary(m3.mod1)
m3.est1<-m3.sum1$coefficients
m3.est1

###bootstrapping
set.seed(1234)
m3.est1_boot<-bootMer(m3.mod1, FUN = fixef, nsim = 1000)
m3.est1_boot#effect size delta P
m3.est1_ci1<-boot.ci(m3.est1_boot, index =1, type = "perc")#delta P
m3.est1_ci1#CI delta P


m3.mod3<-lmer(LRR_betaS~delta.P1_rescale+MAP_rescale+
                (1|study:site:SiteBlock:treat.duration.m)-1,REML=T,data=Div.data_local, na.action="na.fail")

m3.sum3<-summary(m3.mod3)
m3.est3<-m3.sum3$coefficients

###bootstrapping
set.seed(1234)
m3.est3_boot<-bootMer(m3.mod3, FUN = fixef, nsim = 1000)
m3.est3_ci1<-boot.ci(m3.est3_boot, index =1, type = "perc")#delta P
m3.est3_ci1<-boot.ci(m3.est3_boot, index =2, type = "perc")#MAP


######################
####### ENS PIE#######
######################

m4.beta<-lmer(LRR_betaSPie~MAT_rescale*delta.P1_rescale+MAP_rescale*delta.P1_rescale+
                +(1|study:site:SiteBlock:treat.duration.m),REML=F,data=Div.data_local, na.action="na.fail")


##########model selection

d4.beta<-dredge(m4.beta)
print(d4.beta)

d4.AIC<-model.avg(d4.beta, subset = delta < 2)
summary(d4.AIC)


m4.best<-lmer(LRR_betaSPie~MAP_rescale*delta.P1_rescale+
                +(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best)
Anova(m4.best)
r.squaredGLMM(m4.best)
qqnorm(resid(m4.best))
qqline(resid(m4.best))
hist(resid(m4.best))


m4.best1<-lmer(LRR_betaSPie~MAP_rescale*delta.P1_rescale+MAT_rescale+
                 +(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best1)
Anova(m4.best1)
r.squaredGLMM(m4.best1) 
qqnorm(resid(m4.best1))
hist(resid(m4.best1))

m4.best3<-lmer(LRR_betaSPie~delta.P1_rescale+
                 +(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best3)
Anova(m4.best3)
r.squaredGLMM(m4.best3) 
qqnorm(resid(m4.best3))
hist(resid(m4.best3))


###estimates 


m4.mod<-lmer(LRR_betaSPie~MAP_rescale*delta.P1_rescale+
                (1|study:site:SiteBlock:treat.duration.m)-1,REML=T,data=Div.data_local, na.action="na.fail")

m4.sum<-summary(m4.mod)
m4.est<-m4.sum$coefficients

##bootstrapping 
set.seed(1234)
m4.est_boot<-bootMer(m4.mod, FUN = fixef, nsim = 1000)
m4.est_boot# t1 = MAP, t2 = delta P, t3 = delta P * MAP
m4.est_ci1<-boot.ci(m4.est_boot, index =1, type = "perc")
m4.est_ci1#MAP
m4.est_ci2<-boot.ci(m4.est_boot, index =2, type = "perc")
m4.est_ci2#delta P
m4.est_ci3<-boot.ci(m4.est_boot, index =3, type = "perc")
m4.est_ci3#delta P * MAP



m4.mod2<-lmer(LRR_betaSPie~delta.P1_rescale+
                (1|study:site:SiteBlock:treat.duration.m)-1,REML=T,data=Div.data_local, na.action="na.fail")

m4.sum2<-summary(m4.mod2)
m4.est2<-m4.sum2$coefficients

##bootstrapping
set.seed(1234)
m4.est2_boot<-bootMer(m4.mod2, FUN = fixef, nsim = 1000)
m4.est2_boot#effect size delat P
m4.est2_ci1<-boot.ci(m4.est2_boot, index =1, type = "perc")
m4.est2_ci1#CI delta P


####################
######gamma#########
####################

####################
##### cum abs ######
####################

m7.cum.abs<-lmer(LRR_cum.abs_gamma~MAT_rescale*delta.P1_rescale+MAP_rescale*delta.P1_rescale+
                   +(1|study:site:treat.duration.m),REML=F,data=Cum.abs_gamma, na.action="na.fail")

d7.beta<-dredge(m7.cum.abs)
print(m7.cum.abs)

d7.AIC<-model.avg(d7.beta, subset = delta < 2)
summary(d7.AIC)

###best models 

m7.cum.abs1<-lmer(LRR_cum.abs_gamma~delta.P1_rescale+
                    (1|study:site:treat.duration.m), REML=T,data=Cum.abs_gamma, na.action="na.fail")
anova(m7.cum.abs1)
Anova(m7.cum.abs1)
r.squaredGLMM(m7.cum.abs1)

m7.cum.abs1<-lmer(LRR_cum.abs_gamma~delta.P1_rescale+MAT_rescale+
                    (1|study:site:treat.duration.m), REML=T,data=Cum.abs_gamma, na.action="na.fail")
anova(m7.cum.abs1)
Anova(m7.cum.abs1)
r.squaredGLMM(m7.cum.abs1)

m7.cum.abs2<-lmer(LRR_cum.abs_gamma~delta.P1_rescale+MAP_rescale+(1|study:site:treat.duration.m), REML=F,data=Cum.abs_gamma, na.action="na.fail")
anova(m7.cum.abs2)
Anova(m7.cum.abs2)
r.squaredGLMM(m7.cum.abs2)


###for estimates 

m7.mod<-lmer(LRR_cum.abs_gamma~delta.P1_rescale+
               (1|study:site:treat.duration.m)-1, REML=F,data=Cum.abs_gamma, na.action="na.fail")

m7.sum<-summary(m7.mod)
m7.est<-m7.sum$coefficients

set.seed(1234)
m7.est_boot<-bootMer(m7.mod, FUN = fixef, nsim = 1000)
m7.est_boot#effect size delta P
m7.est_ci1<-boot.ci(m7.est_boot, index =1, type = "perc")#delta P
m7.est_ci1#CI delta P

##############
##### S ######
#############


####Precip####

m5.beta<-lmer(LRR_gammaS~MAT_rescale*delta.P1_rescale+MAP_rescale*delta.P1_rescale+
                + (1|study:site:treat.duration.m),  REML=F,data=Div.data_gamma, na.action="na.fail")


##########---- multimodel inference----###########

d5.beta<-dredge(m5.beta)
print(d5.beta)

d5.AIC<-model.avg(d5.beta, subset = delta < 2)
summary(d5.AIC)

####only one model, no averaging possible!


m5.best<-lmer(LRR_gammaS~delta.P1_rescale+
                +(1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m5.best)
Anova(m5.best)
r.squaredGLMM(m5.best) 
qqnorm(resid(m5.best))
qqline(resid(m5.best))
hist(resid(m5.best))


m5.best2<-lmer(LRR_gammaS~MAP_rescale*delta.P1_rescale+
                 (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m5.best2)
Anova(m5.best2)
r.squaredGLMM(m5.best2) 
qqnorm(resid(m5.best2))
qqline(resid(m5.best2))
hist(resid(m5.best2))


m5.best3<-lmer(LRR_gammaS~delta.P1_rescale+MAT_rescale+
                 (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m5.best3)
Anova(m5.best3)
r.squaredGLMM(m5.best3) 
qqnorm(resid(m5.best3))
qqline(resid(m5.best3))
hist(resid(m5.best3))

m5.best4<-lmer(LRR_gammaS~MAT_rescale+
                 (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m5.best4)
Anova(m5.best4)
r.squaredGLMM(m5.best4) 
qqnorm(resid(m5.best4))
qqline(resid(m5.best4))
hist(resid(m5.best4))


m5.best5<-lmer(LRR_gammaS~delta.P1_rescale+MAP_rescale+
                 (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m5.best5)
Anova(m5.best5)
r.squaredGLMM(m5.best5) 
qqnorm(resid(m5.best5))
qqline(resid(m5.best5))
hist(resid(m5.best5))

###estimates 


m5.mod<-lmer(LRR_gammaS~delta.P1_rescale+
                (1|study:site:treat.duration.m)-1,REML=T,data=Div.data_gamma, na.action="na.fail")

m5.sum<-summary(m5.mod)
m5.est<-m5.sum$coefficients

##bootstrapping 
set.seed(1234)
m5.est_boot<-bootMer(m5.mod, FUN = fixef, nsim = 1000)
m5.est_boot#delta P
m5.est_ci1<-boot.ci(m5.est_boot, index =1, type = "perc")
m5.est_ci1#CI delta P


m5.mod1<-lmer(LRR_gammaS~delta.P1_rescale*MAP_rescale+
               (1|study:site:treat.duration.m)-1,REML=T,data=Div.data_gamma,  na.action="na.fail")

m5.sum1<-summary(m5.mod1)
m5.est1<-m5.sum1$coefficients

##bootstrapping 
set.seed(1234)
m5.est1_boot<-bootMer(m5.mod1, FUN = fixef, nsim = 1000)
m5.est1_boot#t1= delta P, t2 = MAP, t3 = delta P*MAP
m5.est1_ci1<-boot.ci(m5.est1_boot, index =1, type = "perc")
m5.est1_ci1# delta P
m5.est1_ci2<-boot.ci(m5.est1_boot, index =2, type = "perc")
m5.est1_ci2#MAP
m5.est1_ci3<-boot.ci(m5.est1_boot, index =3, type = "perc")
m5.est1_ci3#MAP*delta P



####################
##### ENS PIE ######
####################


####Precip####

m6.beta<-lmer(LRR_gammaSPie~MAT_rescale*delta.P1_rescale+MAP_rescale*delta.P1_rescale+
                (1|study:site:treat.duration.m),REML=F, data=Div.data_gamma, na.action="na.fail")
qqnorm(resid(m6.beta))
hist(resid(m6.beta))

##########---- multimodel inference----###########

d6.beta<-dredge(m6.beta)
print(d6.beta)

d6.AIC<-model.avg(d6.beta, subset = delta < 2)
summary(d6.AIC)

###best models

m6.best<-lmer(LRR_gammaSPie~MAP_rescale+
                (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m6.best)
Anova(m6.best)
r.squaredGLMM(m6.best)
qqnorm(resid(m6.best))
qqline(resid(m6.best))
hist(resid(m6.best))

m6.best2<-lmer(LRR_gammaSPie~MAT_rescale+
                 (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m6.best2)
Anova(m6.best2)
r.squaredGLMM(m6.best2) 
qqnorm(resid(m6.best2))
hist(resid(m6.best2))

m6.best3<-lmer(LRR_gammaSPie~delta.P1_rescale+
                 (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m6.best3)
Anova(m6.best3)
r.squaredGLMM(m6.best3) 
qqnorm(resid(m6.best3))
hist(resid(m6.best3))


###estimates 


m6.mod3<-lmer(LRR_gammaSPie~delta.P1_rescale+
                (1|study:site:treat.duration.m)-1,REML=T,data=Div.data_gamma, na.action="na.fail")

m6.sum3<-summary(m6.mod3)
m6.est3<-m6.sum3$coefficients


##bootstrapping effect size and CI
set.seed(1234)
m6.est3_boot<-bootMer(m6.mod3, FUN = fixef, nsim = 1000)
m6.est3_boot#delta P
m6.est3_ci1<-boot.ci(m6.est3_boot, index =1, type = "perc")
m6.est3_ci1#delta P

