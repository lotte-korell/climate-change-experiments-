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

levels(as.factor(Div.data_local$study))
levels(as.factor(Div.data_local$site))
levels(as.factor(Div.data_local$site))


########################
######grand mean########
########################

## Total abundance at local and gamma scale 

grand_mean.cum.abs_P <- lmer(LRR_cum.abs ~ 1 + (1|study:site:SiteBlock), data=Cum.abs_local) 
summary(grand_mean.cum.abs_P)
r.squaredGLMM(grand_mean.cum.abs_P)

set.seed(1234)
m1.grand_mean_boot<-bootMer(grand_mean.cum.abs_P, FUN = fixef, nsim = 1000)
m1.grand_mean_boot # estimate

grand_mean.gamma.cum.abs_P <- lmer(LRR_cum.abs_gamma ~ 1 + (1|study:site), data=Cum.abs_gamma) 
summary(grand_mean.gamma.cum.abs_P)
r.squaredGLMM(grand_mean.gamma.cum.abs_P)

set.seed(1234)
m7.grand_mean_boot<-bootMer(grand_mean.gamma.cum.abs_P, FUN = fixef, nsim = 1000)
m7.grand_mean_boot # estimate

## Species richness at local, turnover and gamma scale 

grand_mean.S_P <- lmer(LRR_Sn ~ 1 + (1|study:site:SiteBlock), data=Div.data_local) 
summary(grand_mean.S_P)
r.squaredGLMM(grand_mean.S_P)

set.seed(1234)
m1S.grand_mean_boot<-bootMer(grand_mean.S_P, FUN = fixef, nsim = 1000)
m1S.grand_mean_boot # estimate


grand_mean.betaS_P <- lmer(LRR_betaSn ~ 1 + (1|study:site:SiteBlock),  data=Div.data_local) 
summary(grand_mean.betaS_P)
r.squaredGLMM(grand_mean.betaS_P)

set.seed(1234)
m3.grand_mean_boot<-bootMer(grand_mean.betaS_P, FUN = fixef, nsim = 1000)
m3.grand_mean_boot # estimate


grand_mean.gammaSn_P <- lmer(LRR_gammaSn ~ 1 + (1|study:site), data=Div.data_gamma) 
summary(grand_mean.gammaSn_P)
r.squaredGLMM(grand_mean.gammaSn_P)

set.seed(1234)
m5.grand_mean_boot<-bootMer(grand_mean.gammaSn_P, FUN = fixef, nsim = 1000)
m5.grand_mean_boot # estimate

## Evenness at local, turnover and gamma scale 

grand_mean.SPie_P <- lmer(LRR_SPie ~ 1 + (1|study:site:SiteBlock),  data=Div.data_local)  
summary(grand_mean.SPie_P)
r.squaredGLMM(grand_mean.SPie_P)

set.seed(1234)
m2.grand_mean_boot<-bootMer(grand_mean.SPie_P, FUN = fixef, nsim = 1000)
m2.grand_mean_boot # estimate


grand_mean.betaSPie_P <- lmer(LRR_betaSPie ~ 1 + (1|study:site:SiteBlock),  data=Div.data_local) 
summary(grand_mean.betaSPie_P)
r.squaredGLMM(grand_mean.betaSPie_P)

set.seed(1234)
m4.grand_mean_boot<-bootMer(grand_mean.betaSPie_P, FUN = fixef, nsim = 1000)
m4.grand_mean_boot # estimate


grand_mean.gammaSPie_P <- lmer(LRR_gammaSPie ~ 1 + (1|study:site), data=Div.data_gamma) 
summary(grand_mean.gammaSPie_P)
r.squaredGLMM(grand_mean.gammaSPie_P)

set.seed(1234)
m6.grand_mean_boot<-bootMer(grand_mean.gammaSPie_P, FUN = fixef, nsim = 1000)
m6.grand_mean_boot # estimate

#######################################################################################
#################################  delta P  ###########################################
#######################################################################################

##############################~~~~~~~total abundance~~~~~~~~###########################
#local
m1.cum<-lmer(LRR_cum.abs~delta.P1_rescale+
               (1|study:site:SiteBlock)-1,REML=T, data=Cum.abs_local,na.action="na.fail")

m1.cum.sum<-summary(m1.cum)
m1.cum.est<-m1.cum.sum$coefficients

set.seed(1234)
m1.cum_boot<-bootMer(m1.cum, FUN = fixef, nsim = 1000)
m1.cum_boot # estimate
m1.cum_boot_ci1<-boot.ci(m1.cum_boot, index =1, type = "perc")
m1.cum_boot_ci1 #CI

#site
m7.mod<-lmer(LRR_cum.abs_gamma~delta.P1_rescale+
               (1|study:site)-1, REML=F,data=Cum.abs_gamma, na.action="na.fail")

m7.sum<-summary(m7.mod)
m7.est<-m7.sum$coefficients

set.seed(1234)
m7.est_boot<-bootMer(m7.mod, FUN = fixef, nsim = 1000)
m7.est_boot#effect size delta P
m7.est_ci1<-boot.ci(m7.est_boot, index =1, type = "perc")#delta P
m7.est_ci1#CI delta P


##############################~~~~~~~species richness~~~~~~~~###########################

#local 
m1.mod<-lmer(LRR_S~delta.P1_rescale+
               (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m1.sum<-summary(m1.mod)
m1.est<-m1.sum$coefficients
m1.est

##bootstrapping effect size and CI
set.seed(1234)
m1.est_boot<-bootMer(m1.mod, FUN = fixef, nsim = 1000)
m1.est_boot#estimate
m1.est_ci1<-boot.ci(m1.est_boot, index =1, type = "perc")#deltaP
m1.est_ci1#CI 

#turnover 
m3a.mod<-lmer(LRR_betaSn~delta.P1_rescale+
               (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m3a.sum<-summary(m3a.mod)
m3a.est<-m3a.sum$coefficients
m3a.est

###bootstrapping

set.seed(1234)
m3a.est_boot<-bootMer(m3a.mod, FUN = fixef, nsim = 1000)
m3a.est_boot#effect size delta P
m3a.est_ci1<-boot.ci(m3a.est_boot, index =1, type = "perc")
m3a.est_ci1#CI #delta P



# site 
m7.mod<-lmer(LRR_gammaSn~delta.P1_rescale+
               (1|study:site)-1,REML=T,data=Div.data_gamma, na.action="na.fail")

m7.sum<-summary(m7.mod)
m7.est<-m7.sum$coefficients


##bootstrapping effect size and CI
set.seed(1234)
m7.est_boot<-bootMer(m7.mod, FUN = fixef, nsim = 1000)
m7.est_boot#delta P
m7.est_ci1<-boot.ci(m7.est_boot, index =1, type = "perc")
m7.est_ci1#



##############################~~~~~~~~~~evenness~~~~~~~~################################
#local
m2.mod<-lmer(LRR_SPie~delta.P1_rescale+
               (1|study:site:SiteBlock)-1,REML=T, data=Div.data_local,na.action="na.fail")

m2.sum<-summary(m2.mod)
m2.est<-m2.sum$coefficients

###bootstrapping
set.seed(1234)
m2.est_boot<-bootMer(m2.mod, FUN = fixef, nsim = 1000)
m2.est_boot #effect size
m2.est_ci1<-boot.ci(m2.est_boot, index =1, type = "perc")#deltaP
m2.est_ci1 #CI 

#turnover
m4.mod3<-lmer(LRR_betaSPie~delta.P1_rescale+
                (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m4.sum3<-summary(m4.mod3)
m4.est3<-m4.sum3$coefficients

##bootstrapping
set.seed(1234)
m4.est3_boot<-bootMer(m4.mod3, FUN = fixef, nsim = 1000)
m4.est3_boot#effect size delta P
m4.est3_ci1<-boot.ci(m4.est3_boot, index =1, type = "perc")
m4.est3_ci1#CI delta P

#site
m6.mod<-lmer(LRR_gammaSPie~delta.P1_rescale+
               (1|study:site)-1,REML=T,data=Div.data_gamma, na.action="na.fail")

m6.sum<-summary(m6.mod)
m6.est<-m6.sum$coefficients


##bootstrapping effect size and CI
set.seed(1234)
m6.est_boot<-bootMer(m6.mod, FUN = fixef, nsim = 1000)
m6.est_boot#delta P
m6.est_ci1<-boot.ci(m6.est_boot, index =1, type = "perc")
m6.est_ci1#


#######################################################################################
#################################  MAP  ###########################################
#######################################################################################

##############################~~~~~~~total abundance~~~~~~~~###########################
#local
m1.cum<-lmer(LRR_cum.abs~MAP_rescale+
               (1|study:site:SiteBlock)-1,REML=T, data=Cum.abs_local,na.action="na.fail")

m1.cum.sum<-summary(m1.cum)
m1.cum.est<-m1.cum.sum$coefficients

set.seed(1234)
m1.cum_boot<-bootMer(m1.cum, FUN = fixef, nsim = 1000)
m1.cum_boot # estimate
m1.cum_boot_ci1<-boot.ci(m1.cum_boot, index =1, type = "perc")
m1.cum_boot_ci1 #CI

#site
m7.mod<-lmer(LRR_cum.abs_gamma~MAP_rescale+
               (1|study:site)-1, REML=F,data=Cum.abs_gamma, na.action="na.fail")

m7.sum<-summary(m7.mod)
m7.est<-m7.sum$coefficients

set.seed(1234)
m7.est_boot<-bootMer(m7.mod, FUN = fixef, nsim = 1000)
m7.est_boot#effect size 
m7.est_ci1<-boot.ci(m7.est_boot, index =1, type = "perc")#delta P
m7.est_ci1#CI 


##############################~~~~~~~species richness~~~~~~~~###########################

#local 
m1.mod<-lmer(LRR_S~MAP_rescale+
               (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m1.sum<-summary(m1.mod)
m1.est<-m1.sum$coefficients
m1.est

##bootstrapping effect size and CI
set.seed(1234)
m1.est_boot<-bootMer(m1.mod, FUN = fixef, nsim = 1000)
m1.est_boot#estimate
m1.est_ci1<-boot.ci(m1.est_boot, index =1, type = "perc")#deltaP
m1.est_ci1#CI 

#turnover 
m3a.mod<-lmer(LRR_betaSn~MAP_rescale+
                (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m3a.sum<-summary(m3a.mod)
m3a.est<-m3a.sum$coefficients
m3a.est

###bootstrapping

set.seed(1234)
m3a.est_boot<-bootMer(m3a.mod, FUN = fixef, nsim = 1000)
m3a.est_boot#effect size delta P
m3a.est_ci1<-boot.ci(m3a.est_boot, index =1, type = "perc")
m3a.est_ci1#CI #delta P


#site
m7.mod<-lmer(LRR_gammaSn~MAP_rescale+
               (1|study:site)-1,REML=T,data=Div.data_gamma, na.action="na.fail")

m7.sum<-summary(m7.mod)
m7.est<-m7.sum$coefficients


##bootstrapping effect size and CI
set.seed(1234)
m7.est_boot<-bootMer(m7.mod, FUN = fixef, nsim = 1000)
m7.est_boot#estimate
m7.est_ci1<-boot.ci(m7.est_boot, index =1, type = "perc")
m7.est_ci1#CI



##############################~~~~~~~~~~evenness~~~~~~~~################################
#local
m2.mod<-lmer(LRR_SPie~MAP_rescale+
               (1|study:site:SiteBlock)-1,REML=T, data=Div.data_local,na.action="na.fail")

m2.sum<-summary(m2.mod)
m2.est<-m2.sum$coefficients

###bootstrapping
set.seed(1234)
m2.est_boot<-bootMer(m2.mod, FUN = fixef, nsim = 1000)
m2.est_boot #effect size
m2.est_ci1<-boot.ci(m2.est_boot, index =1, type = "perc")#deltaP
m2.est_ci1 #CI 

#turnover
m4.mod3<-lmer(LRR_betaSPie~MAP_rescale+
                (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m4.sum3<-summary(m4.mod3)
m4.est3<-m4.sum3$coefficients

##bootstrapping
set.seed(1234)
m4.est3_boot<-bootMer(m4.mod3, FUN = fixef, nsim = 1000)
m4.est3_boot#effect size 
m4.est3_ci1<-boot.ci(m4.est3_boot, index =1, type = "perc")
m4.est3_ci1#CI 

#site
m6.mod<-lmer(LRR_gammaSPie~MAP_rescale+
               (1|study:site)-1,REML=T,data=Div.data_gamma, na.action="na.fail")

m6.sum<-summary(m6.mod)
m6.est<-m6.sum$coefficients


##bootstrapping effect size and CI
set.seed(1234)
m6.est_boot<-bootMer(m6.mod, FUN = fixef, nsim = 1000)
m6.est_boot#effect size
m6.est_ci1<-boot.ci(m6.est_boot, index =1, type = "perc")
m6.est_ci1#CI


#######################################################################################
#################################    PET    ###########################################
#######################################################################################

##############################~~~~~~~total abundance~~~~~~~~###########################
#local
m1.cum<-lmer(LRR_cum.abs~PET_rescale+
               (1|study:site:SiteBlock)-1,REML=T, data=Cum.abs_local,na.action="na.fail")

m1.cum.sum<-summary(m1.cum)
m1.cum.est<-m1.cum.sum$coefficients

set.seed(1234)
m1.cum_boot<-bootMer(m1.cum, FUN = fixef, nsim = 1000)
m1.cum_boot # estimate
m1.cum_boot_ci1<-boot.ci(m1.cum_boot, index =1, type = "perc")
m1.cum_boot_ci1 #CI

#site
m7.mod<-lmer(LRR_cum.abs_gamma~PET_rescale+
               (1|study:site)-1, REML=F,data=Cum.abs_gamma, na.action="na.fail")

m7.sum<-summary(m7.mod)
m7.est<-m7.sum$coefficients

set.seed(1234)
m7.est_boot<-bootMer(m7.mod, FUN = fixef, nsim = 1000)
m7.est_boot#effect size 
m7.est_ci1<-boot.ci(m7.est_boot, index =1, type = "perc")#delta P
m7.est_ci1#CI 


##############################~~~~~~~species richness~~~~~~~~###########################

#local 
m1.mod<-lmer(LRR_S~PET_rescale+
               (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m1.sum<-summary(m1.mod)
m1.est<-m1.sum$coefficients
m1.est

##bootstrapping effect size and CI
set.seed(1234)
m1.est_boot<-bootMer(m1.mod, FUN = fixef, nsim = 1000)
m1.est_boot#estimate
m1.est_ci1<-boot.ci(m1.est_boot, index =1, type = "perc")#deltaP
m1.est_ci1#CI 

#turnover 
m3a.mod<-lmer(LRR_betaSn~PET_rescale+
                (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m3a.sum<-summary(m3a.mod)
m3a.est<-m3a.sum$coefficients
m3a.est

###bootstrapping

set.seed(1234)
m3a.est_boot<-bootMer(m3a.mod, FUN = fixef, nsim = 1000)
m3a.est_boot#effect size delta P
m3a.est_ci1<-boot.ci(m3a.est_boot, index =1, type = "perc")
m3a.est_ci1#CI #delta P



# site 
m7.mod<-lmer(LRR_gammaSn~PET_rescale+
               (1|study:site)-1,REML=T,data=Div.data_gamma, na.action="na.fail")

m7.sum<-summary(m7.mod)
m7.est<-m7.sum$coefficients


##bootstrapping effect size and CI
set.seed(1234)
m7.est_boot<-bootMer(m7.mod, FUN = fixef, nsim = 1000)
m7.est_boot#estimate
m7.est_ci1<-boot.ci(m7.est_boot, index =1, type = "perc")
m7.est_ci1#CI



##############################~~~~~~~~~~evenness~~~~~~~~################################
#local
m2.mod<-lmer(LRR_SPie~PET_rescale+
               (1|study:site:SiteBlock)-1,REML=T, data=Div.data_local,na.action="na.fail")

m2.sum<-summary(m2.mod)
m2.est<-m2.sum$coefficients

###bootstrapping
set.seed(1234)
m2.est_boot<-bootMer(m2.mod, FUN = fixef, nsim = 1000)
m2.est_boot #effect size
m2.est_ci1<-boot.ci(m2.est_boot, index =1, type = "perc")#deltaP
m2.est_ci1 #CI 

#turnover
m4.mod3<-lmer(LRR_betaSPie~PET_rescale+
                (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m4.sum3<-summary(m4.mod3)
m4.est3<-m4.sum3$coefficients

##bootstrapping
set.seed(1234)
m4.est3_boot<-bootMer(m4.mod3, FUN = fixef, nsim = 1000)
m4.est3_boot#effect size 
m4.est3_ci1<-boot.ci(m4.est3_boot, index =1, type = "perc")
m4.est3_ci1#CI 

#site
m6.mod<-lmer(LRR_gammaSPie~PET_rescale+
               (1|study:site)-1,REML=T,data=Div.data_gamma, na.action="na.fail")

m6.sum<-summary(m6.mod)
m6.est<-m6.sum$coefficients


##bootstrapping effect size and CI
set.seed(1234)
m6.est_boot<-bootMer(m6.mod, FUN = fixef, nsim = 1000)
m6.est_boot#effect size
m6.est_ci1<-boot.ci(m6.est_boot, index =1, type = "perc")
m6.est_ci1#CI

#######################################################################################
#################################  duration  ###########################################
#######################################################################################

##############################~~~~~~~total abundance~~~~~~~~###########################
#local
m1.cum<-lmer(LRR_cum.abs~duration_rescale+
               (1|study:site:SiteBlock)-1,REML=T, data=Cum.abs_local,na.action="na.fail")

m1.cum.sum<-summary(m1.cum)
m1.cum.est<-m1.cum.sum$coefficients

set.seed(1234)
m1.cum_boot<-bootMer(m1.cum, FUN = fixef, nsim = 1000)
m1.cum_boot # estimate
m1.cum_boot_ci1<-boot.ci(m1.cum_boot, index =1, type = "perc")
m1.cum_boot_ci1 #CI

#site
m7.mod<-lmer(LRR_cum.abs_gamma~duration_rescale+
               (1|study:site)-1, REML=F,data=Cum.abs_gamma, na.action="na.fail")

m7.sum<-summary(m7.mod)
m7.est<-m7.sum$coefficients

set.seed(1234)
m7.est_boot<-bootMer(m7.mod, FUN = fixef, nsim = 1000)
m7.est_boot#effect size delta P
m7.est_ci1<-boot.ci(m7.est_boot, index =1, type = "perc")#delta P
m7.est_ci1#CI delta P


##############################~~~~~~~species richness~~~~~~~~###########################

#local 
m1.mod<-lmer(LRR_S~duration_rescale+
               (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m1.sum<-summary(m1.mod)
m1.est<-m1.sum$coefficients
m1.est

##bootstrapping effect size and CI
set.seed(1234)
m1.est_boot<-bootMer(m1.mod, FUN = fixef, nsim = 1000)
m1.est_boot#estimate
m1.est_ci1<-boot.ci(m1.est_boot, index =1, type = "perc")#deltaP
m1.est_ci1#CI 

#turnover 
m3a.mod<-lmer(LRR_betaSn~duration_rescale+
                (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m3a.sum<-summary(m3a.mod)
m3a.est<-m3a.sum$coefficients
m3a.est

###bootstrapping

set.seed(1234)
m3a.est_boot<-bootMer(m3a.mod, FUN = fixef, nsim = 1000)
m3a.est_boot#effect size delta P
m3a.est_ci1<-boot.ci(m3a.est_boot, index =1, type = "perc")
m3a.est_ci1#CI #delta P



# site 
m7.mod<-lmer(LRR_gammaSn~duration_rescale+
               (1|study:site)-1,REML=T,data=Div.data_gamma, na.action="na.fail")

m7.sum<-summary(m7.mod)
m7.est<-m7.sum$coefficients


##bootstrapping effect size and CI
set.seed(1234)
m7.est_boot<-bootMer(m7.mod, FUN = fixef, nsim = 1000)
m7.est_boot#delta P
m7.est_ci1<-boot.ci(m7.est_boot, index =1, type = "perc")
m7.est_ci1#



##############################~~~~~~~~~~evenness~~~~~~~~################################
#local
m2.mod<-lmer(LRR_SPie~duration_rescale+
               (1|study:site:SiteBlock)-1,REML=T, data=Div.data_local,na.action="na.fail")

m2.sum<-summary(m2.mod)
m2.est<-m2.sum$coefficients

###bootstrapping
set.seed(1234)
m2.est_boot<-bootMer(m2.mod, FUN = fixef, nsim = 1000)
m2.est_boot #effect size
m2.est_ci1<-boot.ci(m2.est_boot, index =1, type = "perc")#deltaP
m2.est_ci1 #CI 

#turnover
m4.mod3<-lmer(LRR_betaSPie~duration_rescale+
                (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m4.sum3<-summary(m4.mod3)
m4.est3<-m4.sum3$coefficients

##bootstrapping
set.seed(1234)
m4.est3_boot<-bootMer(m4.mod3, FUN = fixef, nsim = 1000)
m4.est3_boot#effect size delta P
m4.est3_ci1<-boot.ci(m4.est3_boot, index =1, type = "perc")
m4.est3_ci1#CI delta P

#site
m6.mod<-lmer(LRR_gammaSPie~duration_rescale+
               (1|study:site)-1,REML=T,data=Div.data_gamma, na.action="na.fail")

m6.sum<-summary(m6.mod)
m6.est<-m6.sum$coefficients


##bootstrapping effect size and CI
set.seed(1234)
m6.est_boot<-bootMer(m6.mod, FUN = fixef, nsim = 1000)
m6.est_boot#delta P
m6.est_ci1<-boot.ci(m6.est_boot, index =1, type = "perc")
m6.est_ci1#



#######################################################################################
################################# treatment direction #################################
#######################################################################################

#~~~~~~total abundnace 
#local
m1.mod2<-lmer(LRR_cum.abs~treatment.direction:delta.P1_rescale+
                (1|study:site:SiteBlock)-1,REML=T,data=Cum.abs_local, na.action="na.fail")

m1.sum2<-summary(m1.mod2)
m1.est2<-m1.sum2$coefficients
m1.est2

##bootstrapping effect size and CI
set.seed(1234)
m1.est2_boot<-bootMer(m1.mod2, FUN = fixef, nsim = 1000)
m1.est2_boot#estimates direction 
m1.est2_ci1<-boot.ci(m1.est2_boot, index =1, type = "perc")#deltaP
m1.est2_ci1#treatment.direction:delta.P
m1.est2_ci2<-boot.ci(m1.est2_boot, index =2, type = "perc")#deltaP
m1.est2_ci2##treatment.direction:delta P

m7.mod1<-lmer(LRR_cum.abs_gamma~treatment.direction:delta.P1_rescale+
                (1|study:site)-1, REML=F,data=Cum.abs_gamma, na.action="na.fail")

m7.sum1<-summary(m7.mod1)
m7.est2<-m7.sum1$coefficients

set.seed(1234)
m7.est2_boot<-bootMer(m7.mod1, FUN = fixef, nsim = 1000)
m7.est2_boot#effect size delta P
m7.est2_ci1<-boot.ci(m7.est2_boot, index =1, type = "perc")#delta P
m7.est2_ci1##treatment decr.P 
m7.est2_ci2<-boot.ci(m7.est2_boot, index =2, type = "perc")#delta P
m7.est2_ci2##treatment incr.P 

#######~~~~~~species richness~~~~~~~######
#local 

m1.mod3<-lmer(LRR_S~treatment.direction:delta.P1_rescale+
                (1|study:site:SiteBlock)-1,REML=T, data=Div.data_local,na.action="na.fail")

m1.sum3<-summary(m1.mod3)
m1.est3<-m1.sum3$coefficients

###bootstrapping
set.seed(1234)
m1.est3_boot<-bootMer(m1.mod3, FUN = fixef, nsim = 1000)
m1.est3_boot #estimates direction 
m1.est3_ci1<-boot.ci(m1.est3_boot, index =1, type = "perc")#deltaP
m1.est3_ci1 #treatment decr.P 
m1.est3_ci2<-boot.ci(m1.est3_boot, index =2, type = "perc")#deltaP
m1.est3_ci2 #treatment incr.P

#turnover
m3a.mod3<-lmer(LRR_betaSn~treatment.direction:delta.P1_rescale+
                 (1|study:site:SiteBlock)-1,REML=T, data=Div.data_local,na.action="na.fail")

m3a.sum3a<-summary(m3a.mod3)
m3a.est3<-m3a.sum3a$coefficients

###bootstrapping
set.seed(1234)
m3a.est3_boot<-bootMer(m3a.mod3, FUN = fixef, nsim = 1000)
m3a.est3_boot #estimates direction 
m3a.est3_ci1<-boot.ci(m3a.est3_boot, index =1, type = "perc")#deltaP
m3a.est3_ci1 #treatment decr.P 
m3a.est3_ci2<-boot.ci(m3a.est3_boot, index =2, type = "perc")#deltaP
m3a.est3_ci2 #treatment incr.P

#site 
m5a.mod2<-lmer(LRR_gammaSn~treatment.direction:delta.P1_rescale+
                 (1|study:site)-1,REML=T,data=Div.data_gamma,  na.action="na.fail")

m5a.sum2<-summary(m5a.mod2)
m5a.est2<-m5a.sum2$coefficients

##bootstrapping 
set.seed(1234)
m5a.est2_boot<-bootMer(m5a.mod2, FUN = fixef, nsim = 1000)
m5a.est2_boot##estimates direction 
m5a.est2_ci1<-boot.ci(m5a.est2_boot, index =1, type = "perc")
m5a.est2_ci1#treatment decr.P 
m5a.est2_ci2<-boot.ci(m5a.est2_boot, index =2, type = "perc")
m5a.est2_ci2#treatment incr.P 


#######~~~~~~evenness~~~~~~~######

m2.mod1<-lmer(LRR_SPie~treatment.direction:delta.P1_rescale+
                (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m2.sum1<-summary(m2.mod1)
m2.est1<-m2.sum1$coefficients

###bootstrapping
set.seed(1234)
m2.est1_boot<-bootMer(m2.mod1, FUN = fixef, nsim = 1000)
m2.est1_boot##effect sizes treatment direction
m2.est1_ci1<-boot.ci(m2.est1_boot, index =1, type = "perc")#MAP
m2.est1_ci2#CI treatment decr.P
m2.est1_ci2<-boot.ci(m2.est1_boot, index =2, type = "perc")#deltaP
m2.est1_ci2#CI treatmen incr.P

#turover
m4.mod2<-lmer(LRR_betaSPie~treatment.direction:delta.P1_rescale+
                (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m4.sum2<-summary(m4.mod2)
m4.est2<-m4.sum2$coefficients

##bootstrapping
set.seed(1234)
m4.est2_boot<-bootMer(m4.mod2, FUN = fixef, nsim = 1000)
m4.est2_boot#effect sizes treatment direction 
m4.est2_ci1<-boot.ci(m4.est2_boot, index =1, type = "perc")
m4.est2_ci1#CI treatment decr.P
m4.est2_ci2<-boot.ci(m4.est2_boot, index =2, type = "perc")
m4.est2_ci2#CI treatmen incr.P

#site 
m6.mod3<-lmer(LRR_gammaSPie~treatment.direction:delta.P1_rescale+
                (1|study:site)-1,REML=T,data=Div.data_gamma, na.action="na.fail")

m6.sum3<-summary(m6.mod3)
m6.est3<-m6.sum3$coefficients


##bootstrapping effect size and CI
set.seed(1234)
m6.est3_boot<-bootMer(m6.mod3, FUN = fixef, nsim = 1000)
m6.est3_boot#effect sizes treatment direction
m6.est3_ci1<-boot.ci(m6.est3_boot, index =1, type = "perc")
m6.est3_ci1#CI treatment decr.P
m6.est3_ci2<-boot.ci(m6.est3_boot, index =2, type = "perc")
m6.est3_ci2#CI treatmen incr.P


#######################################################################################
#################################  delta P : PET ##################################
#######################################################################################


##############################~~~~~~~total abundance~~~~~~~~###########################


m1.cum.best6<-lmer(LRR_cum.abs~delta.P1_rescale*PET_rescale+
                     (1|study:site:SiteBlock)-1,REML=T, data=Cum.abs_local, na.action="na.fail")

m1.cum.best6.sum<-summary(m1.cum.best6)
m1.cum.best6.est<-m1.cum.best6.sum$coefficients

set.seed(1234)
m1.cum.best6_boot<-bootMer(m1.cum.best6, FUN = fixef, nsim = 1000)
m1.cum.best6_boot#t1 = MAP, t2 = delta P, t3 = delta P * MAP
m1.cum.best6_ci1<-boot.ci(m1.cum.best6_boot, index =1, type = "perc")
m1.cum.best6_ci1#delta P
m1.cum.best6_ci2<-boot.ci(m1.cum.best6_boot, index =2, type = "perc")
m1.cum.best6_ci2#MAP
m1.cum.best6_ci3<-boot.ci(m1.cum.best6_boot, index =3, type = "perc")
m1.cum.best6_ci3#delta P*MAP

m7.mod1<-lmer(LRR_cum.abs_gamma~delta.P1_rescale*PET_rescale+
                (1|study:site)-1, REML=F,data=Cum.abs_gamma, na.action="na.fail")

m7.sum1<-summary(m7.mod1)
m7.est2<-m7.sum1$coefficients

set.seed(1234)
m7.est2_boot<-bootMer(m7.mod1, FUN = fixef, nsim = 1000)
m7.est2_boot#effect size delta P
m7.est2_ci1<-boot.ci(m7.est2_boot, index =1, type = "perc")#delta P
m7.est2_ci1#delta P
m7.est2_ci2<-boot.ci(m7.est2_boot, index =2, type = "perc")#delta P
m7.est2_ci2##PET
m7.est2_ci3<-boot.ci(m7.est2_boot, index =2, type = "perc")#delta P
m7.est2_ci3##delta P*PET


##############################~~~~~~~species richness~~~~~~~~###########################

# local 
m1.mod1<-lmer(LRR_S~PET_rescale*delta.P1_rescale+
                (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local,  na.action="na.fail")

m1.sum1<-summary(m1.mod1)
m1.est1<-m1.sum1$coefficients
m1.est1

##bootstrapping effect size and CI
set.seed(1234)
m1.est1_boot<-bootMer(m1.mod1, FUN = fixef, nsim = 1000)
m1.est1_boot#t1 = PET, t2 = delta P, t3 = elta P * PET
m1.est1_ci1<-boot.ci(m1.est1_boot, index =1, type = "perc")
m1.est1_ci1#PET
m1.est1_ci2<-boot.ci(m1.est1_boot, index =2, type = "perc")
m1.est1_ci2#deltaP
m1.est1_ci3<-boot.ci(m1.est1_boot, index =3, type = "perc")
m1.est1_ci3#PET*delta P


# turover 
m3a.mod1<-lmer(LRR_betaSn~PET_rescale*delta.P1_rescale+
                 (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m3a.sum1<-summary(m3a.mod1)
m3a.est1<-m3a.sum1$coefficients
m3a.est1

###bootstrapping
set.seed(1234)
m3a.est1_boot<-bootMer(m3a.mod1, FUN = fixef, nsim = 1000)
m3a.est1_boot#effect size delta P
m3a.est1_ci1<-boot.ci(m3a.est1_boot, index =1, type = "perc")#delta P
m3a.est1_ci1#CI PET
m3a.est1_ci2<-boot.ci(m3a.est1_boot, index =2, type = "perc")#delta P
m3a.est1_ci2#CI delta
m3a.est1_ci3<-boot.ci(m3a.est1_boot, index =3, type = "perc")#delta P
m3a.est1_ci3#CI delta : PET

#site 
m7.mod1<-lmer(LRR_gammaSn~PET_rescale*delta.P1_rescale+
                (1|study:site)-1,REML=T,data=Div.data_gamma, na.action="na.fail")

m7.sum1<-summary(m7.mod1)
m7.est1<-m7.sum1$coefficients


##bootstrapping effect size and CI
set.seed(1234)
m7.est1_boot<-bootMer(m7.mod1, FUN = fixef, nsim = 1000)
m7.est1_boot#delta P
m7.est1_ci1<-boot.ci(m7.est1_boot, index =1, type = "perc")
m7.est1_ci1#
m7.est1_ci2<-boot.ci(m7.est1_boot, index =2, type = "perc")
m7.est1_ci2#
m7.est1_ci3<-boot.ci(m7.est1_boot, index =3, type = "perc")
m7.est1_ci3#


##############################~~~~~~~~~~evenness~~~~~~~~################################

#local
m2.mod1<-lmer(LRR_SPie~delta.P1_rescale*PET_rescale+
                (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m2.sum1<-summary(m2.mod1)
m2.est1<-m2.sum1$coefficients

###bootstrapping
set.seed(1234)
m2.est1_boot<-bootMer(m2.mod1, FUN = fixef, nsim = 1000)
m2.est1_boot#t1 = PET, t2 = delta P, t3 = delta P * PET
m2.est1_ci1<-boot.ci(m2.est1_boot, index =1, type = "perc")#PET
m2.est1_ci2#CI delta P
m2.est1_ci2<-boot.ci(m2.est1_boot, index =2, type = "perc")#deltaP
m2.est1_ci2#CI PET
m2.est1_ci3<-boot.ci(m2.est1_boot, index =3, type = "perc")#deltaP*PET
m2.est1_ci3#CI delta P:PET


#turnover
m4.mod<-lmer(LRR_betaSPie~PET_rescale*delta.P1_rescale+
               (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m4.sum<-summary(m4.mod)
m4.est<-m4.sum$coefficients

##bootstrapping 
set.seed(1234)
m4.est_boot<-bootMer(m4.mod, FUN = fixef, nsim = 1000)
m4.est_boot# t1 = PET, t2 = delta P, t3 = delta P * PET
m4.est_ci1<-boot.ci(m4.est_boot, index =1, type = "perc")
m4.est_ci1#PET
m4.est_ci2<-boot.ci(m4.est_boot, index =2, type = "perc")
m4.est_ci2#delta P
m4.est_ci3<-boot.ci(m4.est_boot, index =3, type = "perc")
m4.est_ci3#delta P * PET

#site
m6.mod<-lmer(LRR_gammaSPie~PET_rescale*delta.P1_rescale+
               (1|study:site)-1,REML=T,data=Div.data_gamma, na.action="na.fail")

m6.sum<-summary(m6.mod)
m6.est<-m6.sum$coefficients

##bootstrapping 
set.seed(1234)
m6.est_boot<-bootMer(m6.mod, FUN = fixef, nsim = 1000)
m6.est_boot# t1 = PET, t2 = delta P, t3 = delta P * PET
m6.est_ci1<-boot.ci(m6.est_boot, index =1, type = "perc")
m6.est_ci1#PET
m6.est_ci2<-boot.ci(m6.est_boot, index =2, type = "perc")
m6.est_ci2#delta P
m6.est_ci3<-boot.ci(m6.est_boot, index =3, type = "perc")
m6.est_ci3#delta P * PET


#######################################################################################
#################################  delta P : MAP ##################################
#######################################################################################


##############################~~~~~~~total abundance~~~~~~~~###########################


m1.cum.best6<-lmer(LRR_cum.abs~delta.P1_rescale*MAP_rescale+
                     (1|study:site:SiteBlock)-1,REML=T, data=Cum.abs_local, na.action="na.fail")

m1.cum.best6.sum<-summary(m1.cum.best6)
m1.cum.best6.est<-m1.cum.best6.sum$coefficients

set.seed(1234)
m1.cum.best6_boot<-bootMer(m1.cum.best6, FUN = fixef, nsim = 1000)
m1.cum.best6_boot#t1 = MAP, t2 = delta P, t3 = delta P * MAP
m1.cum.best6_ci1<-boot.ci(m1.cum.best6_boot, index =1, type = "perc")
m1.cum.best6_ci1#delta P
m1.cum.best6_ci2<-boot.ci(m1.cum.best6_boot, index =2, type = "perc")
m1.cum.best6_ci2#MAP
m1.cum.best6_ci3<-boot.ci(m1.cum.best6_boot, index =3, type = "perc")
m1.cum.best6_ci3#delta P*MAP


m7.mod1<-lmer(LRR_cum.abs_gamma~delta.P1_rescale*MAP_rescale+
                (1|study:site)-1, REML=F,data=Cum.abs_gamma, na.action="na.fail")

m7.sum1<-summary(m7.mod1)
m7.est2<-m7.sum1$coefficients

set.seed(1234)
m7.est2_boot<-bootMer(m7.mod1, FUN = fixef, nsim = 1000)
m7.est2_boot#effect size delta P
m7.est2_ci1<-boot.ci(m7.est2_boot, index =1, type = "perc")#delta P
m7.est2_ci1#delta P
m7.est2_ci2<-boot.ci(m7.est2_boot, index =2, type = "perc")#delta P
m7.est2_ci2##MAP
m7.est2_ci3<-boot.ci(m7.est2_boot, index =2, type = "perc")#delta P
m7.est2_ci3##delta P*MAP


##############################~~~~~~~species richness~~~~~~~~###########################

# local 
m1.mod1<-lmer(LRR_S~MAP_rescale*delta.P1_rescale+
                (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local,  na.action="na.fail")

m1.sum1<-summary(m1.mod1)
m1.est1<-m1.sum1$coefficients
m1.est1

##bootstrapping effect size and CI
set.seed(1234)
m1.est1_boot<-bootMer(m1.mod1, FUN = fixef, nsim = 1000)
m1.est1_boot#t1 = MAP, t2 = delta P, t3 = elta P * MAP
m1.est1_ci1<-boot.ci(m1.est1_boot, index =1, type = "perc")
m1.est1_ci1#MAP
m1.est1_ci2<-boot.ci(m1.est1_boot, index =2, type = "perc")
m1.est1_ci2#deltaP
m1.est1_ci3<-boot.ci(m1.est1_boot, index =3, type = "perc")
m1.est1_ci3#MAP*delta P


# turover 
m3a.mod1<-lmer(LRR_betaSn~MAP_rescale*delta.P1_rescale+
                 (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m3a.sum1<-summary(m3a.mod1)
m3a.est1<-m3a.sum1$coefficients
m3a.est1

###bootstrapping
set.seed(1234)
m3a.est1_boot<-bootMer(m3a.mod1, FUN = fixef, nsim = 1000)
m3a.est1_boot#effect size delta P
m3a.est1_ci1<-boot.ci(m3a.est1_boot, index =1, type = "perc")#delta P
m3a.est1_ci1#CI MAP
m3a.est1_ci2<-boot.ci(m3a.est1_boot, index =2, type = "perc")#delta P
m3a.est1_ci2#CI delta
m3a.est1_ci3<-boot.ci(m3a.est1_boot, index =3, type = "perc")#delta P
m3a.est1_ci3#CI delta : MAP

#site
m7.mod1<-lmer(LRR_gammaSn~MAP_rescale*delta.P1_rescale+
                (1|study:site)-1,REML=T,data=Div.data_gamma1, na.action="na.fail")

m7.sum1<-summary(m7.mod1)
m7.est1<-m7.sum1$coefficients


##bootstrapping effect size and CI
set.seed(1234)
m7.est1_boot<-bootMer(m7.mod1, FUN = fixef, nsim = 1000)
m7.est1_boot#delta P
m7.est1_ci1<-boot.ci(m7.est1_boot, index =1, type = "perc")
m7.est1_ci1#
m7.est1_ci2<-boot.ci(m7.est1_boot, index =2, type = "perc")
m7.est1_ci2#
m7.est1_ci3<-boot.ci(m7.est1_boot, index =3, type = "perc")
m7.est1_ci3#


##############################~~~~~~~~~~evenness~~~~~~~~################################

#local
m2.mod1<-lmer(LRR_SPie~delta.P1_rescale*MAP_rescale+
                (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m2.sum1<-summary(m2.mod1)
m2.est1<-m2.sum1$coefficients

###bootstrapping
set.seed(1234)
m2.est1_boot<-bootMer(m2.mod1, FUN = fixef, nsim = 1000)
m2.est1_boot#t1 = MAT, t2 = delta P, t3 = delta P * MAT
m2.est1_ci1<-boot.ci(m2.est1_boot, index =1, type = "perc")#MAT
m2.est1_ci2#CI delta P
m2.est1_ci2<-boot.ci(m2.est1_boot, index =2, type = "perc")#deltaP
m2.est1_ci2#CI MAT
m2.est1_ci3<-boot.ci(m2.est1_boot, index =3, type = "perc")#deltaP*MAT
m2.est1_ci3#CI delta P:MAT


#turnover
m4.mod<-lmer(LRR_betaSPie~MAP_rescale*delta.P1_rescale+
               (1|study:site:SiteBlock)-1,REML=T,data=Div.data_local, na.action="na.fail")

m4.sum<-summary(m4.mod)
m4.est<-m4.sum$coefficients

##bootstrapping 
set.seed(1234)
m4.est_boot<-bootMer(m4.mod, FUN = fixef, nsim = 1000)
m4.est_boot# t1 = MAT, t2 = delta P, t3 = delta P * MAT
m4.est_ci1<-boot.ci(m4.est_boot, index =1, type = "perc")
m4.est_ci1#MAT
m4.est_ci2<-boot.ci(m4.est_boot, index =2, type = "perc")
m4.est_ci2#delta P
m4.est_ci3<-boot.ci(m4.est_boot, index =3, type = "perc")
m4.est_ci3#delta P * MAT

#site
m6.mod<-lmer(LRR_gammaSPie~MAP_rescale*delta.P1_rescale+
               (1|study:site)-1,REML=T,data=Div.data_gamma, na.action="na.fail")

m6.sum<-summary(m6.mod)
m6.est<-m6.sum$coefficients

##bootstrapping 
set.seed(1234)
m6.est_boot<-bootMer(m6.mod, FUN = fixef, nsim = 1000)
m6.est_boot# t1 = MAT, t2 = delta P, t3 = delta P * MAT
m6.est_ci1<-boot.ci(m6.est_boot, index =1, type = "perc")
m6.est_ci1#MAT
m6.est_ci2<-boot.ci(m6.est_boot, index =2, type = "perc")
m6.est_ci2#delta P
m6.est_ci3<-boot.ci(m6.est_boot, index =3, type = "perc")
m6.est_ci3#delta P * MAT




