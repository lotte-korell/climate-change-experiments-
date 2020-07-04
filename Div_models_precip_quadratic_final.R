getwd()
setwd("C:/Users/korell/R/Meta.analysis/results")
dir()


library(lme4)
library(MASS)
library(car)
library(ggplot2)
library(MuMIn)
library(arm)



Div.data_local <- read.csv("Precip_data_final.csv" ,dec=".", sep=";",h=T)
Div.data_gamma <- read.csv("Precip_data_gamma_final.csv" ,dec=".", sep=";",h=T)
Cum.abs_local <- read.csv("ES.cum.abs_final.csv" ,dec=".", sep=";",h=T)
Cum.abs_gamma<- read.csv("ES.cum.abs.gamma_final.csv" ,dec=".", sep=";",h=T)

##convert duration into factor to include it as random effect in the model
Div.data_local$treat.duration.m<-as.factor(Div.data_local$treat.duration.m)
Div.data_gamma$treat.duration.m<-as.factor(Div.data_gamma$treat.duration.m)
Cum.abs_local$treat.duration.m<-as.factor(Cum.abs_local$treat.duration.m)
Cum.abs_gamma$treat.duration.m<-as.factor(Cum.abs_gamma$treat.duration.m)


###########################
####### local scale #######
###########################

###############################
########## cum abs ######
###############################


##model selection

m1.cum.abs_alt<-lmer(LRR_cum.abs~MAT_rescale*delta.P1_rescale+MAP_rescale*delta.P1_rescale
                     +MAT_rescale*I(delta.P1_rescale^2)+MAP_rescale*I(delta.P1_rescale^2) +
                   (1|study:site:SiteBlock:treat.duration.m), REML=F,data=Cum.abs_local, na.action="na.fail")

d1.cum.abs_alt<-dredge(m1.cum.abs_alt)#, m.lim=c(0,3))
print(d1.cum.abs_alt)

d1.cum.absAIC_alt<-model.avg(d1.cum.abs_alt, subset = delta < 2)
summary(d1.cum.absAIC_alt)

###best models

m1.cum.abs<-lmer(LRR_cum.abs~delta.P1_rescale+
                   (1|study:site:SiteBlock:treat.duration.m), REML=T,data=Cum.abs_local, na.action="na.fail")
summary(m1.cum.abs)
anova(m1.cum.abs)
Anova(m1.cum.abs)
r.squaredGLMM(m1.cum.abs)

#Response: LRR_cum.abs
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale 28.27  1  1.055e-07 ***

m2.cum.abs<-lmer(LRR_cum.abs~delta.P1_rescale+I(delta.P1_rescale^2) +MAT_rescale*I(delta.P1_rescale^2)+
                   (1|study:site:SiteBlock:treat.duration.m), REML=T,data=Cum.abs_local,na.action="na.fail")
summary(m2.cum.abs)
anova(m2.cum.abs)
Anova(m2.cum.abs)
r.squaredGLMM(m2.cum.abs)

#Response: LRR_cum.abs
#Chisq Df Pr(>Chisq)   
#delta.P1_rescale                  7.2466  1   0.007104 **
#I(delta.P1_rescale^2)             1.5720  1   0.209912   
#MAT_rescale                       0.2207  1   0.638496   
#I(delta.P1_rescale^2):MAT_rescale 4.3470  1   0.037074 * 

m3.cum.abs<-lmer(LRR_cum.abs~delta.P1_rescale+I(delta.P1_rescale^2)+
                   (1|study:site:SiteBlock:treat.duration.m), REML=T,data=Cum.abs_local,na.action="na.fail")
summary(m3.cum.abs)
anova(m3.cum.abs)
Anova(m3.cum.abs)
r.squaredGLMM(m3.cum.abs)

#Response: LRR_cum.abs
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale      22.7668  1  1.829e-06 ***
#I(delta.P1_rescale^2)  1.4086  1     0.2353         


m4.cum.abs<-lmer(LRR_cum.abs~delta.P1_rescale+MAP_rescale+
                   (1|study:site:SiteBlock:treat.duration.m), REML=T,data=Cum.abs_local,na.action="na.fail")
summary(m4.cum.abs)
anova(m4.cum.abs)
Anova(m4.cum.abs)
r.squaredGLMM(m4.cum.abs)

#Response: LRR_cum.abs
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale 28.4576  1  9.577e-08 ***
#MAP_rescale       0.4431  1     0.5056     



m5.cum.abs<-lmer(LRR_cum.abs~delta.P1_rescale+I(delta.P1_rescale^2)+MAP_rescale+
                   (1|study:site:SiteBlock:treat.duration.m), REML=T,data=Cum.abs_local,na.action="na.fail")
summary(m5.cum.abs)
anova(m5.cum.abs)
Anova(m5.cum.abs)
r.squaredGLMM(m5.cum.abs)

#Response: LRR_cum.abs
#Chisq Df Pr(>Chisq)   
#MAT_rescale                       0.2207  1   0.638496   
#delta.P1_rescale                  7.2466  1   0.007104 **
#I(delta.P1_rescale^2)             1.5720  1   0.209912   
#MAT_rescale:I(delta.P1_rescale^2) 4.3470  1   0.037074 * 


m6.cum.abs<-lmer(LRR_cum.abs~MAT_rescale+delta.P1_rescale+I(delta.P1_rescale^2) +MAT_rescale:I(delta.P1_rescale^2)+MAT_rescale:I(delta.P1_rescale^2)+
                   (1|study:site:SiteBlock:treat.duration.m), REML=T, data=Cum.abs_local, na.action="na.fail")
summary(m6.cum.abs)
summary(m6.cum.abs)
anova(m6.cum.abs)
Anova(m6.cum.abs)
r.squaredGLMM(m5.cum.abs)


#Response: LRR_cum.abs
#Chisq Df Pr(>Chisq)   
#MAT_rescale                       0.2207  1   0.638496   
#delta.P1_rescale                  7.2466  1   0.007104 **
#I(delta.P1_rescale^2)             1.5720  1   0.209912   
#MAT_rescale:I(delta.P1_rescale^2) 4.3470  1   0.037074 * 

m7.cum.abs<-lmer(LRR_cum.abs~delta.P1_rescale+MAT_rescale+
                   (1|study:site:SiteBlock:treat.duration.m), REML=T,data=Cum.abs_local,na.action="na.fail")
summary(m7.cum.abs)
anova(m7.cum.abs)
Anova(m7.cum.abs)
r.squaredGLMM(m7.cum.abs)

#Response: LRR_cum.abs
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale 25.9557  1  3.493e-07 ***
#MAT_rescale       0.1036  1     0.7475    




###############################
####### Species richness ######
###############################


####### Precip ######


#### model selection
m1.beta_alt<-lmer(LRR_S~ MAT_rescale*delta.P1_rescale+MAP_rescale*delta.P1_rescale
                  +MAT_rescale*I(delta.P1_rescale^2)+MAP_rescale*I(delta.P1_rescale^2) +
                  (1|study:site:SiteBlock:treat.duration.m),REML=F,data=Div.data_local, na.action="na.fail")

d1.beta_alt<-dredge(m1.beta_alt)
print(m1.beta_alt)

d1.AIC_alt<-model.avg(d1.beta_alt, subset = delta < 2)
summary(d1.AIC_alt)


#####best models 

m1.best1<-lmer(LRR_S~delta.P1_rescale*MAP_rescale+
                 (1|study:site:SiteBlock:treat.duration.m), REML=T,data=Div.data_local, na.action="na.fail")
anova(m1.best1)
Anova(m1.best1)
r.squaredGLMM(m1.best1)
qqnorm(resid(m1.best1))
qqline(resid(m1.best1))
hist(resid(m1.best1))

#Response: LRR_S
#Chisq Df Pr(>Chisq)   
#delta.P1_rescale             10.4838  1   0.001204 **
#MAP_rescale                   0.0084  1   0.927068   
#delta.P1_rescale:MAP_rescale  4.7139  1   0.029919 * 


m1.best2<-lmer(LRR_S~delta.P1_rescale+
                 (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m1.best2)
Anova(m1.best2)
r.squaredGLMM(m1.best2)
qqnorm(resid(m1.best2))
qqline(resid(m1.best2))
hist(resid(m1.best2))

#Response: LRR_S
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale 10.918  1  0.0009524 ***
  

m1.best3<-lmer(LRR_S~delta.P1_rescale*MAP_rescale+I(delta.P1_rescale^2)+
                 (1|study:site:SiteBlock:treat.duration.m), REML=T,data=Div.data_local, na.action="na.fail")
anova(m1.best3)
Anova(m1.best3)
r.squaredGLMM(m1.best3)
qqnorm(resid(m1.best3))
qqline(resid(m1.best3))
hist(resid(m1.best3))

#Response: LRR_S
#Chisq Df Pr(>Chisq)   
#elta.P1_rescale             10.3684  1   0.001282 **
#MAP_rescale                   0.0075  1   0.931075   
#I(delta.P1_rescale^2)         0.5864  1   0.443808   
#delta.P1_rescale:MAP_rescale  5.2952  1   0.021384 * 

 

###############################
#######     Evenness    #######
###############################

##model selection


m2.beta_alt<-lmer(LRR_SPie~ MAT_rescale*delta.P1_rescale+MAP_rescale*delta.P1_rescale+MAT_rescale*I(delta.P1_rescale^2)+MAP_rescale*I(delta.P1_rescale^2) +
                (1|study:site:SiteBlock:treat.duration.m),REML=F,data=Div.data_local, na.action="na.fail")

##########---- multimodel inference----###########

d2.beta_alt<-dredge(m2.beta_alt)
print(d2.beta_alt)

d2.AIC_alt<-model.avg(d2.beta_alt, subset = delta < 2)
summary(d2.AIC_alt)

####best models

m2.best<-lmer(LRR_SPie~delta.P1_rescale+
                (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m2.best)
Anova(m2.best)
r.squaredGLMM(m2.best) 
qqnorm(resid(m2.best))
qqline(resid(m2.best))
hist(resid(m2.best))

#Response: LRR_SPie
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale 13.914  1  0.0001914 ***

m2.best1<-lmer(LRR_SPie~MAP_rescale*delta.P1_rescale+
                (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m2.best1)
Anova(m2.best1)
r.squaredGLMM(m2.best1) 
qqnorm(resid(m2.best1))
qqline(resid(m2.best1))
hist(resid(m2.best1))

#RResponse: LRR_SPie
#Chisq Df Pr(>Chisq)    
#MAP_rescale                   0.3688  1  0.5436798    
#delta.P1_rescale             13.8606  1  0.0001969 ***
#MAP_rescale:delta.P1_rescale  3.6135  1  0.0573127 .  


m2.best2<-lmer(LRR_SPie~delta.P1_rescale+I(delta.P1_rescale^2)+
                 (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m2.best2)
Anova(m2.best2)
r.squaredGLMM(m2.best2)
qqnorm(resid(m2.best2))
qqline(resid(m2.best2))
hist(resid(m2.best2))


#Response: LRR_SPie
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale      11.7522  1  0.0006077 ***
#I(delta.P1_rescale^2)  0.9937  1  0.3188514    


m2.best3<-lmer(LRR_SPie~delta.P1_rescale+MAP_rescale*I(delta.P1_rescale^2)+MAT_rescale+
                 (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m2.best3)
Anova(m2.best3)
r.squaredGLMM(m2.best3)
qqnorm(resid(m2.best3))
qqline(resid(m2.best3))
hist(resid(m2.best3))

#Response: LRR_SPie
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale                  14.6123  1   0.000132 ***
#MAP_rescale                        0.6845  1   0.408045    
#I(delta.P1_rescale^2)              1.1853  1   0.276281    
#MAT_rescale                        0.1674  1   0.682448    
#MAP_rescale:I(delta.P1_rescale^2)  2.7688  1   0.096119 .  



m2.best4<-lmer(LRR_SPie~MAP_rescale+delta.P1_rescale+
                 (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local,na.action="na.fail")
anova(m2.best4)
Anova(m2.best4)
r.squaredGLMM(m2.best4)
qqnorm(resid(m2.best4))
qqline(resid(m2.best4))
hist(resid(m2.best4))

#RResponse: LRR_SPie
#Chisq Df Pr(>Chisq)    
#MAP_rescale       0.4072  1  0.5234126    
#delta.P1_rescale 14.2807  1  0.0001575 ***


m2.best5<-lmer(LRR_SPie~MAP_rescale*delta.P1_rescale+I(delta.P1_rescale^2)+
                 (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local,na.action="na.fail")
anova(m2.best5)
Anova(m2.best5)
r.squaredGLMM(m2.best5)
qqnorm(resid(m2.best5))
qqline(resid(m2.best5))
hist(resid(m2.best5))


#Response: LRR_SPie
#Chisq Df Pr(>Chisq)    
#MAP_rescale                   0.5913  1  0.4419274    
#delta.P1_rescale             12.1767  1  0.0004839 ***
#I(delta.P1_rescale^2)         0.2450  1  0.6205882    
#MAP_rescale:delta.P1_rescale  2.6505  1  0.1035159    

m2.best7<-lmer(LRR_SPie~MAT_rescale+delta.P1_rescale+
                 (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local,na.action="na.fail")
anova(m2.best7)
Anova(m2.best7)
r.squaredGLMM(m2.best7)
qqnorm(resid(m2.best7))
qqline(resid(m2.best7))
hist(resid(m2.best7))

#Response: LRR_SPie
#Chisq Df Pr(>Chisq)    
#MAT_rescale       0.0845  1  0.7712976    
#delta.P1_rescale 13.4285  1  0.0002478 ***



##############################
######turnover sclae #########
##############################

###############################
##### Species richness ######
##############################


###model selection

m3.beta_alt<-lmer(LRR_betaS~MAT_rescale * delta.P1_rescale + MAT_rescale * I(delta.P1_rescale^2) +  
                    MAP_rescale * delta.P1_rescale + MAP_rescale * I(delta.P1_rescale^2) + 
                (1|study:site:SiteBlock:treat.duration.m),REML=F,data=Div.data_local, na.action="na.fail")


##########---- multimodel inference----###########

d3.beta_alt<-dredge(m3.beta_alt)
print(d3.beta_alt)

d3.AIC_alt<-model.avg(d3.beta_alt, subset = delta < 2)
summary(d3.AIC_alt)


###best  models

m3.best<-lmer(LRR_betaS~delta.P1_rescale+
                (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m3.best)
Anova(m3.best)
r.squaredGLMM(m3.best) 
qqnorm(resid(m3.best))
qqline(resid(m3.best))
hist(resid(m3.best))

#Response: LRR_betaS
#Chisq Df Pr(>Chisq)  
#delta.P1_rescale 2.8673  1     0.0904 .


m3.best2<-lmer(LRR_betaS~delta.P1_rescale+MAP_rescale+
                 (1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m3.best2)
Anova(m3.best2)
r.squaredGLMM(m3.best2) 
qqnorm(resid(m3.best2))
qqline(resid(m3.best2))
hist(resid(m3.best2))

#Response: LRR_betaS
#Chisq Df Pr(>Chisq)
#delta.P1_rescale 2.3685  1     0.1238
#MAP_rescale      0.2820  1     0.5954

######################
####### ENS PIE#######
#####################

###model selection


m4.beta_alt<-lmer(LRR_betaSPie ~ MAT_rescale * delta.P1_rescale + MAT_rescale * I(delta.P1_rescale^2) +
                MAP_rescale * delta.P1_rescale + MAP_rescale * I(delta.P1_rescale^2) 
                +(1|study:site:SiteBlock:treat.duration.m),REML=F,data=Div.data_local, na.action="na.fail")


##########---- multimodel inference----###########

d4.beta_alt<-dredge(m4.beta_alt)
print(d4.beta_alt)

d4.AIC_alt<-model.avg(d4.beta_alt, subset = delta < 3)
summary(d4.AIC_alt)

###best models 

m4.best<-lmer(LRR_betaSPie~ delta.P1_rescale +
                MAP_rescale * I(delta.P1_rescale^2) + 
                MAT_rescale* I(delta.P1_rescale^2)+
                +(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best)
Anova(m4.best)
r.squaredGLMM(m4.best)
qqnorm(resid(m4.best))
qqline(resid(m4.best))
hist(resid(m4.best))

#Response: LRR_betaSPie
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale                  18.0306  1  2.174e-05 ***
# MAP_rescale                        1.4174  1  0.2338373    
#I(delta.P1_rescale^2)              0.2795  1  0.5970476    
#MAT_rescale                        0.9702  1  0.3246423    
#MAP_rescale:I(delta.P1_rescale^2) 11.7560  1  0.0006065 ***
#I(delta.P1_rescale^2):MAT_rescale  1.3131  1  0.2518340    


m4.best<-lmer(LRR_betaSPie~ delta.P1_rescale + MAP_rescale * I(delta.P1_rescale^2)   +
                +(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best)
Anova(m4.best)
r.squaredGLMM(m4.best)
qqnorm(resid(m4.best))
qqline(resid(m4.best))
hist(resid(m4.best))


#Response: LRR_betaSPie
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale                  19.4893  1  1.012e-05 ***
#MAP_rescale                        1.0410  1   0.307584    
#I(delta.P1_rescale^2)              0.2051  1   0.650661    
#MAP_rescale:I(delta.P1_rescale^2)  9.9931  1   0.001571 ** 



m4.best1<-lmer(LRR_betaSPie~ MAP_rescale *delta.P1_rescale + MAP_rescale * I(delta.P1_rescale^2) +
                 +(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best1)
Anova(m4.best1)
r.squaredGLMM(m4.best1)
qqnorm(resid(m4.best1))
qqline(resid(m4.best1))
hist(resid(m4.best1))

#Response: LRR_betaSPie
#Chisq Df Pr(>Chisq)    
#MAP_rescale                        1.1890  1    0.27554    
#delta.P1_rescale                  16.7613  1  4.239e-05 ***
#I(delta.P1_rescale^2)              1.6180  1    0.20337    
#MAT_rescale                        2.2189  1    0.13633    
#MAP_rescale:delta.P1_rescale       2.0488  1    0.15233    
#MAP_rescale:I(delta.P1_rescale^2)  4.7814  1    0.02877 *  

m4.best2<-lmer(LRR_betaSPie~ delta.P1_rescale + MAP_rescale * I(delta.P1_rescale^2) +MAT_rescale+
                 +(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best2)
Anova(m4.best2)
r.squaredGLMM(m4.best2)
qqnorm(resid(m4.best2))
qqline(resid(m4.best2))
hist(resid(m4.best2))

#Response: LRR_betaSPie
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale                  20.3794  1  6.351e-06 ***
#MAP_rescale                        1.0546  1  0.3044423    
#I(delta.P1_rescale^2)              0.2803  1  0.5964938    
#MAT_rescale                        0.9704  1  0.3245765    
#MAP_rescale:I(delta.P1_rescale^2) 10.9237  1  0.0009494 ***


m4.best3<-lmer(LRR_betaSPie~ delta.P1_rescale + MAP_rescale * I(delta.P1_rescale^2)  +MAT_rescale
                +(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best3)
Anova(m4.best3)
r.squaredGLMM(m4.best3)
qqnorm(resid(m4.best3))
qqline(resid(m4.best3))
hist(resid(m4.best3))

#Response: LRR_betaSPie
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale                  20.3794  1  6.351e-06 ***
#MAP_rescale                        1.0546  1  0.3044423    
#I(delta.P1_rescale^2)              0.2803  1  0.5964938    
#MAT_rescale                        0.9704  1  0.3245765    
#MAP_rescale:I(delta.P1_rescale^2) 10.9237  1  0.0009494 ***


m4.best4<-lmer(LRR_betaSPie~ MAT_rescale + MAP_rescale *delta.P1_rescale + MAP_rescale * I(delta.P1_rescale^2)
                 +(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best4)
Anova(m4.best4)
r.squaredGLMM(m4.best4)
qqnorm(resid(m4.best4))
qqline(resid(m4.best4))
hist(resid(m4.best4))

#Response: LRR_betaSPie
#Chisq Df Pr(>Chisq)    
#MAT_rescale                        2.1706  1    0.14068    
#MAP_rescale                        0.9725  1    0.32405    
#delta.P1_rescale                  20.1054  1  7.329e-06 ***
#I(delta.P1_rescale^2)              2.6602  1    0.10289    
#MAP_rescale:delta.P1_rescale       2.4162  1    0.12008    
#MAP_rescale:I(delta.P1_rescale^2)  5.0278  1    0.02494 *  
  ---



m4.best6<-lmer(LRR_betaSPie~MAT_rescale * I(delta.P1_rescale^2) +MAT_rescale * delta.P1_rescale+
                 MAP_rescale * delta.P1_rescale + MAP_rescale * I(delta.P1_rescale^2) 
                 +(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best6)
Anova(m4.best6)
r.squaredGLMM(m4.best6) 
qqnorm(resid(m4.best6))
hist(resid(m4.best6))

#Response: LRR_betaSPie
#Chisq Df Pr(>Chisq)    
#MAT_rescale                        2.1504  1   0.142530    
#I(delta.P1_rescale^2)              2.2998  1   0.129395    
#delta.P1_rescale                  17.8098  1  2.441e-05 ***
#MAP_rescale                        0.6025  1   0.437637    
#MAT_rescale:I(delta.P1_rescale^2)  3.3759  1   0.066155 .  
#MAT_rescale:delta.P1_rescale       2.1622  1   0.141446    
#delta.P1_rescale:MAP_rescale       0.3266  1   0.567665    
#I(delta.P1_rescale^2):MAP_rescale  7.5495  1   0.006003 ** 


m4.best7<-lmer(LRR_betaSPie~delta.P1_rescale+MAT_rescale * I(delta.P1_rescale^2) + MAP_rescale * I(delta.P1_rescale^2)+ 
               +(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best7)
Anova(m4.best7)
r.squaredGLMM(m4.best7) 
qqnorm(resid(m4.best7))
hist(resid(m4.best7))

#Response: LRR_betaSPie
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale                  18.0306  1  2.174e-05 ***
#MAT_rescale                        0.9702  1  0.3246423    
#I(delta.P1_rescale^2)              0.2795  1  0.5970476    
#MAP_rescale                        1.4174  1  0.2338373    
#MAT_rescale:I(delta.P1_rescale^2)  1.3131  1  0.2518340    
#I(delta.P1_rescale^2):MAP_rescale 11.7560  1  0.0006065 ***



m4.best8<-lmer(LRR_betaSPie~MAT_rescale * delta.P1_rescale +
                  MAP_rescale * I(delta.P1_rescale^2) 
                +(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best8)
Anova(m4.best8)
r.squaredGLMM(m4.best8)
qqnorm(resid(m4.best8))
qqline(resid(m4.best8))
hist(resid(m4.best8))

#Response: LRR_betaSPie
#Chisq Df Pr(>Chisq)    
#MAT_rescale                        0.9642  1  0.3261311    
#delta.P1_rescale                  20.3508  1  6.447e-06 ***
#MAP_rescale                        0.4541  1  0.5003728    
#I(delta.P1_rescale^2)              0.0090  1  0.9245465    
#MAT_rescale:delta.P1_rescale       1.5933  1  0.2068624    
#MAP_rescale:I(delta.P1_rescale^2) 12.0454  1  0.0005192 ***


m4.best9<-lmer(LRR_betaSPie~ MAT_rescale * I(delta.P1_rescale^2) + MAP_rescale * I(delta.P1_rescale^2) 
               +(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
anova(m4.best9)
Anova(m4.best9)
r.squaredGLMM(m4.best9)
qqnorm(resid(m4.best9))
qqline(resid(m4.best9))
hist(resid(m4.best9))

#Response: LRR_betaSPie
#Chisq Df Pr(>Chisq)  
#MAT_rescale                       0.1118  1     0.7381  
#I(delta.P1_rescale^2)             0.0366  1     0.8482  
#MAP_rescale                       0.0147  1     0.9035  
#MAT_rescale:I(delta.P1_rescale^2) 3.2165  1     0.0729 .
#I(delta.P1_rescale^2):MAP_rescale 2.6141  1     0.1059  


####################
######gamma#########
####################

####################
##### cum abs ######
####################


##model selection

m7.cum.abs_alt<-lmer(LRR_cum.abs_gamma~MAT_rescale * delta.P1_rescale + MAT_rescale * I(delta.P1_rescale^2) +  
                       MAP_rescale * delta.P1_rescale + MAP_rescale * I(delta.P1_rescale^2)
                     +(1|study:site:treat.duration.m),REML=F,data=Cum.abs_gamma, na.action="na.fail")

d7.beta_alt<-dredge(m7.cum.abs_alt)
print(m7.cum.abs_alt)

d7.AIC_alt<-model.avg(d7.beta_alt, subset = delta < 2)
summary(d7.AIC_alt)


###best models 

m7.cum.abs1<-lmer(LRR_cum.abs_gamma~delta.P1_rescale
                  +(1|study:site:treat.duration.m), REML=T,data=Cum.abs_gamma, na.action="na.fail")
anova(m7.cum.abs1)
Anova(m7.cum.abs1)
r.squaredGLMM(m7.cum.abs1)

#Response: LRR_cum.abs_gamma
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale 17.08  1  3.584e-05 ***

m7.cum.abs2<-lmer(LRR_cum.abs_gamma~MAT_rescale + delta.P1_rescale 
                  + (1|study:site:treat.duration.m), REML=T,data=Cum.abs_gamma, na.action="na.fail")
anova(m7.cum.abs2)
Anova(m7.cum.abs2)
r.squaredGLMM(m7.cum.abs2)

#Response: LRR_cum.abs_gamma
#Chisq Df Pr(>Chisq)    
#MAT_rescale       0.6125  1   0.433862    
#delta.P1_rescale 12.9552  1   0.000319 ***

m7.cum.abs3<-lmer(LRR_cum.abs_gamma~delta.P1_rescale+MAP_rescale+
                    (1|study:site:treat.duration.m), REML=F,data=Cum.abs_gamma, na.action="na.fail")
anova(m7.cum.abs3)
Anova(m7.cum.abs3)
r.squaredGLMM(m7.cum.abs3)

#Response: LRR_cum.abs_gamma
#Chisq Df Pr(>Chisq)    
#delta.P1_rescale 18.6710  1  1.553e-05 ***
#MAP_rescale       0.4514  1     0.5017    

##############
##### S ######
#############



####model selection

m5.beta_alt<-lmer(LRR_gammaS~MAT_rescale * delta.P1_rescale + MAT_rescale * I(delta.P1_rescale^2) +
                  MAP_rescale * delta.P1_rescale + MAP_rescale * I(delta.P1_rescale^2)+
                    (1|study:site:treat.duration.m),  REML=F,data=Div.data_gamma, na.action="na.fail")


##########---- multimodel inference----###########

d5.beta_alt<-dredge(m5.beta_alt)
print(d5.beta_alt)

d5.AIC_alt<-model.avg(d5.beta_alt, subset = delta < 2)
summary(d5.AIC_alt)


####only one model, no averaging possible!


m5.best<-lmer(LRR_gammaS~delta.P1_rescale+
                +(1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m5.best)
Anova(m5.best)
r.squaredGLMM(m5.best) 
qqnorm(resid(m5.best))
qqline(resid(m5.best))
hist(resid(m5.best))

#Response: LRR_gammaS
#Chisq Df Pr(>Chisq)  
#delta.P1_rescale 4.354  1    0.03692 *

m5.best2<-lmer(LRR_gammaS~MAP_rescale*delta.P1_rescale+
                 (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m5.best2)
Anova(m5.best2)
r.squaredGLMM(m5.best2) 
qqnorm(resid(m5.best2))
qqline(resid(m5.best2))
hist(resid(m5.best2))

#Response: LRR_gammaS
#Chisq Df Pr(>Chisq)  
#MAP_rescale                  1.3980  1    0.23707  
#delta.P1_rescale             5.3504  1    0.02072 *
#MAP_rescale:delta.P1_rescale 3.0850  1    0.07902 .


m5.best3<-lmer(LRR_gammaS~MAT_rescale+delta.P1_rescale+
                +(1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m5.best3)
Anova(m5.best3)
r.squaredGLMM(m5.best3) 
qqnorm(resid(m5.best3))
qqline(resid(m5.best3))
hist(resid(m5.best3))

#Response: LRR_gammaS
#Chisq Df Pr(>Chisq)  
#      1.8019  1     0.1795  
#delta.P1_rescale 2.8165  1     0.0933 .
  ---

m5.best4<-lmer(LRR_gammaS~MAT_rescale+
                 (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m5.best4)
Anova(m5.best4)
r.squaredGLMM(m5.best4) 
qqnorm(resid(m5.best4))
qqline(resid(m5.best4))
hist(resid(m5.best4))

#Response: LRR_gammaS
#Chisq Df Pr(>Chisq)  
#MAT_rescale 3.1903  1    0.07407 .


m5.best5<-lmer(LRR_gammaS~MAP_rescale+delta.P1_rescale+
                 (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m5.best5)
Anova(m5.best5)
r.squaredGLMM(m5.best5) 
qqnorm(resid(m5.best5))
qqline(resid(m5.best5))
hist(resid(m5.best5))


#Response: LRR_gammaS
#Chisq Df Pr(>Chisq)  
#MAP_rescale      1.3050  1    0.25330  
#delta.P1_rescale 5.0619  1    0.02446 *

m5.best5<-lmer(LRR_gammaS~MAP_rescale * I(delta.P1_rescale^2)+delta.P1_rescale+
                 (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m5.best5)
Anova(m5.best5)
r.squaredGLMM(m5.best5) 
qqnorm(resid(m5.best5))
qqline(resid(m5.best5))
hist(resid(m5.best5))


#Response: LRR_gammaS
#Chisq Df Pr(>Chisq)   
#MAP_rescale                       1.2826  1   0.257408   
#I(delta.P1_rescale^2)             0.0461  1   0.830002   
#delta.P1_rescale                  7.9784  1   0.004734 **
# MAP_rescale:I(delta.P1_rescale^2) 3.9059  1   0.048117 * 


####################
##### ENS PIE ######
####################


#####model selection

m6.beta_alt<-lmer(LRR_gammaSPie~MAT_rescale * delta.P1_rescale + MAT_rescale * I(delta.P1_rescale^2) +
                    MAP_rescale * delta.P1_rescale + MAP_rescale * I(delta.P1_rescale^2)+
                (1|study:site:treat.duration.m),REML=F, data=Div.data_gamma, na.action="na.fail")
qqnorm(resid(m6.beta_alt))
hist(resid(m6.beta_alt))

##########---- multimodel inference----###########

d6.beta_alt<-dredge(m6.beta_alt)
print(d6.beta_alt)

d6.AIC_alt<-model.avg(d6.beta_alt, subset = delta < 2)
summary(d6.AIC_alt)


###best models

m6.best<-lmer(LRR_gammaSPie~  I(delta.P1_rescale^2)+
                (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m6.best)
Anova(m6.best)
r.squaredGLMM(m6.best)
qqnorm(resid(m6.best))
qqline(resid(m6.best))
hist(resid(m6.best))

#Response: LRR_gammaSPie
#Chisq Df Pr(>Chisq)  
#I(delta.P1_rescale^2) 4.673  1    0.03064 *


m6.best2<-lmer(LRR_gammaSPie~ MAP_rescale * I(delta.P1_rescale^2)+
                 (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m6.best2)
Anova(m6.best2)
r.squaredGLMM(m6.best2) 
qqnorm(resid(m6.best2))
hist(resid(m6.best2))

#MAP_rescale                       0.3362  1    0.56200  
#I(delta.P1_rescale^2)             4.2960  1    0.03820 *
#MAP_rescale:I(delta.P1_rescale^2) 2.9845  1    0.08406 .

m6.best3<-lmer(LRR_gammaSPie~delta.P1_rescale + MAP_rescale * I(delta.P1_rescale^2)+
                 (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m6.best3)
Anova(m6.best3)
r.squaredGLMM(m6.best3) 
qqnorm(resid(m6.best3))
hist(resid(m6.best3))

#Response: LRR_gammaSPie
#Chisq Df Pr(>Chisq)  
#delta.P1_rescale                  2.2006  1    0.13795  
#MAP_rescale                       0.4816  1    0.48770  
#I(delta.P1_rescale^2)             4.4567  1    0.03477 *
#MAP_rescale:I(delta.P1_rescale^2) 4.5474  1    0.03297 *

m6.best4<-lmer(LRR_gammaSPie~delta.P1_rescale+I(delta.P1_rescale^2)+
                 (1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
anova(m6.best4)
Anova(m6.best4)
r.squaredGLMM(m6.best4) 
qqnorm(resid(m6.best4))
hist(resid(m6.best4))

#Response: LRR_gammaSPie
#Chisq Df Pr(>Chisq)  
#delta.P1_rescale      0.5566  1    0.45562  
#I(delta.P1_rescale^2) 4.7608  1    0.02912 *

