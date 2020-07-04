getwd()
setwd("C:/Users/korell/R/Meta.analysis/results")
dir()

library(ggplot2)
library(lme4)
library(effects)
library(grid)
library(gridExtra)
library(cowplot)
library(dplyr)


Div.data_local <- read.csv("Precip_data_final.csv" ,dec=".", sep=";",h=T)
Div.data_gamma <- read.csv("Precip_data_gamma_final.csv" ,dec=".", sep=";",h=T)
Cum.abs_local <- read.csv("ES.cum.abs_final.csv" ,dec=".", sep=";",h=T)
Cum.abs_gamma<- read.csv("ES.cum.abs.gamma_final.csv" ,dec=".", sep=";",h=T)



S_mod<-lmer(LRR_S~delta.P.year*MAP_bioclim+(1|study:site:SiteBlock:treat.duration.m), REML=T,data=Div.data_local, na.action="na.fail")
SPie_mod<-lmer(LRR_SPie~delta.P.year*MAP_bioclim+(1|study:site:SiteBlock:treat.duration.m), REML=T,data=Div.data_local, na.action="na.fail")
beta.S_mod<-lmer(LRR_betaS~delta.P.year*MAP_bioclim+(1|study:site:treat.duration.m), REML=T,data=Div.data_local, na.action="na.fail")
beta.SPie_mod<-lmer(LRR_betaSPie~delta.P.year*MAP_bioclim+(1|study:site:treat.duration.m), REML=T,data=Div.data_local, na.action="na.fail")
gamma.S_mod<-lmer(LRR_gammaS~delta.P.year*MAP_bioclim+(1|study:site:treat.duration.m), REML=T,data=Div.data_gamma, na.action="na.fail")
gamma.SPie_mod<-lmer(LRR_gammaSPie~delta.P.year*MAP_bioclim+(1|study:site:treat.duration.m), REML=T,data=Div.data_gamma, na.action="na.fail")


ef_S<-effect("delta.P.year*MAP_bioclim",S_mod, xlevels=list(MAP_bioclim=c(400,800, 1200), delta.P.year=c(-60,-25,0,25,50,75,100)))
ef_S_Df<-as.data.frame(ef_S)

ef_S<-ggplot(ef_S_Df, aes(x=delta.P.year, y=fit)) +
  labs(y = "Species richness", x = "Precipitation manipulation [%]")+
  scale_y_continuous(limits = c(-0.6, 0.6), breaks=c(-0.5,-0.25, 0, 0.25 ,0.50 ))+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=factor(MAP_bioclim)),alpha = 0.6) +
  geom_line() +
  facet_grid(. ~ MAP_bioclim)+
  scale_fill_manual(values = c("yellow", "limegreen", "dodgerblue"), guide = FALSE)+
  theme_classic()+
  theme(legend.title = element_text(colour="black", size=10, face="bold" ),
        axis.text.x=element_text(size=16, color="black"), 
        axis.text.y=element_text(size=16, color="black"),
        axis.title.y = element_text(size=18, color="black"),
        axis.title.x = element_text(size=18))
ef_S

png(filename="Y:/Home/korell/Postdoc/1_Projekte/Meta-analysis/Manuskript/Metaanalysis/Results/Figures/effects_S.png" , width = 11, height= 4, units="in", res=600)
ef_S
dev.off()


ef_S<-effect("delta.P.year*MAP_bioclim",S_mod, xlevels=list(MAP_bioclim=c(400,800), delta.P.year=c(-60,-25,0,25,50,75,100)))
ef_S_Df<-as.data.frame(ef_S)

ef_S<-ggplot(ef_S_Df, aes(x=delta.P.year, y=fit)) +
  labs(y = "Species richness", x = "Precipitation manipulation [%]")+
  scale_y_continuous(limits = c(-0.6, 0.6), breaks=c(-0.5,-0.25, 0, 0.25 ,0.50 ))+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=factor(MAP_bioclim)),alpha = 0.6) +
  geom_line() +
  facet_grid(. ~ MAP_bioclim)+
  scale_fill_manual(values = c("yellow", "dodgerblue"), guide = FALSE)+
  theme_classic()+
  theme(legend.title = element_text(colour="black", size=10, face="bold" ),
        axis.text.x=element_text(size=16, color="black"), 
        axis.text.y=element_text(size=16, color="black"),
        axis.title.y = element_text(size=18, color="black"),
        axis.title.x = element_text(size=18))
ef_S


ef_SPie<-effect("delta.P.year*MAP_bioclim",SPie_mod, xlevels=list(MAP_bioclim=c(400,800,1200), delta.P.year=c(-60,-25,0,25,50,75,100)))
ef_SPie_Df<-as.data.frame(ef_SPie)

ef_SPie<-ggplot(ef_SPie_Df, aes(x=delta.P.year, y=fit)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks=c(-0.5,-0.25, 0, 0.25 ,0.50 ))+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  labs(y = "Evenness", x = "Precipitation manipulation")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=factor(MAP_bioclim)),alpha = 0.6) +
  geom_line() +
  facet_grid(. ~ MAP_bioclim)+
  scale_fill_manual(values = c("yellow", "limegreen", "dodgerblue"), guide = FALSE)+
  theme_classic()+
  theme(legend.title = element_text(colour="black", size=10, face="bold" ),
        axis.text.x=element_text(size=16, color="black"), 
        axis.text.y=element_text(size=16, color="black"),
        axis.title.y = element_text(size=18, color="black"),
        axis.title.x = element_text(size=18))


png(filename="Y:/Home/korell/Postdoc/1_Projekte/Meta-analysis/Manuskript/Metaanalysis/Results/Figures/effects_SPie.png" , width = 11, height= 4, units="in", res=600)
ef_SPie
dev.off()

###beta scale 

ef_betaS<-effect("delta.P.year*MAP_bioclim",beta.S_mod, xlevels=list(MAP_bioclim=c(400,800,1200), delta.P.year=c(-60,-25,0,25,50,75,100)))
ef_betaS_Df<-as.data.frame(ef_betaS)

ef_betaS<-ggplot(ef_betaS_Df, aes(x=delta.P.year, y=fit)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks=c(-0.5,-0.25, 0, 0.25 ,0.50 ))+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-70,-25, 0, 25 ,50 ,75, 100))+
  labs(y = "Species richness", x = "Precipitation manipulation")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=factor(MAP_bioclim)),alpha = 0.6) +
  geom_line() +
  facet_grid(. ~ MAP_bioclim)+
  scale_fill_manual(values = c("yellow", "limegreen", "dodgerblue"), guide = FALSE)+
  theme_classic()+
  theme(legend.title = element_text(colour="black", size=10, face="bold" ),
        axis.text.x=element_text(size=16, color="black"), 
        axis.text.y=element_text(size=16, color="black"),
        axis.title.y = element_text(size=18, color="black"),
        axis.title.x = element_text(size=18))

ef_betaS

png(filename="Y:/Home/korell/Postdoc/1_Projekte/Meta-analysis/Manuskript/Metaanalysis/Results/Figures/effects_betaS.png" , width = 11, height= 4, units="in", res=600)
ef_betaS
dev.off()


ef_betaSPie<-effect("delta.P.year*MAP_bioclim",beta.SPie_mod, xlevels=list(MAP_bioclim=c(400,800,1200), delta.P.year=c(-60,-25,0,25,50,75,100)))
ef_betaSPie_Df<-as.data.frame(ef_betaSPie)

ef_betaSPie<-ggplot(ef_betaSPie_Df, aes(x=delta.P.year, y=fit)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks=c(-0.5,-0.25, 0, 0.25 ,0.50 ))+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  labs(y = "Evenness", x = "Precipitation manipulation [%]")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=factor(MAP_bioclim)),alpha = 0.5) +
  geom_line() +
  facet_grid(. ~ MAP_bioclim)+
  scale_fill_manual(values = c("yellow", "limegreen", "dodgerblue"), guide = FALSE)+
  theme_classic()+
  theme(legend.title = element_text(colour="black", size=10, face="bold" ),
        axis.text.x=element_text(size=16, color="black"), 
        axis.text.y=element_text(size=16, color="black"),
        axis.title.y = element_text(size=18, color="black"),
        axis.title.x = element_text(size=18))

ef_betaSPie

png(filename="Y:/Home/korell/Postdoc/1_Projekte/Meta-analysis/Manuskript/Metaanalysis/Results/Figures/effects_betaSPie.png" , width = 11, height= 4, units="in", res=600)
ef_betaSPie
dev.off()

###gamma scale

ef_gammaS<-effect("delta.P.year*MAP_bioclim",gamma.S_mod, xlevels=list(MAP_bioclim=c(400,800,1200), delta.P.year=c(-60,-25,0,25,50,75,100)))
ef_gammaS_Df<-as.data.frame(ef_gammaS)

ef_gammaS<-ggplot(ef_gammaS_Df, aes(x=delta.P.year, y=fit)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks=c(-0.5,-0.25, 0, 0.25 ,0.50 ))+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  labs(y = "Species richness", x = "Precipitation manipulation")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=factor(MAP_bioclim)),alpha = 0.5) +
  geom_line() +
  facet_grid(. ~ MAP_bioclim)+
  scale_fill_manual(values = c("yellow", "limegreen", "dodgerblue"), guide = FALSE)+
  theme_classic()+
  theme(legend.title = element_text(colour="black", size=10, face="bold" ),
        axis.text.x=element_text(size=16, color="black"), 
        axis.text.y=element_text(size=16, color="black"),
        axis.title.y = element_text(size=18, color="black"),
        axis.title.x = element_text(size=18))
ef_gammaS

png(filename="Y:/Home/korell/Postdoc/1_Projekte/Meta-analysis/Manuskript/Metaanalysis/Results/Figures/effects_gammaS.png" , width = 11, height= 4, units="in", res=600)
ef_gammaS
dev.off()


ef_gammaSPie<-effect("delta.P.year*MAP_bioclim",gamma.SPie_mod, xlevels=list(MAP_bioclim=c(400,800,1200), delta.P.year=c(-60,-25,0,25,50,75,100)))
ef_gammaSPie_Df<-as.data.frame(ef_gammaSPie)


ef_gammaSPie<-ggplot(ef_gammaSPie_Df, aes(x=delta.P.year, y=fit)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks=c(-0.5,-0.25, 0, 0.25 ,0.50 ))+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  labs(y = "Evenness", x = "Precipitation manipulation [%]")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=factor(MAP_bioclim)),alpha = 0.5) +
  geom_line() +
  facet_grid(. ~ MAP_bioclim)+
  scale_fill_manual(values = c("yellow", "springgreen3", "blue4"), guide = FALSE)+
  theme_classic()+
  theme(legend.title = element_text(colour="black", size=10, face="bold" ),
        axis.text.x=element_text(size=10, color="black"), 
        axis.text.y=element_text(size=10, color="black"),
        axis.title.y = element_text(size=10, color="black"),
        axis.title.x = element_text(size=10))


png(filename="Y:/Home/korell/Postdoc/1_Projekte/Meta-analysis/Manuskript/Metaanalysis/Results/Figures/effects_gammaSPie.png" , width = 11, height= 4, units="in", res=600)
ef_gammaSPie
dev.off()


