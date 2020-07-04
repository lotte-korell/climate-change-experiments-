getwd()
setwd("C:/Users/korell/R/Meta.analysis/results")
dir()

library(ggplot2)
library(lme4)
library(MuMIn)
library(gridExtra)
library(cowplot)
library(ggrepel)
library(dplyr)
library(effects)


Div.data_local <- read.csv("Precip_data_final.csv" ,dec=".", sep=";",h=T)
Div.data_gamma <- read.csv("Precip_data_gamma_final.csv" ,dec=".", sep=";",h=T)
Cum.abs_local <- read.csv("ES.cum.abs_final.csv" ,dec=".", sep=";",h=T)
Cum.abs_gamma<- read.csv("ES.cum.abs.gamma_final.csv" ,dec=".", sep=";",h=T)



#----------------------------Species richness-----------------------

m1.deltaP<-lmer(LRR_S~delta.P.year+(1|study:site:SiteBlock:treat.duration.m), REML=T,data=Div.data_local, na.action="na.fail")


# create a dataset with delta P and the sequence for which to compute value
mm_P_LRR_S<- expand.grid(delta.P.year = seq(min(Div.data_local$delta.P.year), max(Div.data_local$delta.P.year), by = .1),
                       LRR_S = 0)
# Create matrix of relevant effect sizes
mm_P_S <- model.matrix(terms(m1.deltaP), mm_P_LRR_S)

# Calculate LRR based on the relevant effect sizes standard error of coefficient
mm_P_LRR_S$LRR_S<- mm_P_S %*% fixef(m1.deltaP)  
pvar.mm_P_LRR_S <- diag(mm_P_S %*% tcrossprod(vcov(m1.deltaP), mm_P_S))

# calculate upper and lower CI (sqrt var = se * 1.96)
mm_P_LRR_S<- data.frame(mm_P_LRR_S, 
                      plo_P = mm_P_LRR_S$LRR_S - 1.96*sqrt(pvar.mm_P_LRR_S),
                      phi_P = mm_P_LRR_S$LRR_S + 1.96*sqrt(pvar.mm_P_LRR_S)) 



deltaP_S<-ggplot(Div.data_local,  
                 aes(x = delta.P.year, y = LRR_S)) +
  labs(y = "Species richness", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-2.5, 2), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2))+
  geom_hline(yintercept = 0, linetype =2, size = 0.5)+
  geom_point(size=2.5, shape=21, aes(fill = MAP_bioclim), stroke=0.3, alpha=0.4)+ 
  geom_ribbon(data = mm_P_LRR_S, mapping = aes(x = delta.P.year, ymin = plo_P, ymax = phi_P),
              fill = "grey", alpha=0.4) +  
  geom_line(data = mm_P_LRR_S, mapping = aes(x = delta.P.year,LRR_S ),
            colour = "black", size = 0.5) +  
  annotate(geom="text", x=-30, y=2, label="Local",
           color="black", size=6,  fontface="bold")+
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 800)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="none",
        aspect.ratio = 1) 
deltaP_S 




m6.deltaP<-lmer(LRR_gammaS~delta.P.year +(1|study:site:treat.duration.m),REML=T,data=Div.data_gamma, na.action="na.fail")
mm_P_LRR_gammaS<- expand.grid(delta.P.year = seq(min(Div.data_gamma$delta.P.year), max(Div.data_gamma$delta.P.year), by = .1),
                            LRR_gammaS = 0)
# Create matrix of relevant effect sizes
mm_P_gammaS <- model.matrix(terms(m6.deltaP), mm_P_LRR_gammaS)

# Calculate LRR based on the relevant effect sizes
mm_P_LRR_gammaS$LRR_gammaS<- mm_P_gammaS  %*% fixef(m6.deltaP)  
pvar.mm_P_LRR_gammaS <- diag(mm_P_gammaS  %*% tcrossprod(vcov(m6.deltaP), mm_P_gammaS ))

# calculate upper and lower CI (sqrt var = se * 1.96)
mm_P_LRR_gammaS<- data.frame(mm_P_LRR_gammaS, 
                           plo_P_gammaS = mm_P_LRR_gammaS$LRR_gammaS - 1.96*sqrt(pvar.mm_P_LRR_gammaS),
                           phi_P_gammaS = mm_P_LRR_gammaS$LRR_gammaS + 1.96*sqrt(pvar.mm_P_LRR_gammaS)) 


deltaP_gammaS<- ggplot(Div.data_gamma, 
                       aes(x = delta.P.year, y = LRR_gammaS)) +
  labs(y = "Species richness", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-2.5, 2), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2))+
  geom_hline(yintercept = 0, linetype =2, size=0.5)+
  geom_point(size=2.5, shape=21, aes(fill = MAP_bioclim), stroke=0.3, alpha=0.4)+ 
  geom_ribbon(data = mm_P_LRR_gammaS, mapping = aes(x = delta.P.year, ymin = plo_P_gammaS, ymax = phi_P_gammaS),
              fill = "grey", alpha = 0.6) +  
  geom_line(data = mm_P_LRR_gammaS, mapping = aes(x = delta.P.year,LRR_gammaS ),
            colour = "black",  size = 0.5) +  
  annotate(geom="text", x=-15, y=2, label="Regional",
           color="black", size=6,  fontface="bold")+
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 850)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="none",
        aspect.ratio = 1) 
deltaP_gammaS


legendP <- get_legend(deltaP_S)



#--------------------------------S Pie--------------------------------------


m2.deltaP<-lmer(LRR_SPie~delta.P.year+(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")

# create a dataset with delta P and the sequence for which to compute values  
mm_P_LRR_SPie<- expand.grid(delta.P.year = seq(min(Div.data_local$delta.P.year), max(Div.data_local$delta.P.year), by = .1),
                          LRR_SPie = 0)
# Create matrix of relevant effect sizes
mm_P_SPie <- model.matrix(terms(m2.deltaP), mm_P_LRR_SPie)

# Calculate LRR based on the relevant effect sizes
mm_P_LRR_SPie$LRR_SPie<- mm_P_SPie %*% fixef(m2.deltaP)  
pvar.mm_P_LRR_SPie <- diag(mm_P_SPie %*% tcrossprod(vcov(m2.deltaP), mm_P_SPie))

# calculate upper and lower CI (sqrt var = se * 1.96)
mm_P_LRR_SPie<- data.frame(mm_P_LRR_SPie, 
                         plo_P_SPie = mm_P_LRR_SPie$LRR_SPie - 1.96*sqrt(pvar.mm_P_LRR_SPie),
                         phi_P_SPie = mm_P_LRR_SPie$LRR_SPie + 1.96*sqrt(pvar.mm_P_LRR_SPie)) 

deltaP_SPie<-ggplot(Div.data_local,  
                    aes(x = delta.P.year, y = LRR_SPie)) +
  labs(y = "Evenness", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-2.5, 2), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2))+
  geom_hline(yintercept = 0, linetype =2, size = 0.4)+
  geom_point(size=2.5, shape=21, aes(fill = MAP_bioclim), stroke=0.3 , alpha=0.5)+ 
  geom_ribbon(data = mm_P_LRR_SPie, mapping = aes(x = delta.P.year, ymin = plo_P_SPie, ymax = phi_P_SPie),
              fill = "grey", alpha = 0.6) +  
  geom_line(data = mm_P_LRR_SPie, mapping = aes(x =delta.P.year,LRR_SPie ),
            colour = "black",  size = 0.5) +  
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 850)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="none",
        aspect.ratio = 1) 

deltaP_SPie 

#delta P was not significant

deltaP_gammaSPie<-ggplot(Div.data_gamma,  
                         aes(x = delta.P.year, y = LRR_gammaSPie)) +
  labs(y = "Evenness", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-2.5, 2), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2))+
  geom_hline(yintercept = 0, linetype =2, size = 0.4)+
  geom_point(size=2.5, shape=21, aes(fill = MAP_bioclim), stroke=0.3, alpha=0.5)+ 
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 850)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="none",
        aspect.ratio = 1) 
deltaP_gammaSPie


png(filename="Y:/Home/korell/Postdoc/1_Projekte/Meta-analysis/Manuskript/Metaanalysis/Results/Figures/Fig1.png" , width = 12, height= 8, units="in", res=600)

plot_grid(deltaP_S+theme(legend.position="none"),deltaP_gammaS+theme(legend.position="none"),legendP,
          deltaP_SPie+theme(legend.position="none"),deltaP_gammaSPie+theme(legend.position="none"),"",
          labels=c("a","b","","c","d",""),label_size = 18,ncol = 3)

dev.off()




#-----------------------------------cum abs-------------------------------------------------



m3.deltaP<-lmer(LRR_cum.abs~delta.P.year+(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Cum.abs_local, na.action="na.fail")

mm_P_LRR_cum.abs<- expand.grid(delta.P.year = seq(min(Cum.abs_local$delta.P.year), max(Cum.abs_local$delta.P.year), by = .1),
                             LRR_cum.abs = 0)
# Create matrix of relevant effect sizes
mm_P_cum.abs <- model.matrix(terms(m3.deltaP), mm_P_LRR_cum.abs)

# Calculate LRR based on the relevant effect sizes
mm_P_LRR_cum.abs$LRR_cum.abs<- mm_P_cum.abs %*% fixef(m3.deltaP)  
pvar.mm_P_LRR_cum.abs <- diag(mm_P_cum.abs %*% tcrossprod(vcov(m3.deltaP), mm_P_cum.abs))

# calculate upper and lower CI (sqrt var = se * 1.96)
mm_P_LRR_cum.abs<- data.frame(mm_P_LRR_cum.abs, 
                            plo_P_cum.abs = mm_P_LRR_cum.abs$LRR_cum.abs - 1.96*sqrt(pvar.mm_P_LRR_cum.abs),
                            phi_P_cum.abs = mm_P_LRR_cum.abs$LRR_cum.abs + 1.96*sqrt(pvar.mm_P_LRR_cum.abs)) 

deltaP_cum.abs<-ggplot(Cum.abs_local,  
                       aes(x = delta.P.year, y = LRR_cum.abs)) +
  labs(y = "Total abundance", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-4.1,4), breaks=c(-4,-3,-2,-1, 0, 1,2,3,4,5))+
  geom_hline(yintercept = 0, linetype =2, size = 0.5 )+
  geom_point(size=3, shape=21, aes(fill = MAP_bioclim), stroke=0.3,  alpha=0.5)+ 
  
  geom_ribbon(data = mm_P_LRR_cum.abs, mapping = aes(x = delta.P.year, ymin = plo_P_cum.abs, ymax =  phi_P_cum.abs),
              fill = "grey", alpha = 0.6) +  
  geom_line(data = mm_P_LRR_cum.abs, mapping = aes(x = delta.P.year,LRR_cum.abs ),
            colour = "black",  size = 0.5) +  
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 850)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="none",
        aspect.ratio = 1) 
deltaP_cum.abs 


m8.deltaP<-lmer(LRR_cum.abs_gamma~delta.P.year+(1|study:site:treat.duration.m),REML=T,data=Cum.abs_gamma, na.action="na.fail")
mm_P_LRR_cum.abs_gamma<- expand.grid(delta.P.year = seq(min(Cum.abs_gamma$delta.P.year), max(Cum.abs_gamma$delta.P.year), by = .1),
                                   LRR_cum.abs_gamma = 0)
# Create matrix of relevant effect sizes
mm_P_cum.abs_gamma <- model.matrix(terms(m8.deltaP), mm_P_LRR_cum.abs_gamma)

# Calculate LRR based on the relevant effect sizes
mm_P_LRR_cum.abs_gamma$LRR_cum.abs_gamma<- mm_P_cum.abs_gamma %*% fixef(m8.deltaP)  
pvar.mm_P_LRR_cum.abs_gamma <- diag(mm_P_cum.abs_gamma  %*% tcrossprod(vcov(m8.deltaP), mm_P_cum.abs_gamma))

# calculate upper and lower CI (sqrt var = se * 1.96)
mm_P_LRR_cum.abs_gamma<- data.frame(mm_P_LRR_cum.abs_gamma, 
                                  plo_P_cum.abs_gamma = mm_P_LRR_cum.abs_gamma$LRR_cum.abs_gamma - 1.96*sqrt(pvar.mm_P_LRR_cum.abs_gamma),
                                  phi_P_cum.abs_gamma = mm_P_LRR_cum.abs_gamma$LRR_cum.abs_gamma + 1.96*sqrt(pvar.mm_P_LRR_cum.abs_gamma)) 


deltaP_gammaCum.abs<-ggplot(Cum.abs_gamma,  
                            aes(x = delta.P.year, y = LRR_cum.abs_gamma)) +
  labs(y = "Total abundance", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-4.1,4), breaks=c(-4,-3,-2,-1, 0, 1,2,3,4,5))+
  geom_hline(yintercept = 0, linetype =2, size = 0.5)+
  geom_point(size=3, shape=21, aes(fill = MAP_bioclim), stroke=0.3, alpha=0.5)+ 
  geom_ribbon(data = mm_P_LRR_cum.abs_gamma, mapping = aes(x = delta.P.year, ymin = plo_P_cum.abs_gamma, ymax = phi_P_cum.abs_gamma),
              fill = "grey", alpha = 0.6) +  
  geom_line(data = mm_P_LRR_cum.abs_gamma, mapping = aes(x = delta.P.year,LRR_cum.abs_gamma ),
            colour = "black", size = 0.5) +  
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 850)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="none",
        aspect.ratio = 1) 

deltaP_gammaCum.abs

png(filename="Y:/Home/korell/Postdoc/1_Projekte/Meta-analysis/Manuskript/Metaanalysis/Results/Figures/FigS_cum.abs.png" , width = 10, height= 8, units="in", res=600)

plot_grid(deltaP_cum.abs+theme(legend.position="none"),deltaP_gammaCum.abs+theme(legend.position="none"),legendP,
          labels=c("a","b",""))

dev.off()


#---------------------------------beta turnover------------------------------------


deltaP_betaS<- ggplot(Div.data_local, 
                      aes(x = delta.P.year, y = LRR_betaS)) +
  labs(y = "Species richness", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-2.5, 2.5), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2, 2.5))+
  geom_hline(yintercept = 0, linetype =2, size = 0.5)+
  geom_point(size=3, shape=21, aes(fill = MAP_bioclim), stroke=0.3, alpha=0.5)+ 
  #geom_ribbon(data = mm_LRR_betaS, mapping = aes(x = mm_LRR_betaS$delta.P.year, ymin = plo, ymax = phi),
  #fill = "grey", alpha = 0.6) +  
  #geom_line(data = mm_LRR_betaS, mapping = aes(x = mm_LRR_betaS$delta.P.year,mm_LRR_betaS$LRR_betaS ),
  #colour = "black", size = 0.7) +  
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 850)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="none",
        aspect.ratio = 1) 
deltaP_betaS

#----------------------------LRR beta SPie-----------------------
m5.deltaP<-lmer(LRR_betaSPie~delta.P.year+(1|study:site:SiteBlock:treat.duration.m),REML=T,data=Div.data_local, na.action="na.fail")
mm_P_LRR_betaSPie<- expand.grid(delta.P.year = seq(min(Div.data_local$delta.P.year), max(Div.data_local$delta.P.year), by = .1),
                              LRR_betaSPie = 0)
# Create matrix of relevant effect sizes
mm_P_betaSPie <- model.matrix(terms(m5.deltaP), mm_P_LRR_betaSPie)

# Calculate LRR based on the relevant effect sizes
mm_P_LRR_betaSPie$LRR_betaSPie<- mm_P_betaSPie%*% fixef(m5.deltaP)  
pvar.mm_P_LRR_betaSPie <- diag(mm_P_betaSPie%*% tcrossprod(vcov(m5.deltaP), mm_P_betaSPie))

# calculate upper and lower CI (sqrt var = se * 1.96)
mm_P_LRR_betaSPie<- data.frame(mm_P_LRR_betaSPie, 
                             plo_P_betaSPie = mm_P_LRR_betaSPie$LRR_betaSPie - 1.96*sqrt(pvar.mm_P_LRR_betaSPie),
                             phi_P_betaSPie = mm_P_LRR_betaSPie$LRR_betaSPie + 1.96*sqrt(pvar.mm_P_LRR_betaSPie)) 

deltaP_betaSPie<-ggplot(Div.data_local,  
                        aes(x = delta.P.year, y = LRR_betaSPie)) +
  labs(y = "Evenness", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-2.5, 2.5), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2, 2.5))+
  geom_hline(yintercept = 0, linetype =2)+
  geom_point(size=3, shape=21, aes(fill = MAP_bioclim), stroke=0.3, alpha=0.5)+ 
  geom_ribbon(data = mm_P_LRR_betaSPie, mapping = aes(x = delta.P.year, ymin = plo_P_betaSPie, ymax = phi_P_betaSPie),
              fill = "grey", alpha = 0.6) +  
  geom_line(data = mm_P_LRR_betaSPie, mapping = aes(x = delta.P.year,LRR_betaSPie ),
            colour = "black",  size = 0.5) +  
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 850)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="right",
        aspect.ratio = 1) 
deltaP_betaSPie 


png(filename="Y:/Home/korell/Postdoc/1_Projekte/Meta-analysis/Manuskript/Metaanalysis/Results/Figures/Fig2_beta.png" , width = 10, height= 8, units="in", res=600)

plot_grid(deltaP_betaS+theme(legend.position="none"),deltaP_betaSPie+theme(legend.position="none"),legendP,
          labels=c("a","b",""))

dev.off()