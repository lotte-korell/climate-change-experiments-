getwd()
setwd("C:/Users/korell/R/Meta.analysis/results")
dir()

library(ggplot2)
library(lme4)
library(gridExtra)
library(cowplot)
library(ggrepel)
library(dplyr)
library(effects)
library(car)

Div.data_local <- read.csv("Precip_data_final.csv" ,dec=".", sep=";",h=T)
Div.data_gamma <- read.csv("Precip_data_gamma_final.csv" ,dec=".", sep=";",h=T)
Cum.abs_local <- read.csv("ES.cum.abs_final.csv" ,dec=".", sep=";",h=T)
Cum.abs_gamma<- read.csv("ES.cum.abs.gamma_final.csv" ,dec=".", sep=";",h=T)

Div.data_local$gr_map <- cut(Div.data_local$MAP_chelsa, c(-Inf,675,1125,Inf), c("dry", "mes","wet"))
Div.data_gamma$gr_map <- cut(Div.data_gamma$MAP_chelsa, c(-Inf,675,1125,Inf), c("dry", "mes","wet"))
Cum.abs_local$gr_map <- cut(Cum.abs_local$MAP_chelsa, c(-Inf,675,1125,Inf), c("dry", "mes","wet"))

#-------------------------------------------------Figure 1---------------------------------------------------------------------

m1.deltaP<-lmer(LRR_S~delta.P.year+(1|study:site:SiteBlock), REML=T,data=Div.data_local, na.action="na.fail")


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
  labs(y = "LRR species richness", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-2.5, 2), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2))+
  geom_hline(yintercept = 0, linetype =2, size = 0.5)+
  geom_point(size=2.5, shape=21, aes(fill = MAP_chelsa), stroke=0.3, alpha=0.4)+ 
  geom_ribbon(data = mm_P_LRR_S, mapping = aes(x = delta.P.year, ymin = plo_P, ymax = phi_P),
              fill = "grey", alpha=0.4) +  
  geom_line(data = mm_P_LRR_S, mapping = aes(x = delta.P.year,LRR_S ),
            colour = "black", size = 0.5) +  
  annotate(geom="text", x=-30, y=2, label="Local",
           color="black", size=6,  fontface="bold")+
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 900)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="none",
        aspect.ratio = 1) 
deltaP_S 

legendP <- get_legend(deltaP_S)


#gamma Sn

m6a.deltaP<-lmer(LRR_gammaSn~delta.P.year +(1|study:site),REML=T,data=Div.data_gamma, na.action="na.fail")
mm_P_LRR_gammaSn<- expand.grid(delta.P.year = seq(min(Div.data_gamma$delta.P.year), max(Div.data_gamma$delta.P.year), by = .1),
                              LRR_gammaSn = 0)
# Create matrix of relevant effect sizes
mm_P_gammaSn <- model.matrix(terms(m6a.deltaP), mm_P_LRR_gammaSn)

# Calculate LRR based on the relevant effect sizes
mm_P_LRR_gammaSn$LRR_gammaSn<- mm_P_gammaSn  %*% fixef(m6a.deltaP)  
pvar.mm_P_LRR_gammaSn <- diag(mm_P_gammaSn  %*% tcrossprod(vcov(m6a.deltaP), mm_P_gammaSn ))

# calculate upper and lower CI (sqrt var = se * 1.96)
mm_P_LRR_gammaSn<- data.frame(mm_P_LRR_gammaSn, 
                             plo_P_gammaSn = mm_P_LRR_gammaSn$LRR_gammaSn - 1.96*sqrt(pvar.mm_P_LRR_gammaSn),
                             phi_P_gammaSn = mm_P_LRR_gammaSn$LRR_gammaSn + 1.96*sqrt(pvar.mm_P_LRR_gammaSn)) 


deltaP_gammaSn<- ggplot(Div.data_gamma, 
                       aes(x = delta.P.year, y = LRR_gammaSn)) +
  labs(y = "LRR species richness", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-2.5, 2), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2))+
  geom_hline(yintercept = 0, linetype =2, size=0.5)+
  geom_point(size=2.5, shape=21, aes(fill = MAP_chelsa), stroke=0.3, alpha=0.4)+ 
  geom_ribbon(data = mm_P_LRR_gammaSn, mapping = aes(x = delta.P.year, ymin = plo_P_gammaSn, ymax = phi_P_gammaSn),
              fill = "grey", alpha = 0.6) +  
  geom_line(data = mm_P_LRR_gammaSn, mapping = aes(x = delta.P.year,LRR_gammaSn ),
            colour = "black",  size = 0.5) +  
  annotate(geom="text", x=-40, y=2, label="Site",
           color="black", size=6,  fontface="bold")+
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 900)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="none",
        aspect.ratio = 1) 
deltaP_gammaSn



m2.deltaP<-lmer(LRR_SPie~delta.P.year+(1|study:site:SiteBlock),REML=T,data=Div.data_local, na.action="na.fail")

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
  labs(y = "LRR evenness", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-2.5, 2), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2))+
  geom_hline(yintercept = 0, linetype =2, size = 0.4)+
  geom_point(size=2.5, shape=21, aes(fill = MAP_chelsa), stroke=0.3 , alpha=0.5)+ 
  geom_ribbon(data = mm_P_LRR_SPie, mapping = aes(x = delta.P.year, ymin = plo_P_SPie, ymax = phi_P_SPie),
              fill = "grey", alpha = 0.6) +  
  geom_line(data = mm_P_LRR_SPie, mapping = aes(x =delta.P.year,LRR_SPie ),
            colour = "black",  size = 0.5) +  
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 900)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="none",
        aspect.ratio = 1) 

deltaP_SPie 



deltaP_gammaSPie<-ggplot(Div.data_gamma,  
                         aes(x = delta.P.year, y = LRR_gammaSPie)) +
  labs(y = "LRR evenness", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-2.5, 2), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2))+
  geom_hline(yintercept = 0, linetype =2, size = 0.4)+
  geom_point(size=2.5, shape=21, aes(fill = MAP_chelsa), stroke=0.3, alpha=0.5)+ 
   scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 900)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="none",
        aspect.ratio = 1) 
deltaP_gammaSPie

#save 9 x 13

plot_grid(deltaP_S+theme(legend.position="none"),deltaP_gammaSn+theme(legend.position="none"),legendP,
          deltaP_SPie+theme(legend.position="none"),deltaP_gammaSPie+theme(legend.position="none"),"",
          labels=c("a","b","","c","d",""),label_size = 18,ncol = 3)



#--------------------------------Fig 2------------------------------------


deltaP_betaSn<- ggplot(Div.data_local, 
                      aes(x = delta.P.year, y = LRR_betaSn)) +
  labs(y = "LRR species richness", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-2.5, 2.5), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2, 2.5))+
  geom_hline(yintercept = 0, linetype =2, size = 0.5)+
  geom_point(size=3, shape=21, aes(fill = MAP_chelsa), stroke=0.3, alpha=0.5)+ 
  #geom_ribbon(data = mm_LRR_betaSn, mapping = aes(x = mm_LRR_betaSn$delta.P.year, ymin = plo, ymax = phi),
  #fill = "grey", alpha = 0.6) +  
  #geom_line(data = mm_LRR_betaSn, mapping = aes(x = mm_LRR_betaSn$delta.P.year,mm_LRR_betaSn$LRR_betaSn ),
  #colour = "black", size = 0.7) +  
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 900)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="none",
        aspect.ratio = 1) 
deltaP_betaSn


#----------------------------LRR beta SPie-----------------------

m5.deltaP<-lmer(LRR_betaSPie~delta.P.year+(1|study:site:SiteBlock),REML=T,data=Div.data_local, na.action="na.fail")
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
  labs(y = "LRR evenness", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-2.5, 2.5), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2, 2.5))+
  geom_hline(yintercept = 0, linetype =2)+
  geom_point(size=3, shape=21, aes(fill = MAP_chelsa), stroke=0.3, alpha=0.5)+ 
  geom_ribbon(data = mm_P_LRR_betaSPie, mapping = aes(x = delta.P.year, ymin = plo_P_betaSPie, ymax = phi_P_betaSPie),
              fill = "grey", alpha = 0.6) +  
  geom_line(data = mm_P_LRR_betaSPie, mapping = aes(x = delta.P.year,LRR_betaSPie ),
            colour = "black",  size = 0.5) +  
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 900)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="right",
        aspect.ratio = 1) 
deltaP_betaSPie 


#save 8 x 10

plot_grid(deltaP_betaSn+theme(legend.position="none"),deltaP_betaSPie+theme(legend.position="none"),legendP,
          labels=c("a","b",""))


#---------------------------------------------------------Predictor effect plots---------------------------------------------------
Cum.abs_local$gr_map <- cut(Cum.abs_local$MAP_chelsa, c(-Inf,675,1125,Inf), c("dry", "mes","wet"))
Div.data_local$gr_map <- cut(Div.data_local$MAP_chelsa, c(-Inf,675,1125,Inf), c("dry", "mes","wet"))
Div.data_gamma$gr_map <- cut(Div.data_gamma$MAP_chelsa, c(-Inf,675,1125,Inf), c("dry", "mes","wet"))


#---------------------------------------------------Figure 3---------------------------------------------------------------------
S_mod<-lmer(LRR_S~delta.P.year*MAP_chelsa+(1|study:site:SiteBlock), REML=T,data=Div.data_local, na.action="na.fail")

ef_S<-effect("delta.P.year*MAP_chelsa",S_mod, xlevels=list(MAP_chelsa=c(450,900,1350), delta.P.year=c(-60,-25,0,25,50,75,100)))
ef_S_Df<-as.data.frame(ef_S)
ef_S_Df$gr_map <- cut(ef_S_Df$MAP_chelsa, c(-Inf,675,1125,Inf), c("dry", "mes","wet"))


ef_S_plot<-ggplot(ef_S_Df, aes(x=delta.P.year, y=fit)) +

  geom_point(Div.data_local, mapping=aes(x = delta.P.year, y = LRR_S),shape=21,size=3, alpha=0.3)+
  scale_y_continuous(limits = c(-2.5, 2), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2))+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  geom_hline(yintercept = 0, linetype =2)+
  labs(y = "LRR Species richness", x = "Precipitation manipulation [%]")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=factor(MAP_chelsa)),alpha = 0.6) +
  geom_line() +
  facet_grid(. ~ gr_map)+
  scale_fill_manual(values = c("yellow", "limegreen", "dodgerblue"), guide = FALSE)+
  theme_classic()+
  theme(legend.title = element_text(colour="black", size=10, face="bold" ),
        axis.text.x=element_text(size=14, color="black"), 
        axis.text.y=element_text(size=14, color="black"),
        axis.title.y = element_text(size=18, color="black"),
        axis.title.x = element_text(size=18))
ef_S_plot


gamma.Sn_mod<-lmer(LRR_gammaSn~delta.P.year*MAP_chelsa+(1|study:site), REML=T,data=Div.data_gamma, na.action="na.fail")

ef_gammaSn<-effect("delta.P.year*MAP_chelsa",gamma.Sn_mod, xlevels=list(MAP_chelsa=c(450,900,1350), delta.P.year=c(-60,-25,0,25,50,75,100)))
ef_gammaSn_Df<-as.data.frame(ef_gammaSn)
ef_gammaSn_Df$gr_map <- cut(ef_gammaSn_Df$MAP_chelsa, c(-Inf,675,1125,Inf), c("dry", "mes","wet"))


ef_gammaSn<-ggplot(ef_gammaSn_Df, aes(x=delta.P.year, y=fit)) +
  geom_point(Div.data_gamma, mapping=aes(x = delta.P.year, y = LRR_gammaSn), shape=21,size=3, alpha=0.3)+
  scale_y_continuous(limits = c(-2.5, 2), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2))+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  geom_hline(yintercept = 0, linetype =2)+
  labs(y = "LRR species richness", x = "Precipitation manipulation [%]")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=factor(MAP_chelsa)),alpha = 0.6) +
  geom_line() +
  facet_grid(. ~ gr_map)+
  scale_fill_manual(values = c("yellow", "limegreen", "dodgerblue"), guide = FALSE)+
  theme_classic()+
  theme(legend.title = element_text(colour="black", size=10, face="bold" ),
        axis.text.x=element_text(size=14, color="black"), 
        axis.text.y=element_text(size=14, color="black"),
        axis.title.y = element_text(size=18, color="black"),
        axis.title.x = element_text(size=18))
ef_gammaSn

#save 9 x 9

plot_grid(ef_S_plot,ef_gammaSn,
          labels=c("a","b"), nrow=2)




#-------------------------------------------------Figure S2---------------------------------------------------------------------
#width = 10, height= 8

m3.deltaP<-lmer(LRR_cum.abs~delta.P.year+(1|study:site:SiteBlock),REML=T,data=Cum.abs_local, na.action="na.fail")

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
  geom_point(size=3, shape=21, aes(fill = MAP_chelsa), stroke=0.3,  alpha=0.5)+ 
  
  geom_ribbon(data = mm_P_LRR_cum.abs, mapping = aes(x = delta.P.year, ymin = plo_P_cum.abs, ymax =  phi_P_cum.abs),
              fill = "grey", alpha = 0.6) +  
  geom_line(data = mm_P_LRR_cum.abs, mapping = aes(x = delta.P.year,LRR_cum.abs ),
            colour = "black",  size = 0.5) +  
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 900)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="none",
        aspect.ratio = 1) 
deltaP_cum.abs 


legendP <- get_legend(deltaP_cum.abs )


m8.deltaP<-lmer(LRR_cum.abs_gamma~delta.P.year+(1|study:site),REML=T,data=Cum.abs_gamma, na.action="na.fail")
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
  geom_point(size=3, shape=21, aes(fill = MAP_chelsa), stroke=0.3, alpha=0.5)+ 
  geom_ribbon(data = mm_P_LRR_cum.abs_gamma, mapping = aes(x = delta.P.year, ymin = plo_P_cum.abs_gamma, ymax = phi_P_cum.abs_gamma),
              fill = "grey", alpha = 0.6) +  
  geom_line(data = mm_P_LRR_cum.abs_gamma, mapping = aes(x = delta.P.year,LRR_cum.abs_gamma ),
            colour = "black", size = 0.5) +  
  scale_fill_gradient2(low = "yellow",mid="springgreen3", high = "blue4", midpoint = 900)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=16, colour = "black"), 
        axis.text.y=element_text(size=16, colour = "black"),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        #legend.position="none",
        aspect.ratio = 1) 

deltaP_gammaCum.abs

#save 8 x 10

plot_grid(deltaP_cum.abs+theme(legend.position="none"),deltaP_gammaCum.abs+theme(legend.position="none"),legendP,
          labels=c("a","b",""))
#-------------------------------------------------Fig S3-----------------------------------------------------------

Cum.abs_mod<-lmer(LRR_cum.abs~delta.P.year*MAP_chelsa+(1|study:site:SiteBlock), REML=T,data=Cum.abs_local, na.action="na.fail")

ef_Cum.abs<-effect("delta.P.year*MAP_chelsa",Cum.abs_mod, xlevels=list(MAP_chelsa=c(450,900,1350), delta.P.year=c(-60,-25,0,25,50,75,100)))
ef_Cum.abs_Df<-as.data.frame(ef_Cum.abs)
ef_Cum.abs_Df$gr_map <- cut(ef_Cum.abs_Df$MAP_chelsa, c(-Inf,675,1125,Inf), c("dry", "mes","wet"))


#height=4, width=11

ef_Cum.abs_plot<-ggplot(ef_Cum.abs_Df, aes(x=delta.P.year, y=fit)) +
  geom_point(Cum.abs_local, mapping=aes(x = delta.P.year, y = LRR_cum.abs),shape=21,size=3, alpha=0.3)+
  scale_y_continuous(limits = c(-4, 4), breaks=c(-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2, 2.5, 3, 3.5,4))+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  geom_hline(yintercept = 0, linetype =2)+
  labs(y = "LRR total abundance", x = "Precipitation manipulation [%]")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=factor(MAP_chelsa)),alpha = 0.6) +
  geom_line() +
  facet_grid(. ~ gr_map)+
  scale_fill_manual(values = c("yellow", "limegreen", "dodgerblue"), guide = FALSE)+
  theme_classic()+
  theme(legend.title = element_text(colour="black", size=10, face="bold" ),
        axis.text.x=element_text(size=14, color="black"), 
        axis.text.y=element_text(size=14, color="black"),
        axis.title.y = element_text(size=18, color="black"),
        axis.title.x = element_text(size=18))
ef_Cum.abs_plot

#save 6 x 8

#-------------------------------------------------Fig S4-----------------------------------------------------------


Div.data_local$gr_pet<- cut(Div.data_local$PET_chelsa, c(-Inf,843,1095,Inf), c("low", "med","high"))

SPie_mod<-lmer(LRR_SPie~delta.P.year*PET_chelsa+(1|study:site:SiteBlock), REML=T,data=Div.data_local, na.action="na.fail")

ef_SPie<-effect("delta.P.year*PET_chelsa",SPie_mod, xlevels=list(PET_chelsa=c(717,969,1221), delta.P.year=c(-60,-25,0,25,50,75,100)))
ef_SPie_Df<-as.data.frame(ef_SPie)
ef_SPie_Df$gr_pet <- cut(ef_SPie_Df$PET_chelsa, c(-Inf,843,1095,Inf), c("low", "med","high"))

#height=4, width=11

ef_SPie_plot<-ggplot(ef_SPie_Df, aes(x=delta.P.year, y=fit)) +
  geom_point(Div.data_local, mapping=aes(x = delta.P.year, y = LRR_SPie),shape=21,size=3, alpha=0.3)+
  scale_y_continuous(limits = c(-2.5, 2), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2))+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  geom_hline(yintercept = 0, linetype =2)+
  labs(y = "Evenness", x = "Precipitation manipulation [%]")+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill=factor(PET_chelsa)),alpha = 0.6) +
  geom_line() +
  facet_grid(. ~ gr_pet)+
  scale_fill_manual(values = c("blue4", "deeppink", "red"), guide = FALSE)+
  theme_classic()+
  theme(legend.title = element_text(colour="black", size=10, face="bold" ),
        axis.text.x=element_text(size=14, color="black"), 
        axis.text.y=element_text(size=14, color="black"),
        axis.title.y = element_text(size=18, color="black"),
        axis.title.x = element_text(size=18))
ef_SPie_plot

#save 6 x 8
#-----------------------------------------Fig  S5 -----------------------------------------------
setwd("C:/Users/korell/R/Meta.analysis/results")

life_hist<- read.csv("life_history.csv" ,dec=".", sep=";",h=T)
life_hist$study_site<-paste(life_hist$study, "-", life_hist$site)
Div.data_local$study_site<-paste(Div.data_local$study, "-", Div.data_local$site)

agr<-aggregate(Div.data_local, by=list(Div.data_local$study_site,Div.data_local$study),FUN=mean)#get mean plot scale species richness/ evenness
life_hist1<-merge(life_hist, agr[,-(3:23)], by.x="study_site",by.y="Group.1")
life_hist2<-distinct(life_hist1,site,.keep_all = TRUE)#only keep each data point once
life_hist3<-subset(life_hist2,!is.na(Perennials), !is.na(Annuals.biennials))# remove datasets that do not provide species identity 

mod<-glmer(monocarp~ PET_chelsa + (1|study), family = binomial, data=life_hist3)
coefs <- coef(mod)

PET <-data.frame(PET_chelsa = seq(from=min(life_hist$PET_chelsa), 
                                  to=max(life_hist$PET_chelsa), by=.8))

PET$prop <-predict(mod, PET,type= c("response"), re.form = NA)

ggplot(life_hist3, aes(x=PET_chelsa, y=monocarp)) +
  geom_point(size=3, shape=21, aes(fill = PET_chelsa), stroke=0.3, alpha=0.5)+ 
  geom_line(PET, mapping = aes(x = PET_chelsa, y = prop)) +
  labs(y = "Probability", x = "Potential evapotranspiration [kg a-1]")+
  scale_fill_gradient2(low = "blue4",mid="deeppink", high = "red", midpoint = 969)+
  theme_classic()+
  theme(legend.title = element_text(colour="black", size=10, face="bold" ),
        axis.text.x=element_text(size=14, color="black"), 
        axis.text.y=element_text(size=14, color="black"),
        axis.title.y = element_text(size=18, color="black"),
        axis.title.x = element_text(size=18))


#------------------------------------------Fig S7 ---------------------------------------------

legendP <- get_legend(deltaP_SPie_alt)

m7.deltaP<-lmer(LRR_gammaSPie~I(delta.P.year^2)+(1|study:site),REML=T,data=Div.data_gamma, na.action="na.fail")
mm_P_LRR_gammaSPie<- expand.grid(delta.P.year = seq(min(Div.data_gamma$delta.P.year), max(Div.data_gamma$delta.P.year), by = .1),
                                 LRR_gammaSPie = 0)
# Create matrix of relevant effect sizes
mm_P_gammaSPie <- model.matrix(terms(m7.deltaP), mm_P_LRR_gammaSPie)

# Calculate LRR based on the relevant effect sizes
mm_P_LRR_gammaSPie$LRR_gammaSPie<- mm_P_gammaSPie  %*% fixef(m7.deltaP)  
pvar.mm_P_LRR_gammaSPie <- diag(mm_P_gammaSPie  %*% tcrossprod(vcov(m7.deltaP), mm_P_gammaSPie ))

# calculate upper and lower CI (sqrt var = se * 1.96)
mm_P_LRR_gammaSPie<- data.frame(mm_P_LRR_gammaSPie, 
                                plo_P_gammaSPie = mm_P_LRR_gammaSPie$LRR_gammaSPie - 1.96*sqrt(pvar.mm_P_LRR_gammaSPie),
                                phi_P_gammaSPie = mm_P_LRR_gammaSPie$LRR_gammaSPie + 1.96*sqrt(pvar.mm_P_LRR_gammaSPie)) 

deltaP_gammaSPie_alt<-ggplot(Div.data_gamma,  
                             aes(x = delta.P.year, y = LRR_gammaSPie)) +
  labs(y = "LRR evenness", x = "Precipitation manipulation [%]")+
  scale_x_continuous(limits = c(-70, 110), breaks=c(-50,-25, 0, 25 ,50 ,75, 100))+
  scale_y_continuous(limits = c(-2.5, 2), breaks=c(-2.5,-2,-1.5,-1,-0.5, 0, 0.5,1, 1.5, 2))+
  geom_hline(yintercept = 0, linetype =2, size = 0.4)+
  geom_point(size=2.5, shape=21, aes(fill = MAP_chelsa), stroke=0.3, alpha=0.5)+ 
  geom_ribbon(data = mm_P_LRR_gammaSPie, mapping = aes(x = delta.P.year, ymin = plo_P_gammaSPie, ymax = phi_P_gammaSPie),
              fill = "grey", alpha = 0.6) +  
  geom_line(data = mm_P_LRR_gammaSPie, mapping = aes(x =delta.P.year,LRR_gammaSPie ),
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
deltaP_gammaSPie_alt


plot_grid(deltaP_SPie_alt+theme(legend.position="none"),deltaP_gammaSPie_alt+theme(legend.position="none"),legendP,
          labels=c("a","b",""))


#save 8 x 10

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fig S8 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

site_mean<-aggregate(Div.data_local, by=list(Div.data_local$gr,Div.data_local$study,Div.data_local$delta.P.year, Div.data_local$treatment.level), FUN=mean) ####aggregate to get mean  per site/community
site_mean1<-subset(site_mean, select = c(Group.1, Group.2, size.local, size.gamma, N.replicates, LRR_S, LRR_SPie, LRR_betaS, LRR_betaSn, LRR_betaSPie))
colnames(site_mean1)<-c("study", "site", "size.local", "size.gamma", "N.replicates", "LRR_S", "LRR_SPie","LRR_betaS", "LRR_betaSn", "LRR_betaSPie")

site_var<-aggregate(Div.data_local, by=list(Div.data_local$gr,Div.data_local$study, Div.data_local$delta.P.year, Div.data_local$treatment.level), FUN=var)####aggregate to get variance  per site/community
site_var1<-subset(site_var, select = c(LRR_S, LRR_SPie, LRR_betaS, LRR_betaSn, LRR_betaSPie))

#calculate inverse of variance
site_var$LRR_S<-1/site_var$LRR_S
site_var$LRR_SPie<-1/site_var$LRR_SPie
site_var$LRR_betaSn<-1/site_var$LRR_betaS
site_var$LRR_betaSPie<-1/site_var$LRR_betaSPie

colnames(site_var1)<-c("var_LRR_S", "var_LRR_SPie","var_LRR_betaS", "var_LRR_betaSn", "var_LRR_betaSPie")

site_div<-cbind(site_mean1, site_var1)
#site_div1<-site_div[-c(41:67), ]exclude Miranda 

##co-dependencies

F8a<-ggplot(site_div, aes(x=size.local, y=LRR_S, size = size.local))+
  labs(y = "LRR species richness [mean]", x = "Local-scale plot size  [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=14, colour = "black"), 
        axis.text.y=element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        #legend.position="none",
        aspect.ratio = 1) 
F8a


F8b<-ggplot(site_div, aes(x=size.local, y=var_LRR_S, size = size.local)) +
  labs(y = "LRR species richness [variance]", x = "Local-scale plot size [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=14, colour = "black"), 
        axis.text.y=element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        #legend.position="none",
        aspect.ratio = 1) 
F8b



F8c<-ggplot(site_div, aes(x=size.local, y=LRR_SPie, size = size.local)) +
  labs(y = "LRR evenness [mean]", x = "Local-scale plot size [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=14, colour = "black"), 
        axis.text.y=element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        #legend.position="none",
        aspect.ratio = 1) 
F8c


F8d<-ggplot(site_div, aes(x=size.local, y=var_LRR_SPie, size = size.local)) +
  labs(y = "LRR evenness [variance]", x = "Local-scale plot size [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=14, colour = "black"), 
        axis.text.y=element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        #legend.position="none",
        aspect.ratio = 1) 
F8d

#turnover

F8e<-ggplot(site_div, aes(x=size.local, y=LRR_betaSn, size = size.local))+
  labs(y = "LRR species richness [mean]", x = "Local-scale plot size  [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=14, colour = "black"), 
        axis.text.y=element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        #legend.position="none",
        aspect.ratio = 1) 
F8e



F8f<-ggplot(site_div, aes(x=size.local, y=var_LRR_betaSn, size = size.local)) +
  labs(y = "LRR species richness [variance]", x = "Local-scale plot size [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=14, colour = "black"), 
        axis.text.y=element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        #legend.position="none",
        aspect.ratio = 1) 
F8f



F8g<-ggplot(site_div, aes(x=size.local, y=LRR_betaSPie, size = size.local)) +
  labs(y = "LRR evenness [mean]", x = "Local-scale plot size [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=14, colour = "black"), 
        axis.text.y=element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        #legend.position="none",
        aspect.ratio = 1) 
F8g


F8h<-ggplot(site_div, aes(x=size.local, y=var_LRR_betaSPie, size = size.local)) +
  labs(y = "LRR evenness [variance]", x = "Local-scale plot size [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=14, colour = "black"), 
        axis.text.y=element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        #legend.position="none",
        aspect.ratio = 1) 
F8h



#~~~~~~~~site scale~~~~~~~ 
F8i<-ggplot(Div.data_gamma, aes(x=size.gamma, y=LRR_gammaSn, size = size.gamma)) +
  labs(y = "LRR species richness", x = "Site-scale plot size [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=14, colour = "black"), 
        axis.text.y=element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        #legend.position="none",
        aspect.ratio = 1) 
F8i

F8j<-ggplot(Div.data_gamma, aes(x=size.gamma, y=LRR_gammaSPie, size = size.gamma)) +
  labs(y = "LRR evenness", x = "Site-scale plot size [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=14, colour = "black"), 
        axis.text.y=element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        #legend.position="none",
        aspect.ratio = 1) 
F8j

legendF8_local <- get_legend(F8a)
legendF8_site <- get_legend(F8i)

plot_grid(F8a+theme(legend.position="none"),F8b+theme(legend.position="none"),F8c+theme(legend.position="none"),F8d+theme(legend.position="none"),
          F8e+theme(legend.position="none"),F8f+theme(legend.position="none"),F8g+theme(legend.position="none"),F8h+theme(legend.position="none"),
          F8i+theme(legend.position="none"),F8j+theme(legend.position="none"),legendF8_local,legendF8_site,
          labels=c("a","b","c","d","e","f","g", "h","i","j","", ""),label_size = 18,ncol = 4)



#co-dependence plot size (local),  plot size (site), nr. of replicate with MAP

#local scale 
F9a<-ggplot(Div.data_local, aes(x=size.local, y=MAP_chelsa, size = size.local)) +
  labs(y = "Mean annual precipitation[mm a-1] ", x = "Local-scale plot size  [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=10, colour = "black"), 
        axis.text.y=element_text(size=10, colour = "black"),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        #legend.position="none",
        aspect.ratio = 1) 
F9a


F9b<-ggplot(Div.data_gamma, aes(x=size.gamma, y=MAP_chelsa, size = size.gamma)) +
  labs(y = "Mean annual precipitation[mm a-1] ", x = "Site-scale plot size  [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=10, colour = "black"), 
        axis.text.y=element_text(size=10, colour = "black"),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        #legend.position="none",
        aspect.ratio = 1) 
F9b

F9c<-ggplot(Div.data_local, aes(x=N.replicates, y=MAP_chelsa, size = N.replicates)) +
  labs(y = "Mean annual precipitation[mm a-1] ", x = "Number of replicates")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=10, colour = "black"), 
        axis.text.y=element_text(size=10, colour = "black"),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        #legend.position="none",
        aspect.ratio = 1) 
F9c

#co-dependence plot size (local),  plot size (site), nr. of replicate with PET


F9d<-ggplot(Div.data_local, aes(x=size.local, y=PET_chelsa, size = size.local)) +
  labs(y = "Potential evapotranspiration [kg a-1]", x = "Local-scale plot size  [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=10, colour = "black"), 
        axis.text.y=element_text(size=10, colour = "black"),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        #legend.position="none",
        aspect.ratio = 1) 
F9d

F9e<-ggplot(Div.data_gamma, aes(x=size.gamma, y=PET_chelsa, size = size.gamma)) +
  labs(y = "Potential evapotranspiration [kg a-1]", x = "Site-scale plot size  [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=10, colour = "black"), 
        axis.text.y=element_text(size=10, colour = "black"),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        #legend.position="none",
        aspect.ratio = 1) 
F9e


F9f<-ggplot(Div.data_local, aes(x=N.replicates, y=PET_chelsa, size = N.replicates)) +
  labs(y = "Potential evapotranspiration [kg a-1] ", x = "Number of replicates")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=10, colour = "black"), 
        axis.text.y=element_text(size=10, colour = "black"),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        #legend.position="none",
        aspect.ratio = 1) 
F9f


#co-dependence plot size (local),  plot size (site), nr. of replicate with delta P

F9g<-ggplot(Div.data_local, aes(x=size.local, y=delta.P.year, size = size.local)) +
  labs(y = "Precipitation manipulation [%]", x = "Local-scale plot size  [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=10, colour = "black"), 
        axis.text.y=element_text(size=10, colour = "black"),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        #legend.position="none",
        aspect.ratio = 1) 
F9g

F9h<-ggplot(Div.data_gamma, aes(x=size.gamma, y=delta.P.year, size = size.gamma)) +
  labs(y = "Precipitation manipulation [%] ", x = "Site-scale plot size  [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=10, colour = "black"), 
        axis.text.y=element_text(size=10, colour = "black"),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        #legend.position="none",
        aspect.ratio = 1) 
F9h

F9i<-ggplot(Div.data_local, aes(x=N.replicates, y=delta.P.year, size = N.replicates)) +
  labs(y = "Precipitation manipulation [%]", x = "Number of replicates")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=10, colour = "black"), 
        axis.text.y=element_text(size=10, colour = "black"),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        #legend.position="none",
        aspect.ratio = 1) 
F9i

F9j<-ggplot(Div.data_local, aes(x=size.local, y=treat.duration.y, size = size.local)) +
  labs(y = "Duration [years]", x = "Local-scale plot size  [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=10, colour = "black"), 
        axis.text.y=element_text(size=10, colour = "black"),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        #legend.position="none",
        aspect.ratio = 1) 
F9j

F9k<-ggplot(Div.data_gamma, aes(x=size.gamma, y=treat.duration.y, size = size.gamma)) +
  labs(y = "Duration [years]", x = "Site-scale plot size  [m2]")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=10, colour = "black"), 
        axis.text.y=element_text(size=10, colour = "black"),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        #legend.position="none",
        aspect.ratio = 1) 
F9k

F9l<-ggplot(Div.data_local, aes(x=N.replicates, y=treat.duration.y, size = N.replicates)) +
  labs(y = "Duration [years]", x = "Number of replicates")+
  geom_point(alpha=0.5)+
  theme_classic()+
  theme(plot.background = element_rect(fill="white"), 
        axis.text.x=element_text(size=10, colour = "black"), 
        axis.text.y=element_text(size=10, colour = "black"),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        #legend.position="none",
        aspect.ratio = 1) 
F9l

legendF9_local <- get_legend(F9a)
legendF9_site <- get_legend(F9b)
legendF9_replicates <- get_legend(F9c)


plot_grid(F9a+theme(legend.position="none"),F9b+theme(legend.position="none"),F9c+theme(legend.position="none"),
          F9d+theme(legend.position="none"),F9e+theme(legend.position="none"),F9f+theme(legend.position="none"),
          F9g+theme(legend.position="none"),F9h+theme(legend.position="none"),F9i+theme(legend.position="none"),
          F9j+theme(legend.position="none"),F9k+theme(legend.position="none"),F9l+theme(legend.position="none"),
          legendF9_local,legendF9_site,legendF9_replicates,
          labels=c("a","b","c","d","e","f","g", "h","i","j","k","l","", "", ""),label_size = 14,ncol = 3)
