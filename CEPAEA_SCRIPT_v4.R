########################################
### Behavioural syndromes and phenotype/environment dependence of behaviour in the polymorphic snail Cepaea nemoralis
### Maxime Dahirel, Valentin Gaudu and Armelle Ansart
### project concept: Maxime Dahirel and Armelle Ansart
### methodology: Maxime Dahirel and Armelle Ansart, with contributions from Valentin Gaudu
### data collection : Valentin Gaudu
### data analysis/ script author: Maxime Dahirel, from an initial contribution by Valentin Gaudu (partial     frequentist version of the univariate models)
### version 2; April 2019
### version 3: July 2019
### version 4: August 2019
#######################################

####### STEP 0: loading packages
library(coda)
library(rstan)
library(bayesplot)
library(brms)
library(plotrix)
#library(future)

library(reshape2)
library(plyr) # for adply
library(MASS) # for mvrnor
rstan_options(auto_write = TRUE)
options(mc.cores = 2) 
library(ggplot2)
library(tidyverse)
library(tidybayes)
library(modelr)
library(RColorBrewer)
library(matrixStats)
library(cowplot)
library(QGglmm)
library(patchwork)

N_chains <- 4
N_iter <- 20000 ## default number is 2000
N_warmup <- 10000 ## default number is 1000
N_thin <- 1 ## thinning frequency, in case it is needed to keep models at a manageable memory size ( 1: no thinning)


######## PREPROCESSING PART OF THE SCRIPT
tab <- data
tab <- subset(tab, is.na(lab_hibern) == FALSE) ## remove individuals added for dispersal tests but with no personality test
tab$id=factor(tab$id) ## and remove their levels to avoid problems down the line

### creer nouvelles variables
tab$boldness <- round(tab$boldness) # we are not more precise than 1sec
tab$boldness[tab$censored_possibleB == 1] <- 1200

##fpt
tab$speed[is.na(tab$speed) == TRUE] <- 0
tab$fpt=80/(tab$speed)    ###speed is recorded in mm.s-1
tab$fpt[tab$fpt>1200]=1200
tab$fpt <-round(tab$fpt)  ##we are not more precise than 1 sec

######## END OF MAJOR PREPROCESSING, can save dataset and start from here for publi
tab$censored_fpt=as.numeric(tab$fpt>=1200)
tab$censored_bold <- as.numeric(tab$boldness >= 1200)

tab$stemp <- scale(tab$temperature)
mean_temp <- attr(tab$stemp, "scaled:center")
sd_temp <- attr(tab$stemp, "scaled:scale")

tab$group <- interaction(tab$band_genotype, tab$landscape)
tab$order.tempS <- tab$order.temp - 2.5 ## center without scaling
tab$order.bS <- tab$order.b - 1.5 ## same
tab$hibern.S <- tab$lab_hibern - 1.5 # same

tab$landscape2 <- -0.5 + as.numeric(tab$landscape == "Marais") ## same

####add dummy values to be able to use the subset approach without restructuring the table and without
### setting NA flags (dummy values will be ignored during fitting)
tab$has.valid.bold <- as.numeric(is.na(tab$boldness) == FALSE) ### workaround with the subset() approach to be implemented
tab$censored_bold <- rep(tab$censored_bold[1:720], rep(1,720))
tab$boldness2 <- c(tab$boldness[1:720], rep(1200, 720)) ### workaround with the subset() approach to be implemented

##recoding
tab$is0b=(tab$band_genotype=="0b")-mean(tab$band_genotype=="0b")
tab$is3b=(tab$band_genotype=="3b")-mean(tab$band_genotype=="3b")
tab$is5b=(tab$band_genotype=="5b")-mean(tab$band_genotype=="5b")


##### MULTIVARIATE MODEL
####FORMULAS

bf_explo <- bf(fpt |cens(censored_fpt) ~ ((is3b+is5b) * landscape2) * stemp + order.tempS +
                    (stemp | p | box) + (stemp | q | id), sigma ~ 1, family = lognormal(link_sigma="identity"))

bf_boldness <- bf(boldness2 | subset(has.valid.bold) + cens(censored_bold) ~ (is3b+is5b) * landscape2 + order.bS +
                    (1 | p | box) + (1 | q | id), sigma ~ 1, family = lognormal(link_sigma="identity")) ##when sigma not modeled link is identitiy (half prior?)


### PRIOR
prior_multi <- c(
  set_prior("normal(0, 1)", class = "b",resp=c("boldness2","fpt")), ##### Fixed effects; weakly informative prior for response centred and scaled
  set_prior("normal(log(400),0.5)",class=c("Intercept"),resp=c("boldness2","fpt")),
  set_prior("lkj(2)", class = "cor"),
  set_prior("exponential(1)",class = "sd",resp=c("boldness2","fpt")),
  set_prior("exponential(1)",dpar = "sigma",resp=c("boldness2","fpt"),class="Intercept")
)

### fit models ##NB: the sampler may throw initialization errors at the beginning because it evaluates loglik at log(0)
### but then come out fine
###better prior for intercepts?

###FITTING
mod <- brm(mvbf(bf_explo + bf_boldness, rescor = FALSE),
              data = tab, chains=N_chains, iter = N_iter, warmup = N_warmup, thin = N_thin, prior = prior_multi, seed = 42,
              control=list(adapt_delta = 0.95,max_treedepth=15) #,save_ranef=FALSE
)

#careful, >1Gb with iter =20000 warmup = 10000 if saving everything
# use save_ranef=FALSE to masssively reduce size (if you don't need individual random effects, simply variances, for inference)


### relire vehtari paper on what is a good bulk and tail ess (i've seen 400, 100 by chain)
### this gives you the bulk and tail ess
### monitor(posterior_samples(mod,as.array=TRUE))


#### add estimation and comparison of R2c and R2m ?
######### saving model here


###extract relevant variance components

###Important: VP(obs) will necessarily be < VP(predicted) because censoring artificially limits VP(obs)

#GLOBAL INTERCEPT BETA0

beta0_bold <- (posterior_samples(mod_light,pars="b_boldness2_Intercept")[,1])
beta0_fpt <- (posterior_samples(mod_light,pars="b_fpt_Intercept")[,1])

###FIXED EFFECT VARIANCES

##total
VFfpt_latent<- rowVars(
  posterior_linpred(mod_light, re_formula = NA,scale = "linear",resp="fpt")
)

VFbold_latent<- rowVars(
  posterior_linpred(mod_light, re_formula = NA,scale = "linear",resp="boldness2")
)

##individual state only
newdata=tab %>% mutate(stemp=0,order.tempS=0,order.bS = 0)  #to get a design matrix with only individual state variance
VFstate_fpt_latent<- rowVars(
  posterior_linpred(mod_light, newdata=newdata, re_formula = NA,scale = "linear",resp="fpt")
)

VFstate_bold_latent<- rowVars(
  posterior_linpred(mod_light, newdata=newdata, re_formula = NA,scale = "linear",resp="boldness2")
)


### DISTRIBUTIONAL VARIANCES
VDfpt_latent<- (posterior_samples(mod_light,pars="sigma_fpt_Intercept")[,1])^2
VDbold_latent<- (posterior_samples(mod_light,pars="sigma_boldness2_Intercept")[,1])^2

##RANDOM EFFECTS INTERCEPT ONLY
##INDIVIDUAL LEVEL
VIfptintercept_latent <- (VarCorr(mod_light,summary=F)$id$sd[,"fpt_Intercept"])^2
VIboldintercept_latent <- (VarCorr(mod_light,summary=F)$id$sd[,"boldness2_Intercept"])^2
##BOX LEVEL
VBfptintercept_latent <- (VarCorr(mod_light,summary=F)$box$sd[,"fpt_Intercept"])^2
VBboldintercept_latent <- (VarCorr(mod_light,summary=F)$box$sd[,"boldness2_Intercept"])^2


##RANDOM EFFECTS TOTAL VARIANCE
##INDIVIDUAL LEVEL
###application of Johnson 2014 MEE (http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12225/abstract) to get the total random variance of one trait at one obs level (ie including temperature ranef)
##create a design matrix
design_id=as.matrix(expand.grid(Intercept=rep(1,length(unique(tab$id))),stemp=unique(tab$stemp)))
## extract the (posterior) VCV matrix
VCVfpt_id=VarCorr(mod_light,summary=FALSE)$id$cov[,c("fpt_Intercept","fpt_stemp"),c("fpt_Intercept","fpt_stemp")]
##apply equation (11) over the posterior
VIfpt_latent<-NA
for (i in 1:dim(VCVfpt_id)[1]){
  VIfpt_latent[i] <-   mean(diag(design_id %*% VCVfpt_id[i,,] %*% t(design_id)))  
}
VIbold_latent<-VIboldintercept_latent

##BOX LEVEL (see above for details)
design_box=as.matrix(expand.grid(Intercept=rep(1,length(unique(tab$box))),stemp=unique(tab$stemp)))
VCVfpt_box <- VarCorr(mod_light,summary=FALSE)$box$cov[,c("fpt_Intercept","fpt_stemp"),c("fpt_Intercept","fpt_stemp")]
##apply equation (11) over the posterior
VBfpt_latent<-NA
for (i in 1:dim(VCVfpt_box)[1]){
  VBfpt_latent[i] <-   mean(diag(design_box %*% VCVfpt_box[i,,] %*% t(design_box)))  
}
VBbold_latent<-VBboldintercept_latent


##posterior correlation boldness explo
cor_fpt_bold<-posterior_samples(mod_light,pars="cor_id__fpt_Intercept__boldness2_Intercept")[,1]

varcomps_fpt_latent=data.frame(beta0=beta0_fpt,
                        VF=VFfpt_latent,VFstate=VFstate_fpt_latent,
                        VI=VIfpt_latent,VIint=VIfptintercept_latent,
                        VB=VBfpt_latent,VBint=VBfptintercept_latent,
                        VD=VDfpt_latent,
                        trait="fpt",
                        type="latent"
                        ) 
varcomps_fpt_latent$VP<-rowSums(varcomps_fpt_latent[c("VF","VI","VB","VD")])

varcomps_bold_latent=data.frame(beta0=beta0_bold,
                               VF=VFbold_latent,VFstate=VFstate_bold_latent,
                               VI=VIbold_latent,VIint=VIboldintercept_latent,
                               VB=VBbold_latent,VBint=VBboldintercept_latent,
                               VD=VDbold_latent,
                               trait="boldness",
                               type="latent"
) 
varcomps_bold_latent$VP<-rowSums(varcomps_bold_latent[c("VF","VI","VB","VD")])



QGlognormal<-list(
  inv.link =function(x){exp(x)},
  var.func=function(x){0}, ##no supplementary distributional variance beyond the one in the log residuals (VD)
  d.inv.link=function(x){exp(x)} #derivative of exp is exp
)
##gives _exact_ observed scales and variances for 'simple' lognormal (no fixed effect, no random effect, only residual var)


beta0_obs=NA
VF_obs=NA
VFstate_obs=NA
VI_obs=NA
VIint_obs=NA
VB_obs=NA
VBint_obs=NA
VP_obs=NA
for(i in 1:dim(varcomps_fpt_latent)[1]){
  beta0_obs[i]<-QGicc(mu=varcomps_fpt_latent$beta0[i],var.comp=varcomps_fpt_latent$VF[i],var.p=varcomps_fpt_latent$VP[i],custom.model = QGlognormal,verbose=F)$mean.obs
  VF_obs[i]<-QGicc(mu=varcomps_fpt_latent$beta0[i],var.comp=varcomps_fpt_latent$VF[i],var.p=varcomps_fpt_latent$VP[i],custom.model = QGlognormal,verbose=F)$var.comp.obs
  VFstate_obs[i]<-QGicc(mu=varcomps_fpt_latent$beta0[i],var.comp=varcomps_fpt_latent$VFstate[i],var.p=varcomps_fpt_latent$VP[i],custom.model = QGlognormal,verbose=F)$var.comp.obs
  VI_obs[i]<-QGicc(mu=varcomps_fpt_latent$beta0[i],var.comp=varcomps_fpt_latent$VI[i],var.p=varcomps_fpt_latent$VP[i],custom.model = QGlognormal,verbose=F)$var.comp.obs
  VIint_obs[i]<-QGicc(mu=varcomps_fpt_latent$beta0[i],var.comp=varcomps_fpt_latent$VIint[i],var.p=varcomps_fpt_latent$VP[i],custom.model = QGlognormal,verbose=F)$var.comp.obs
  VB_obs[i]<-QGicc(mu=varcomps_fpt_latent$beta0[i],var.comp=varcomps_fpt_latent$VB[i],var.p=varcomps_fpt_latent$VP[i],custom.model = QGlognormal,verbose=F)$var.comp.obs
  VBint_obs[i]<-QGicc(mu=varcomps_fpt_latent$beta0[i],var.comp=varcomps_fpt_latent$VBint[i],var.p=varcomps_fpt_latent$VP[i],custom.model = QGlognormal,verbose=F)$var.comp.obs
  VP_obs[i]<-QGicc(mu=varcomps_fpt_latent$beta0[i],var.comp=varcomps_fpt_latent$VP[i],var.p=varcomps_fpt_latent$VP[i],custom.model = QGlognormal,verbose=F)$var.comp.obs
}
VD_obs=VP_obs-(VF_obs+VI_obs+VB_obs)

varcomps_fpt_obs=data.frame(beta0=beta0_obs,
                               VF=VF_obs,VFstate=VFstate_obs,
                               VI=VI_obs,VIint=VIint_obs,
                               VB=VB_obs,VBint=VBint_obs,
                               VD=VD_obs,
                               trait="fpt",
                               type="obs",
                            VP=VP_obs
) 

varcomps_fpt_obs$iter=1:dim(varcomps_fpt_obs)[1]
varcomps_fpt_latent$iter=1:dim(varcomps_fpt_latent)[1]
vcfpt=rbind(varcomps_fpt_obs,varcomps_fpt_latent)


beta0_obs=NA
VF_obs=NA
VFstate_obs=NA
VI_obs=NA
VIint_obs=NA
VB_obs=NA
VBint_obs=NA
VP_obs=NA
for(i in 1:dim(varcomps_bold_latent)[1]){
  beta0_obs[i]<-QGicc(mu=varcomps_bold_latent$beta0[i],var.comp=varcomps_bold_latent$VF[i],var.p=varcomps_bold_latent$VP[i],custom.model = QGlognormal,verbose=F)$mean.obs
  VF_obs[i]<-QGicc(mu=varcomps_bold_latent$beta0[i],var.comp=varcomps_bold_latent$VF[i],var.p=varcomps_bold_latent$VP[i],custom.model = QGlognormal,verbose=F)$var.comp.obs
  VFstate_obs[i]<-QGicc(mu=varcomps_bold_latent$beta0[i],var.comp=varcomps_bold_latent$VFstate[i],var.p=varcomps_bold_latent$VP[i],custom.model = QGlognormal,verbose=F)$var.comp.obs
  VI_obs[i]<-QGicc(mu=varcomps_bold_latent$beta0[i],var.comp=varcomps_bold_latent$VI[i],var.p=varcomps_bold_latent$VP[i],custom.model = QGlognormal,verbose=F)$var.comp.obs
  VIint_obs[i]<-QGicc(mu=varcomps_bold_latent$beta0[i],var.comp=varcomps_bold_latent$VIint[i],var.p=varcomps_bold_latent$VP[i],custom.model = QGlognormal,verbose=F)$var.comp.obs
  VB_obs[i]<-QGicc(mu=varcomps_bold_latent$beta0[i],var.comp=varcomps_bold_latent$VB[i],var.p=varcomps_bold_latent$VP[i],custom.model = QGlognormal,verbose=F)$var.comp.obs
  VBint_obs[i]<-QGicc(mu=varcomps_bold_latent$beta0[i],var.comp=varcomps_bold_latent$VBint[i],var.p=varcomps_bold_latent$VP[i],custom.model = QGlognormal,verbose=F)$var.comp.obs
  VP_obs[i]<-QGicc(mu=varcomps_bold_latent$beta0[i],var.comp=varcomps_bold_latent$VP[i],var.p=varcomps_bold_latent$VP[i],custom.model = QGlognormal,verbose=F)$var.comp.obs
}
VD_obs=VP_obs-(VF_obs+VI_obs+VB_obs)

varcomps_bold_obs=data.frame(beta0=beta0_obs,
                            VF=VF_obs,VFstate=VFstate_obs,
                            VI=VI_obs,VIint=VIint_obs,
                            VB=VB_obs,VBint=VBint_obs,
                            VD=VD_obs,
                            trait="boldness",
                            type="obs",
                            VP=VP_obs
) 

varcomps_bold_obs$iter=1:dim(varcomps_bold_obs)[1]
varcomps_bold_latent$iter=1:dim(varcomps_bold_latent)[1]
vcbold=rbind(varcomps_bold_obs,varcomps_bold_latent)

vc_complete=rbind(vcfpt,vcbold)

vc_complete %>%
  select(-beta0) %>% 
  gather(key="VC",value="value",VF,VFstate,VI,VIint,VB,VBint,VD,VP)%>% 
  select(-iter)%>% 
  group_by(trait,type,VC) %>% 
  mean_hdi() %>% 
  print(n=Inf)
#sae but for repeatabilities
vc_complete %>% mutate(Ri=VIint/VP,Rstate=(VIint+VFstate)/VP, ratio = VIint/(VIint+VFstate)) %>% 
  gather(key="RPT",value="value",Rstate,Ri,ratio) %>% 
  group_by(trait,type,RPT) %>% select(RPT,value) %>% mean_hdi()


vc_complete %>% 
  gather(key="which",value="VC",2:8) %>% 
  ggplot()+
  geom_violin(aes(x=which,y=VC/VP,fill=trait))+
  facet_wrap(~type)

#show relative variance components
vc_complete %>% gather(key="which",value="VC",2:8) %>% filter(type=="obs") %>% 
  mutate(which=fct_relevel(which,c("VF","VFstate","VB","VBint","VI","VIint","VD")) )%>% 
  mutate(trait=fct_recode(trait,Boldness="boldness",Exploration="fpt")) %>% 
  ggplot()+geom_halfeyeh(aes(y=reorder(which,desc(which)),x=100*VC/VP),point_interval = mean_hdi,.width = 0.001,scale="width")+facet_wrap(~trait)+
  scale_x_continuous(name="Variance explained (%, observed scale)", lim=c(0,100))+
  scale_y_discrete(name="Variance component",
                   labels=c(expression(italic(V[D])),
                            expression(italic(V[I(intercept)])),
                            expression(italic(V[I])),
                            expression(italic(V[B(intercept)])),
                            expression(italic(V[B])),
                            expression(italic(V[F(state)])),
                            expression(italic(V[F]))
                            )) +
  theme_cowplot()+
  background_grid()
  





####gives mean and HPDintervals for all "main" parameters
lat_fix=mod_light %>% 
  gather_draws(c(`b_.*`),regex=TRUE) %>% 
  mean_hdi() %>% 
  print(n=Inf)

lat_var=mod_light %>% gather_draws(c(`sd_.*`, `b_sigma_.*`),regex=TRUE)%>% mutate(.value=.value^2) %>% mean_hdi() %>% print(n=Inf)

lat_cor=mod_light %>% 
  gather_draws(c(`cor_.*`),regex=TRUE) %>% 
  mean_hdi() %>% 
  print(n=Inf)







### PLOTS
##### representing predictions on plots
## at average temperatures (everything else except phenotype is averaged out)

newdata <- unique(tab[c("band_genotype","is3b","is5b")]) %>% 
  mutate(stemp = 0, order.bS = 0, order.tempS = 0, landscape2 = 0, id = tab$id[1], box = tab$box[1],has.valid.bold=1) %>%
  add_fitted_draws(model = mod_light, re_formula = NA, scale = "linear",resp=c("boldness2")) %>%
  mutate(lboldness = .value) %>% 
  select(-(.value))

newdata2 <- unique(tab[c("band_genotype","is3b","is5b")]) %>% 
  mutate(stemp = 0, order.bS = 0, order.tempS = 0, landscape2 = 0, id = tab$id[1], box = tab$box[1],has.valid.bold=1) %>%
  add_fitted_draws(model = mod_light, re_formula = NA, scale = "linear",resp=c("fpt")) %>%
  mutate(lfpt = .value) %>% 
  select(-(.value))

newdata<-inner_join(newdata,newdata2) ; rm(newdata2)
#### to do: need to separate speed and boldness in two tables, add means, backtransform, mise en page
### for boldness
plot_1 <- filter(tab,is.na(boldness)==FALSE) %>% ggplot() +
  geom_line(aes(x = as.numeric(band_genotype)+0.8*order.bS, y=boldness, group=id),col = "grey", show.legend = FALSE,alpha=0.4) +
  geom_violin(data=newdata,aes(y = exp(lboldness), x = band_genotype, fill = band_genotype),draw_quantiles = 0.5, show.legend = FALSE,alpha=0.8) +
  scale_fill_brewer(palette = "YlOrBr", direction = 1) +
  scale_y_log10(name = "Boldness (latency, sec)", limits = c(3, 1200), breaks = c(5, 10, 20, 50, 100, 200, 500, 1000)) +
  scale_x_discrete(name = "Shell phenotype", labels = c("0 bands", "3 bands", "5 bands")) +
  #geom_line(data = tab, aes(x = band_genotype, y = boldness,group=id),position = position_dodge(width = 0.9), col = "grey", show.legend = FALSE,alpha=0.3)+
  theme_cowplot()


###arrow_bold
arrow_bold=ggplot()+geom_segment(aes(x=c(0,0),xend=c(0,0),y=c(30,500),yend=c(500,30)),arrow=arrow())+
  geom_text(aes(x=c(0,0),y=c(25,580),label=c("bolder","shyer")))+
  scale_y_log10(limits = c(3, 1200))+theme_void()

###arrow_speed
arrow_speed=ggplot()+geom_segment(aes(x=c(0,0),xend=c(0,0),y=c(400,1000),yend=c(1000,400)),arrow=arrow())+
  geom_text(aes(x=c(0,0),y=c(380,1050),label=c("faster","slower")))+
  scale_y_log10(limits = c(300, 1200))+theme_void()

### for speeds
plot_2 <- ggplot(data = tab) + 
  geom_line(aes(x = as.numeric(band_genotype)+0.3*stemp,fpt, group=id),col = "grey", show.legend = FALSE,alpha=0.4) +
  geom_violin(data=newdata, aes(y = exp(lfpt), x = band_genotype, fill = band_genotype),draw_quantiles = 0.5, show.legend = FALSE, alpha= 0.8) +
  scale_fill_brewer(palette = "YlOrBr", direction = 1) +
  scale_y_log10(name = "Exploration (first passage time, sec)", limits = c(250, 1200), breaks = c(250, 500, 750, 1000)) +
  scale_x_discrete(name = "Shell phenotype", labels = c("0 bands", "3 bands", "5 bands")) +
  #geom_line(data = tab, aes(x = band_genotype, y = fpt,group=id), position=position_dodge(width = 0.9),col = "grey", show.legend = FALSE,alpha=0.3)+
  theme_cowplot()



#########PLOT 3 
## at temperatures range (everything else except phenotype and temp is averaged out)
newdata <- unique(tab[c("band_genotype","is3b","is5b")])
newdata<-merge(newdata,expand.grid(stemp = seq_range(tab$stemp, n = 100), order.bS = 0, order.tempS = 0, landscape2 = 0, id = tab$id[1], box = tab$box[1]))
newdata <-newdata %>%
  add_fitted_draws(model = mod_light, re_formula = NA,scale = "linear",resp="fpt") %>%
  mutate(Temperature = (stemp * sd_temp) + mean_temp) %>%
  mutate(lfpt = .value) %>%
  group_by(.draw, Temperature) ###we average out phenotype here
#### to do: need to separate speed and boldness in two tables, add means, backtransform, mise en page
plot_3 <-  ggplot(data = tab, aes(x = temperature, y = fpt))+
  #geom_jitter(col = "grey", show.legend = FALSE,alpha=0.4) +
  geom_line(data = tab, aes(x = temperature, y = fpt,group=id), col = "grey", show.legend = FALSE,alpha=0.3) + 
  stat_lineribbon(data=newdata, aes(x = Temperature,y = exp(lfpt)), .width = 0.95,show.legend=FALSE,fill="pink", alpha= 0.8) +
  scale_y_log10(name = "", limits = c(250, 1200),breaks = c(250, 500, 750, 1000)) +
  scale_x_continuous(name = "Test temperature (°C)", limits = c(14, 26)) +
  theme_cowplot()

plot_grid(arrow_bold,plot_1, NULL,arrow_speed,plot_2, plot_3, labels = c("A","","", "B","", "C"), 
          ncol = 3, nrow = 2,rel_widths=c(0.2,1.5,1.5))

BLUPS=ranef(mod)
BLUPS_id=data.frame(bold=BLUPS$id[,"Estimate","boldness2_Intercept"],fpt=BLUPS$id[,"Estimate","fpt_Intercept"],
                    boldSD=BLUPS$id[,"Est.Error","boldness2_Intercept"],fptSD=BLUPS$id[,"Est.Error","fpt_Intercept"],
                    boldCIlow=BLUPS$id[,"Q2.5","boldness2_Intercept"],boldCIhigh=BLUPS$id[,"Q97.5","boldness2_Intercept"],
                    fptCIlow=BLUPS$id[,"Q2.5","fpt_Intercept"],fptCIhigh=BLUPS$id[,"Q97.5","fpt_Intercept"])

plot_4 <- ggplot(data=BLUPS_id,aes(x=bold,y=fpt))+
  geom_segment(aes(x=boldCIlow,xend=boldCIhigh,y=fpt,yend=fpt),col = "grey", show.legend = FALSE,alpha=0.5)+
  geom_segment(aes(x=bold,xend=bold,y=fptCIlow,yend=fptCIhigh),col = "grey", show.legend = FALSE,alpha=0.5)+
  geom_point(aes(x=bold,y=fpt))+
  scale_x_continuous(name = "Boldness BLUP")+
  scale_y_continuous(name = "Exploration BLUP")+
  theme_cowplot()

plot_cor<-ggplot(data.frame(cor=cor_fpt_bold))+
  geom_density(aes(x=cor),fill="grey",col="white")+
  theme_classic()+
  theme(axis.line.y.left=element_blank(),panel.border = element_rect(colour = "black",fill=NA))+
  scale_x_continuous(name="correlation",limits=c(-1,1),labels=c(-1,"",0,"",1))+
  scale_y_continuous(NULL,breaks=NULL)+
  geom_vline(xintercept = 0,lty=2)+
  theme_cowplot()



###arrow_bold
arrow_bold=ggplot()+geom_segment(aes(y=c(0,0),yend=c(0,0),x=c(-2,2.25),xend=c(2.25,-2)),arrow=arrow())+
  geom_text(aes(y=c(0,0),x=c(-2.5,2.75),label=c("bolder","shyer")))+
  scale_x_continuous(limits = c(-3, 3))+theme_void()

###arrow_speed
arrow_speed=ggplot()+geom_segment(aes(x=c(0,0),xend=c(0,0),y=c(-0.25,0.65),yend=c(0.65,-0.25)),arrow=arrow())+
  geom_text(aes(x=c(0,0),y=c(-0.35,0.75),label=c("faster","slower")))+
  scale_y_continuous(limits = c(-0.5,0.75))+theme_void()

plot_4b=ggdraw(plot_4)+draw_plot(plot_cor,x=0.75,y=0.1,height=0.2,width=0.2)



plot_grid(arrow_speed,plot_4b, NULL,arrow_bold, 
          ncol = 2, nrow = 2,rel_widths=c(0.15,1.5),rel_heights = c(1.5,0.15))


vc_complete %>% filter(type=="obs") %>% 
  mutate(VFstateVI=VP*VFstate/(VIint + VFstate)) %>% 
  gather(key="varcomp",value="var",c(2:8,13)) %>% 
  ggplot()+
  geom_violin(aes(x=varcomp,y=var/VP))+
  facet_wrap(~trait)+
  scale_y_continuous("Variation explained (in % of VP, on observed scale)",limits=c(0,1))+
  scale_x_discrete("Variance component")+coord_flip()

### à faire, plotter les correlations et repetabilites
## looic /kfoldic finaux

### r2m, r2c

## cleaning dataset and script

##corriger bugs structurels dan sles plots qui font planter R plot

##ajouter hibernation aux modèles?
###voir si les bonnes moyennes empiriques pour boldness sont mean(raw) ou exp(mean(log(raw))) (comportement de la grand mean tend à dire option 2)
### cepaea nemoralis mean speed 3.63 cm/ min at 20°C on foam (here on plastic though) so mean fpt should be 2 min if they simply ran straight and left
### but higher here

##.poser une compa entre effet de geno et effet de temp sur speed explo
## à chaque temp sur le gradient testé, le 5 bandes explorent à la même vitesse que les 0 bandes 2 degrés avat