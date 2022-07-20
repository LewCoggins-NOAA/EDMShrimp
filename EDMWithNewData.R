#GP-EDM models of Brown Shrimp Data using GPEDM

rm(list = ls())

setwd("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp")

library("readxl");library("dplyr");library(GPEDM);library(ggplot2);library(tidyr);library(rEDM)

library(TMB) # TMB package is required and C code must be pre-compiled
source("TMB_hGPEDM_example_v1_cht.R") # maintained by Cheng-Han
dyn.load(dynlib("hgpEDMtmb_v1"))


#set embedding dimension
emdim<-3

#Summer SEAMAP
#Read Data

BS.Sumdat<-read.csv("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp\\ChenhanDat\\SummerCPUEbs_zonesCJFAS_remake_forLew.csv",
                 header=T)#;BS.Sumdat
BS.Sumdat=BS.Sumdat[,-dim(BS.Sumdat)[2]]
#BS.Sumdat[,-1]=scale(dat[,-1],center=T)
zones=paste("zone",c(11,seq(14,21)),sep='')



#inpdat <- data.frame(Time=rep(1:nrow(BS.Sumdat),9), Population=zones, Abundance=as.numeric(as.matrix(BS.Sumdat)))
inpdat <- data.frame(Time=rep(1:dim(BS.Sumdat)[1],times=9), Population=rep(zones,each=dim(BS.Sumdat)[1]), Abundance=as.numeric(as.matrix(BS.Sumdat)))
inpdat$Population <- as.factor(inpdat$Population)
BS.Sumfit_01 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=emdim, tau=1) # this E is ED-1 described above
summary(BS.Sumfit_01)
plot(BS.Sumfit_01)
BS.Sumfit_02 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=emdim, tau=1, predictmethod="loo")
summary(BS.Sumfit_02)
plot(BS.Sumfit_02)

for(i in 1:9){plot(BS.Sumdat[,i],type='b',main=zones[i])}

#Fall SEAMAP
#Read Data

BS.Falldat<-read.csv("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp\\ChenhanDat\\FallCPUEbs_zonesCJFAS_remake_forLew.csv",
                 header=T)#;BS.Falldat
BS.Falldat=BS.Falldat[,-dim(BS.Falldat)[2]]
#Falldat[,-1]=scale(dat[,-1],center=T)

inpdat <- data.frame(Time=rep(1:dim(BS.Sumdat)[1],times=9), Population=rep(zones,each=33), Abundance=as.numeric(as.matrix(BS.Falldat)))
inpdat$Population <- as.factor(inpdat$Population)
BS.Fallfit_01 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=emdim, tau=1) # this E is ED-1 described above
summary(BS.Fallfit_01)
BS.Fallfit_02 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=emdim, tau=1, predictmethod="loo")
summary(BS.Fallfit_02)

#zonemean SEAMAP
#Read Data

BS.zonemeandat<-read.csv("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp\\ChenhanDat\\zonemeanCPUEbs_zonesCJFAS_remake_forLew.csv",
                  header=T)#;BS.zonemeandat
BS.zonemeandat=BS.zonemeandat[,-dim(BS.zonemeandat)[2]]
#zonemeandat[,-1]=scale(dat[,-1],center=T)

inpdat <- data.frame(Time=rep(1:dim(BS.Sumdat)[1],times=9), Population=rep(zones,each=33), Abundance=as.numeric(as.matrix(BS.zonemeandat)))
inpdat$Population <- as.factor(inpdat$Population)
BS.zonemeanfit_01 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=emdim, tau=1) # this E is ED-1 described above
summary(BS.zonemeanfit_01)
BS.zonemeanfit_02 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=emdim, tau=1, predictmethod="loo")
summary(BS.zonemeanfit_02)


#GP-EDM models of White Shrimp Data using GPEDM


#Summer SEAMAP
#Read Data

WS.Sumdat<-read.csv("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp\\ChenhanDat\\SummerCPUEws_zonesCJFAS_remake_forLew.csv",
                    header=T)#;WS.Sumdat
WS.Sumdat=WS.Sumdat[,-dim(WS.Sumdat)[2]]
#WS.Sumdat[,-1]=scale(dat[,-1],center=T)

inpdat <- data.frame(Time=rep(1:dim(BS.Sumdat)[1],times=9), Population=rep(zones,each=33), Abundance=as.numeric(as.matrix(WS.Sumdat)))
inpdat$Population <- as.factor(inpdat$Population)
WS.Sumfit_01 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=emdim, tau=1) # this E is ED-1 described above
summary(WS.Sumfit_01)
plot(WS.Sumfit_01)
WS.Sumfit_02 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=emdim, tau=1, predictmethod="loo")
summary(WS.Sumfit_02)
plot(WS.Sumfit_02)

#Fall SEAMAP
#Read Data

WS.Falldat<-read.csv("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp\\ChenhanDat\\FallCPUEws_zonesCJFAS_remake_forLew.csv",
                     header=T)#;WS.Falldat
WS.Falldat=WS.Falldat[,-dim(WS.Falldat)[2]]
#Falldat[,-1]=scale(dat[,-1],center=T)

inpdat <- data.frame(Time=rep(1:dim(BS.Sumdat)[1],times=9), Population=rep(zones,each=33), Abundance=as.numeric(as.matrix(WS.Falldat)))
inpdat$Population <- as.factor(inpdat$Population)
WS.Fallfit_01 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=emdim, tau=1) # this E is ED-1 described above
summary(WS.Fallfit_01)
WS.Fallfit_02 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=emdim, tau=1, predictmethod="loo")
summary(WS.Fallfit_02)

#zonemean SEAMAP
#Read Data

WS.zonemeandat<-read.csv("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp\\ChenhanDat\\zonemeanCPUEws_zonesCJFAS_remake_forLew.csv",
                         header=T)#;WS.zonemeandat
WS.zonemeandat=WS.zonemeandat[,-dim(WS.zonemeandat)[2]]
#zonemeandat[,-1]=scale(dat[,-1],center=T)

inpdat <- data.frame(Time=rep(1:dim(BS.Sumdat)[1],times=9), Population=rep(zones,each=33), Abundance=as.numeric(as.matrix(WS.zonemeandat)))
inpdat$Population <- as.factor(inpdat$Population)
WS.zonemeanfit_01 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=emdim, tau=1) # this E is ED-1 described above
summary(WS.zonemeanfit_01)
WS.zonemeanfit_02 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=emdim, tau=1, predictmethod="loo")
summary(WS.zonemeanfit_02)

#Chenghan Brown Shrimp


cpue.all <- read.csv("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp\\ChenhanDat\\FallCPUEbs_zonesCJFAS_remake_forLew.csv",
                     header=T) #Brown shrimps summer CPUE from SEAMAP
cpue.all<-cpue.all[,-dim(cpue.all)[2]]
str(cpue.all)
#
ED <- 4 # embedding dimension ED-1
ntot <- NROW(embed(cpue.all[,1],ED)) #data size of time-delay embedding
pred.ix <- 1:ntot # prediction and library index for time-delay embedding
lib.ix <- 1:ntot
stack.all <- stackpred.all <- rr <- c()
for(i in 1:ncol(cpue.all)){
  tt <- cpue.all[,i]
  ttemb <- embed(tt,ED)[lib.ix,]
  ttembpred <- embed(tt,ED)[pred.ix,]
  rr <- rbind(rr,apply(ttemb,2,function(x){max(x)-min(x)})) #scaling factor for each spatial block
  stack.all <- rbind(stack.all,ttemb)
  stackpred.all <- rbind(stackpred.all,ttembpred)
}
#
nsite <- ncol(cpue.all) # number of statistical zones
nrblock <- NROW(ttemb)  # data size for each zone
nrblock.pred <- NROW(pred.ix) # predicted data size foe each zone
#
y <- stack.all[,1] # there's a slight difference depending scaled y
X <- stack.all[,-1]
yy <- stackpred.all[,1]
XX <- stackpred.all[,-1]
rr.temp <- rr
rr <- rr[,-1] # scaling factors for inputs
#
# R-TMB and R built-in optimization procedure
# initialize data and parameter settings
data_tmb <-list("X"=as.matrix(X),"y"=as.vector(y),"nsite"=nsite,"nrblock"=nrblock,"rr"=rr,"ysigma2"=var(y)*1.01)
parms_tmb <- list("logW"=log(rep(1,ncol(X))),"logtau"=log(var(y)*0.1),"logg"=log(var(y)*0.1),"logrho"=log(0.5))
obj <- MakeADFun(data=data_tmb,parameters=parms_tmb,DLL="hgpEDMtmb_v1",hessian=TRUE)
out <- optim(obj$par,
             obj$fn,
             obj$gr,
             upper=c(log(rep(10,ncol(X))),log(var(y)),log(var(y)),log(1)),
             method="L-BFGS-B")

print(est1 <- round(exp(out$par),3)) 

Kc.est <- obj$report()$Kc # store useful outputs of covariance/ correlation matrix
R.est <- obj$report()$R 
print(round(R.est,3)) # dynamic correlation estimates
#
pred.gpTMB <- predstackGPsep(est=out$par,X=X,y=y,XX=XX,Kc=Kc.est,R=R.est,nsite=nsite,rr=rr,nrblock.pred=nrblock.pred,nrblock=nrblock)
paste0("In-sample prediction rho: ",round(cor(pred.gpTMB$mup,yy),3)) # in-sample prediction r
plot(pred.gpTMB$mup,yy,xlab="GP prediction",ylab="observation"); abline(0,1)
# comparing zone-scale predictions
predmat <- matrix(pred.gpTMB$mup,ncol=nsite)
yymat <- matrix(yy,ncol=nsite)
corZone <- NULL 
for(i in 1:nsite){
  corZone[i] <- cor(predmat[,i],yymat[,i])
}
paste0("In-sample zone-scale prediction rho: ",round(corZone,3))

for(i in 1:9){
    plot(yymat[,i],type='b',main=zones[i],col='blue')
    lines(predmat[,i],type='l')
}

