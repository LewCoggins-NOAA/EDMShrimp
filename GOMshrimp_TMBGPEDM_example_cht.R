########################################################################################################
# TMB-based GPEDM example for Gulf of Mexico shrimp dynamics [Tsai et al.(In Review) ]
#
# Issues to fix: 1. Some difference appeared when comparing with the GPEDM package maintained by Tanya
#                2. optim(),nlminb(),rprop() all produce the same estimation using compiled TMB code
#                3. The difference may result from the data pre-scaling (or not) of observation y and thus var(y)
#                   for the boundary condition or the pre-scaling for each coordinate vectors
#
#
########################################################################################################
#setwd("~/Desktop")
setwd("C:\\Users\\lewis.coggins\\Documents\\EDM\\ShrimpEDM")
library(TMB) # TMB package is required and C code must be pre-compiled
source("TMB_hGPEDM_example_v1_cht.R") # maintained by Cheng-Han
dyn.load(dynlib("hgpEDMtmb_v1"))
#
setwd("/Users/admin/Desktop/SpatialEDM_shrimps")
setwd("C:\\Users\\lewis.coggins\\Documents\\EDM\\ShrimpEDM")
cpue.all <- read.csv("forLewExample_cpue_nineStatzones_brownshrimp_summer.csv") #Brown shrimps summer CPUE from SEAMAP
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
# Could also try nlminb() which produces the same answer
out2 <- nlminb(obj$par,obj$fn,obj$gr,upper=c(log(rep(10,ncol(X))),log(var(y)),log(var(y)),log(1)))
# Could also try rprop() algorithm which produces similar answer
library(CaDENCE)
out3 <- rprop(obj$par,obj$fn)
#
print(est1 <- round(exp(out$par),3)) # check MAP estimated by "L-BFGS-B"
print(est2 <- round(exp(out2$par),3)) # check MAP estimated by "nlminb"
print(est3 <- round(exp(out3$par),3)) # check MAP estimated by "rprop"
ests <- cbind(est1,est2,est3)
ests <- rbind(ests,c(var(y),var(y),var(y)))
colnames(ests)<- c("L-BFGS-B","nlminb","rprop")
rownames(ests) <- c("phi1","phi2","phi3","tau","g","rho","predetermined_Var(obs)")
ests
#
print(out$convergence) # check convergence
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


####################################################################################################
# Comparing TMB-GPEDM results with the GPEDM version on https://github.com/tanyalrogers/GPEDM
# install.packages("devtools")
# devtools::install_github("tanyalrogers/GPEDM")
inputs_forTanya <- data.frame(Time=rep(1:nrow(cpue.all),9), Population=paste0("Pop",rep(1:9,each=nrow(cpue.all))), Abundance=as.numeric(as.matrix(cpue.all)))
inputs_forTanya$Population <- as.factor(inputs_forTanya$Population)
#gpfit_test_01 <-  fitGP(data=inputs_forTanya, yd="Abundance", pop="Population", E=3, tau=1) # this E is ED-1 described above
gpfit_test_01 <-  fitGP(data=inputs_forTanya, y="Abundance", pop="Population", E=3, tau=1) # this E is ED-1 described above
summary(gpfit_test_01)
#gpfit_test_02 <-  fitGP(data=inputs_forTanya, yd="Abundance", pop="Population", E=3, tau=1, predictmethod="loo")
gpfit_test_02 <-  fitGP(data=inputs_forTanya, y="Abundance", pop="Population", E=3, tau=1, predictmethod="loo")
summary(gpfit_test_02)

plot(gpfit_test_01)
plot(gpfit_test_02)
