
rm(list = ls())

setwd("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp")

library("readxl");library("dplyr");library(GPEDM);library(ggplot2);library(tidyr);library(rEDM)


#Read Data from Chen-han

dat<-read_xlsx("SEAMAPBrownSummer.xlsx")
dat=cbind(dat,rowMeans(dat[,-1]))
colnames(dat)[11]<-"mean"
#dat[,-1]=scale(dat[,-1],center=T)
zones=c("zone_11","zone_14","zone_15","zone_16","zone_17","zone_18","zone_19","zone_20","zone_21","mean")

#fit BigZone
for(i in 1:9){
  
  #i=1
  
  inpdat=as.data.frame(cbind(dat$Time,dat[[zones[i]]]))
  inpdat=as.data.frame(cbind(inpdat[,1],rep(zones[i],dim(inpdat)[1]),inpdat[,2]))
  inpdat[,3]=as.numeric(inpdat[,3]);inpdat[,1]=as.numeric(inpdat[,1])
  names(inpdat)<-c("Time","Zone","Index");
  if(i==1) bigdat=inpdat else bigdat=rbind(bigdat,inpdat)
}

N=nrow(bigdat);N

BigZone=fitGP(data = bigdat, y = "Index", pop = "Zone", E=4, tau=1, scaling = "global", predictmethod = "loo")












lib<-c(1,20)
pred<-c(1,33)
TAU=-1

#BeginLoop
for(i in 1:9){
  
#i=1

inpdat=as.data.frame(cbind(dat$Time,dat[[zones[i]]]))
names(inpdat)<-c("Time",zones[i])


#determine optimum Embedding Dimension
# simpout <- simplex(inpdat, lib, pred, tau=TAU,E=1:8);simpout
# 
# #plot(as.numeric(simpout$rho) ~ as.numeric(simpout$E), type = "l", xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)",main=names(inpdat)[2])
# rho=as.numeric(simpout$rho)[1:5];rho
# OptEmbed=which(rho==max(rho));OptEmbed
# rho=EmbedDimension(dataFrame=inpdat,lib=lib,pred=pred,tau=TAU,columns=names(inpdat)[2],target=names(inpdat)[2],showPlot=F)$rho[1:4]#;rho
# OptEmbed=which(rho==max(rho));OptEmbed
# 
# SimpPreds=Simplex(dataFrame=inpdat,lib=lib,pred=pred,E=OptEmbed,tau=TAU,columns=names(inpdat)[2],target=names(inpdat)[2],showPlot=F)
# #text(1995,max(SimpPreds$Observations,na.rm=T),names(inpdat)[2])
# preds=SimpPreds %>% filter(Time %in% (1990:2019))
# plot(Observations~Time,data=preds,type='p',col='red',main=names(inpdat)[2], ylab='Index',ylim=c(0,max(preds$Observations,na.rm=T)*1.1))
# lines(Predictions~Time,data=preds,type='l',col='blue',lty=2)
# legend('topright',legend=c("Observed","Predicted"),lty=c(1,2),col=c('red','blue'),bty='n',)
# r=cor(preds$Observations, preds$Predictions)
# rho=ComputeError(preds$Observations, preds$Predictions)$rho
# text(1995,max(preds$Observations,na.rm=T),paste('rho=',round(rho,2)))
# text(1995,max(preds$Observations,na.rm=T)*.9,paste('E=',OptEmbed))
# text(1995,max(preds$Observations,na.rm=T)*.8,paste('Simplex'))
# 
# #Now Do SMAP
# e=OptEmbed
# fit=SMap(dataFrame=inpdat,E=e,lib=lib,pred=pred,columns=names(inpdat)[2],target=names(inpdat)[2],showPlot=F);fit
# preds=fit$predictions %>% filter(Time %in% (1990:2019))
# plot(Observations~Time,data=preds,type='l',col='red',main=names(inpdat)[2], ylab='Index',ylim=c(0,max(preds$Observations,na.rm=T)*1.1))
# lines(Predictions~Time,data=preds,type='l',col='blue',lty=2)
# legend('topright',legend=c("Observed","Predicted"),lty=c(1,2),col=c('red','blue'),bty='n',)
# r=cor(preds$Observations, preds$Predictions)
# rho=ComputeError(preds$Observations, preds$Predictions)$rho
# text(1995,max(preds$Observations,na.rm=T),paste('rho=',round(rho,2)))
# text(1995,max(preds$Observations,na.rm=T)*.9,paste('E=',e))
# text(1995,max(preds$Observations,na.rm=T)*.8,paste('SMAP'))



#Now Do GP-EDM
dat<-read_xlsx("SEAMAPBrownSummer.xlsx")
dat=cbind(dat,rowMeans(dat[,-1]))
colnames(dat)[11]<-"mean"
#dat[,-1]=scale(dat[,-1],center=T)
zones=c("zone_11","zone_14","zone_15","zone_16","zone_17","zone_18","zone_19","zone_20","zone_21","mean")
#e=c(3,7,4,4,7,7,7,7,6)#e=rep(7,9)
e=rep(7,9)


inpdat=as.data.frame(cbind(dat$Time,dat[[zones[i]]]))
inpdat=as.data.frame(cbind(inpdat[,1],rep(zones[i],dim(inpdat)[1]),inpdat[,2]))
inpdat[,3]=as.numeric(inpdat[,3]);inpdat[,1]=as.numeric(inpdat[,1])
names(inpdat)<-c("Time","Zone","Index");inpdat
N=nrow(inpdat)

ShrimpZone=fitGP(data = inpdat, y = "Index", pop = "Zone", E=e[i], tau=1, scaling = "global", predictmethod = "loo")

#> Plotting out of sample results.
#plot(ShrimpZone)

plot(inpdat$Time,ShrimpZone$inputs$y,col='black',type='p',ylab='Index',xlab='Year',main=zones[i])
lines(inpdat$Time,ShrimpZone$insampresults$predmean,col='red',lty=2)
lines(inpdat$Time,ShrimpZone$outsampresults$predmean,col='blue')
#legend('topright',legend=c("Observed","InSampPred","OutSampPred"),pch=c("o","",""),lty=c(0,1,2),col=c('black','red','blue'),bty='n')
legend('topright',legend=c("OutSampPred","InSampPred"),lty=c(1,2),col=c('blue','red'),bty='n')
rho=ComputeError(ShrimpZone$outsampresults$obs, ShrimpZone$outsampresults$predmean)$rho
text(1995,max(ShrimpZone$outsampresults$obs,na.rm=T),paste('rho=',round(rho,2)))
text(1995,max(ShrimpZone$outsampresults$obs,na.rm=T)*.94,paste('E=',e[i]))
text(1995,max(ShrimpZone$outsampresults$obs,na.rm=T)*.88,'GP-EDM')

#summary(ShrimpZone)

#getconditionals(ShrimpZone)
#seqpred=predict(ShrimpZone,predictmethod = "sequential")
#plot(seqpred)


plotdat=cbind(filter(BigZone$outsampresults,pop==zones[i]),inpdat$Time)
names(plotdat)[7]='Year'

plot(obs~Year,dat=plotdat,col='black',type='p',ylab='Index',xlab='Year',main=zones[i])
#lines(inpdat$Time,ShrimpZone$insampresults$predmean,col='red',lty=2)
lines(predmean~Year,dat=plotdat,col='blue')
#legend('topright',legend=c("Observed","InSampPred","OutSampPred"),pch=c("o","",""),lty=c(0,1,2),col=c('black','red','blue'),bty='n')
legend('topright',legend=c("Pred"),lty=c(1),col=c('blue'),bty='n')
rho=ComputeError(plotdat$obs, plotdat$predmean)$rho
text(1995,max(plotdat$obs,na.rm=T),paste('rho=',round(rho,2)))
text(1995,max(plotdat$obs,na.rm=T)*.9,"BigZone")


}




#Now Do GP-EDM On Adam's SEAMAPIndex

phimat=matrix(1,8,8)
rhovec=rep(1,8)

dat<-read_xlsx("SEAMAPSUMTotalUnscaled.xlsx")
#dat<-read_xlsx("SEAMAPSUMTotal.xlsx")
dat=cbind(dat,rep("Gulf",33))
dat=dat[,c(1,3,2)]
names(dat)[2]='extent'

for(i in 1:8){
SEAMAPFIT=fitGP(data =dat, y = "Index", pop = "extent", E=i, tau=1, scaling = "global", predictmethod = "loo")

if(i==1){plot(dat$Year,SEAMAPFIT$outsampresults$obs,ylab='Index',xlab='Year',pch=16)}
if(i==8){lines(dat$Year,SEAMAPFIT$outsampresults$predmean,lty=i, col=i,lwd=4)}
lines(dat$Year,SEAMAPFIT$outsampresults$predmean,lty=i, col=i)

rhovec[i]=ComputeError(SEAMAPFIT$outsampresults$obs, SEAMAPFIT$outsampresults$predmean)$rho
if(i==8){
  text(1995,max(SEAMAPFIT$outsampresults$obs,na.rm=T),paste('rho=',round(rhovec[i],2)))
  text(1995,max(SEAMAPFIT$outsampresults$obs,na.rm=T)*.95,paste('E=',i))
  }

getconditionals(SEAMAPFIT)
seqpred=predict(SEAMAPFIT,predictmethod = "sequential")
plot(seqpred)
phimat[i,1:i]=SEAMAPFIT$pars[1:i]
if(i<8){phimat[i,(i+1):8]=NA}

}

phimat
rhovec
matplot((phimat))
        

# for(i in 1:9){
#   
#   #i=1
#   
#   inpdat=as.data.frame(cbind(dat$Time,dat[[zones[i]]]))
#   inpdat=as.data.frame(cbind(inpdat[,1],rep(zones[i],dim(inpdat)[1]),inpdat[,2]))
#   inpdat[,3]=as.numeric(inpdat[,3]);inpdat[,1]=as.numeric(inpdat[,1])
#   names(inpdat)<-c("Time","Zone","Index");
#   if(i==1) bigdat=inpdat else bigdat=rbind(bigdat,inpdat)
# }
# 
# N=nrow(bigdat);N
# 
# BigZone=fitGP(data = bigdat, y = "Index", pop = "Zone", E=5, tau=1, scaling = "global", predictmethod = "loo")
# summary(BigZone)
# op=par()
# Par(mfcol=c(2,2))
# plot(BigZone)
# par(op)
# 




