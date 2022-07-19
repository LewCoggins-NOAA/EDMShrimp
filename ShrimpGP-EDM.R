#code to look at simple EDM models of Brown Shrimp Data using GP-EDM

rm(list = ls())

setwd("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp")

library("readxl");library("dplyr");library(GPEDM);library(ggplot2);library(tidyr);library(rEDM)


#Read Data from Chen-han

dat<-read_xlsx("SEAMAPBrownSummer.xlsx")
dat=cbind(dat,rowMeans(dat[,-1]))
colnames(dat)[11]<-"mean"
#dat[,-1]=scale(dat[,-1],center=T)
zones=c("zone_11","zone_14","zone_15","zone_16","zone_17","zone_18","zone_19","zone_20","zone_21","mean")
e=c(3,7,4,4,7,7,7,7,6)#e=rep(7,9)



for(i in 1:9){
  
#i=9

inpdat=as.data.frame(cbind(dat$Time,dat[[zones[i]]]))
inpdat=as.data.frame(cbind(inpdat[,1],rep(zones[i],dim(inpdat)[1]),inpdat[,2]))
inpdat[,3]=as.numeric(inpdat[,3]);inpdat[,1]=as.numeric(inpdat[,1])
names(inpdat)<-c("Time","Zone","Index");inpdat
N=nrow(inpdat)


ShrimpZone=fitGP(data = inpdat, y = "Index", pop = "Zone", E=e[i], tau=1, scaling = "global", predictmethod = "loo")
summary(ShrimpZone)

#> Plotting out of sample results.
#plot(ShrimpZone)

plot(inpdat$Time,ShrimpZone$inputs$y,col='black',type='p',ylab='Index',xlab='Year',main=zones[i])
lines(inpdat$Time,ShrimpZone$insampresults$predmean,col='red',lty=2)
lines(inpdat$Time,ShrimpZone$outsampresults$predmean,col='blue')
#legend('topright',legend=c("Observed","InSampPred","OutSampPred"),pch=c("o","",""),lty=c(0,1,2),col=c('black','red','blue'),bty='n')
legend('topright',legend=c("OutSampPred","InSampPred"),lty=c(1,2),col=c('blue','red'),bty='n')
rho=ComputeError(ShrimpZone$outsampresults$obs, ShrimpZone$outsampresults$predmean)$rho
text(1995,max(ShrimpZone$outsampresults$obs,na.rm=T),paste('rho=',round(rho,2)))
text(1995,max(ShrimpZone$outsampresults$obs,na.rm=T)*.95,paste('E=',e[i]))

#getconditionals(ShrimpZone)
#seqpred=predict(ShrimpZone,predictmethod = "sequential")
#plot(seqpred)
}


for(i in 1:9){
  
#i=1
  
  inpdat=as.data.frame(cbind(dat$Time,dat[[zones[i]]]))
  inpdat=as.data.frame(cbind(inpdat[,1],rep(zones[i],dim(inpdat)[1]),inpdat[,2]))
  inpdat[,3]=as.numeric(inpdat[,3]);inpdat[,1]=as.numeric(inpdat[,1])
  names(inpdat)<-c("Time","Zone","Index");
  if(i==1) bigdat=inpdat else bigdat=rbind(bigdat,inpdat)
}
  
N=nrow(bigdat);N
  
BigZone=fitGP(data = bigdat, y = "Index", pop = "Zone", E=9, tau=1, scaling = "global", predictmethod = "loo")
summary(BigZone)
plot(BigZone)


i=1
plotdat=cbind(filter(BigZone$outsampresults,pop==zones[i]),inpdat$Time)
names(plotdat)[7]='Year'

plot(obs~Year,dat=plotdat,col='black',type='p',ylab='Index',xlab='Year',main=zones[i])
#lines(inpdat$Time,ShrimpZone$insampresults$predmean,col='red',lty=2)
lines(predmean~Year,dat=plotdat,col='blue')
#legend('topright',legend=c("Observed","InSampPred","OutSampPred"),pch=c("o","",""),lty=c(0,1,2),col=c('black','red','blue'),bty='n')
legend('topright',legend=c("Pred"),lty=c(1),col=c('blue'),bty='n')
rho=ComputeError(plotdat$obs, plotdat$predmean)$rho
text(1995,max(plotdat$obs,na.rm=T),paste('rho=',round(rho,2)))
#text(1995,max(plotdat$obs,na.rm=T)*.95,paste('E=',e[i]))




