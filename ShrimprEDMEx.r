#code to look at simple EDM models of Brown Shrimp Data


rm(list = ls())

setwd("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp")

library("rEDM");library("readxl");library("dplyr")


#Read Data from Chen-han

dat<-read_xlsx("SEAMAPBrownSummer.xlsx")
dat=cbind(dat,rowMeans(dat[,-1]))
colnames(dat)[11]<-"mean"
#dat[,-1]=scale(dat[,-1],center=T)
zones=c("zone_11","zone_14","zone_15","zone_16","zone_17","zone_18","zone_19","zone_20","zone_21","mean")


#Analyze Chen-han zone data
lib<-c(1,20)
pred<-c(1,33)
TAU=-1

for(i in 1:9){
#i=1

inpdat=as.data.frame(cbind(dat$Time,dat[[zones[i]]]))
names(inpdat)<-c("Time",zones[i])

#determine optimum Embedding Dimension
simpout <- simplex(inpdat, lib, pred, tau=TAU,E=1:8);simpout

#plot(as.numeric(simpout$rho) ~ as.numeric(simpout$E), type = "l", xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)",main=names(inpdat)[2])
rho=as.numeric(simpout$rho)[1:5];rho
OptEmbed=which(rho==max(rho));OptEmbed
rho=EmbedDimension(dataFrame=inpdat,lib=lib,pred=pred,tau=TAU,columns=names(inpdat)[2],target=names(inpdat)[2],showPlot=F)$rho[1:4]#;rho
OptEmbed=which(rho==max(rho));OptEmbed

SimpPreds=Simplex(dataFrame=inpdat,lib=lib,pred=pred,E=OptEmbed,tau=TAU,columns=names(inpdat)[2],target=names(inpdat)[2],showPlot=F)
#text(1995,max(SimpPreds$Observations,na.rm=T),names(inpdat)[2])
preds=SimpPreds %>% filter(Time %in% (1990:2019))
plot(Observations~Time,data=preds,type='l',col='red',main=names(inpdat)[2], ylab='Index',ylim=c(0,max(preds$Observations,na.rm=T)*1.1))
lines(Predictions~Time,data=preds,type='l',col='blue',lty=2)
legend('topright',legend=c("Observed","Predicted"),lty=c(1,2),col=c('red','blue'),bty='n',)
r=cor(preds$Observations, preds$Predictions)
rho=ComputeError(preds$Observations, preds$Predictions)$rho
text(1995,max(preds$Observations,na.rm=T),paste('rho=',round(rho,2)))
text(1995,max(preds$Observations,na.rm=T)*.9,paste('E=',OptEmbed))
text(1995,max(preds$Observations,na.rm=T)*.8,paste('Simplex'))


e=OptEmbed
fit=SMap(dataFrame=inpdat,E=e,lib=lib,pred=pred,columns=names(inpdat)[2],target=names(inpdat)[2],showPlot=F);fit
preds=fit$predictions %>% filter(Time %in% (1990:2019))
plot(Observations~Time,data=preds,type='l',col='red',main=names(inpdat)[2], ylab='Index',ylim=c(0,max(preds$Observations,na.rm=T)*1.1))
lines(Predictions~Time,data=preds,type='l',col='blue',lty=2)
legend('topright',legend=c("Observed","Predicted"),lty=c(1,2),col=c('red','blue'),bty='n',)
r=cor(preds$Observations, preds$Predictions)
rho=ComputeError(preds$Observations, preds$Predictions)$rho
text(1995,max(preds$Observations,na.rm=T),paste('rho=',round(rho,2)))
text(1995,max(preds$Observations,na.rm=T)*.9,paste('E=',e))
text(1995,max(preds$Observations,na.rm=T)*.8,paste('SMAP'))
}



#repeat analysis on the Summer Brown Shrimp SEAMAP index
dat<-read_xlsx("SEAMAPSUMTotal.xlsx")

lib<-c(1,20)
pred<-c(1,33)
TAU=-1

simpout <- simplex(dat, lib, pred,tau=TAU);simpout
#plot(as.numeric(simpout$rho) ~ as.numeric(simpout$E), type = "l", xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)")
rho=as.numeric(simpout$rho)[1:5];rho
OptEmbed=which(rho==max(rho));OptEmbed
rho=EmbedDimension(dataFrame=dat,lib=lib,pred=pred,tau=TAU,columns="Index",target="Index",showPlot=F)$rho[1:4]#;rho
OptEmbed=which(rho==max(rho));OptEmbed

SimpPreds=Simplex(dataFrame=dat,lib=lib,pred=pred,E=OptEmbed,tau=TAU,columns="Index",target="Index",showPlot = )
preds=SimpPreds %>% filter(Year %in% (1989:2019))

plot(Observations~Year,data=preds,type='l',col='red',main='SEAMAP Summer', ylab='Index',ylim=c(0,max(preds$Observations,na.rm=T)*1.1))
lines(Predictions~Year,data=preds,type='l',col='blue',lty=2)
legend('topright',legend=c("Observed","Predicted"),lty=c(1,2),col=c('red','blue'),bty='n',)
r=cor(preds$Observations, preds$Predictions)
rho=ComputeError(preds$Observations, preds$Predictions)$rho
text(1995,max(preds$Observations,na.rm=T),paste('rho=',round(rho,4)))


