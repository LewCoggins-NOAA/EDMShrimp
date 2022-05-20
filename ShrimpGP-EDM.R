#code to look at simple EDM models of Brown Shrimp Data using GP-EDM

rm(list = ls())

setwd("C:\\Users\\lewis.coggins\\Documents\\EDM\\ShrimpEDM")

library("readxl");library("dplyr");library(GPEDM);library(ggplot2);library(tidyr)


#Read Data from Chen-han

dat<-read_xlsx("SEAMAPBrownSummer.xlsx")
dat=cbind(dat,rowMeans(dat[,-1]))
colnames(dat)[11]<-"mean"
#dat[,-1]=scale(dat[,-1],center=T)
zones=c("zone_11","zone_14","zone_15","zone_16","zone_17","zone_18","zone_19","zone_20","zone_21","mean")

i=1

inpdat=as.data.frame(cbind(dat$Time,dat[[zones[i]]]))
inpdat=as.data.frame(cbind(inpdat[,1],rep(zones[i],dim(inpdat)[1]),inpdat[,2]))
inpdat[,3]=as.numeric(inpdat[,3]);inpdat[,1]=as.numeric(inpdat[,1])
names(inpdat)<-c("Time","Zone","Index");inpdat
N=nrow(inpdat)


ShrimpZone=fitGP(data = inpdat, y = "Index", pop = "Zone", E=3, tau=1, scaling = "local", predictmethod = "loo")
summary(ShrimpZone)

#> Plotting out of sample results.
plot(ShrimpZone)
plot()

plot(inpdat$Time,ShrimpZone$inputs$y,col='red',type='l')
#lines(inpdat$Time,ShrimpZone$outsampresults$predmean)
lines(inpdat$Time,ShrimpZone$insampresults$predmean)










#Tanya's example of GP-Edm
#rm(list = ls())

#install.packages("devtools") #if required
#devtools::install_github("tanyalrogers/GPEDM")
library(GPEDM);library(ggplot2);library(tidyr)


data("thetalog2pop")

pA=subset(thetalog2pop,Population=="PopA")
pB=subset(thetalog2pop,Population=="PopB")
N=nrow(pA)
par(mfrow=c(2,2),mar=c(4,4,2,1))
plot(Abundance~Time,data=pA,type="l",main="PopA")
plot(Abundance~Time,data=pB,type="l",main="PopB")
plot(pA$Abundance[1:(N-1)],pA$Abundance[2:N],
     xlab="Abundance t",ylab="Abundance t+1",main="PopA")
plot(pB$Abundance[1:(N-1)],pB$Abundance[2:N],
     xlab="Abundance t",ylab="Abundance t+1",main="PopB")


tlogtest=fitGP(data = thetalog2pop, y = "Abundance", pop = "Population", E=1, tau=1, 
               scaling = "local", predictmethod = "loo")
summary(tlogtest)

#> Plotting out of sample results.
plot(tlogtest)



#> Plotting in sample results.
plot(tlogtest, plotinsamp = T)


ggplot(tlogtest$insampresults,aes(x=timestep,y=predmean)) +
  facet_wrap(pop~., scales = "free") +
  geom_line() + geom_ribbon(aes(ymin=predmean-predfsd,ymax=predmean+predfsd), alpha=0.4) +
  geom_point(aes(y=obs)) +
  theme_bw()


#Plot conditional responses
con=getconditionals(tlogtest)

#have to convert conditionals output to long format
#there may be a more concise way to do this

npreds=length(grep("_yMean",colnames(con)))
conlong1=gather(con[,1:(npreds+1)],x,xValue,2:(npreds+1))
conlong2=gather(con[,c(1,(npreds+2):(2*npreds+1))],ym,yMean,2:(npreds+1))
conlong3=gather(con[,c(1,(2*npreds+2):(3*npreds+1))],ys,ySD,2:(npreds+1))
conlong=cbind.data.frame(conlong1,yMean=conlong2$yMean,ySD=conlong3$ySD)
ggplot(conlong,aes(x=xValue,y=yMean)) +
  facet_grid(pop~x, scales = "free") +
  geom_line() + geom_ribbon(aes(ymin=yMean-ySD,ymax=yMean+ySD), alpha=0.4) +
  theme_bw()



predvars=tlogtest$inputs$x_names2;predvars
npreds=length(predvars)
lscales=tlogtest$pars[1:npreds];lscales
par(mar=c(4,4,1,1))
plot(factor(predvars),lscales,xlab="Predictor",ylab="Inverse length scale")

#sequential predictions (they should improve over time)
seqpred=predict(tlogtest,predictmethod = "sequential");seqpred
plot(seqpred)


pAtrain=pA[1:40,]
pAtest=pA[41:50,]
tlogtest2=fitGP(data = pAtrain, y = "Abundance", E=2, tau=1,
                newdata = pAtest)
plot(tlogtest2)

pAlags=makelags(pA, y = "Abundance", E=2, tau=1)
pAdata=cbind(pA,pAlags)
pAtrain=pAdata[1:40,]
pAtest=pAdata[41:50,]
tlogtest3=fitGP(data = pAtrain, y = "Abundance", x=colnames(pAlags),
                newdata = pAtest)
plot(tlogtest3)

lags1=makelags(thetalog2pop,y=c("Abundance"),pop="Population",time="Time",E=3,tau=1)
fore1=makelags(thetalog2pop,y=c("Abundance"),pop="Population",time="Time",E=3,tau=1,forecast = T)
data1=cbind(thetalog2pop, lags1)

tlogfore=fitGP(data = data1, y = "Abundance", x=c("Abundance_1","Abundance_2","Abundance_3"), 
               pop = "Population", time = "Time", scaling = "local", newdata = fore1)

ggplot(tlogfore$insampresults,aes(x=timestep,y=predmean)) +
  facet_wrap(pop~., scales = "free") +
  geom_line() + geom_ribbon(aes(ymin=predmean-predsd,ymax=predmean+predsd), alpha=0.4) +
  geom_point(aes(y=obs)) +
  geom_point(data=tlogfore$outsampresults, aes(y=predmean), color="red") +
  geom_errorbar(data=tlogfore$outsampresults,
                aes(ymin=predmean-predsd,ymax=predmean+predsd),color="red") +
  theme_bw()


