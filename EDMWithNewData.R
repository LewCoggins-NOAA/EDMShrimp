#GP-EDM models of Brown Shrimp Data using GPEDM

rm(list = ls())

setwd("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp")

library("readxl");library("dplyr");library(GPEDM);library(ggplot2);library(tidyr);library(rEDM)

#Summer SEAMAP
#Read Data

BS.Sumdat<-read.csv("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp\\ChenhanDat\\SummerCPUEbs_zonesCJFAS_remake_forLew.csv",
                 header=T)#;BS.Sumdat
BS.Sumdat=BS.Sumdat[,-dim(BS.Sumdat)[2]]
#BS.Sumdat[,-1]=scale(dat[,-1],center=T)
zones=c("zone11","zone14","zone15","zone16","zone17","zone18","zone19","zone20","zone21")

inpdat <- data.frame(Time=rep(1:nrow(BS.Sumdat),9), Population=zones, Abundance=as.numeric(as.matrix(BS.Sumdat)))
inpdat$Population <- as.factor(inpdat$Population)
BS.Sumfit_01 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=3, tau=1) # this E is ED-1 described above
summary(BS.Sumfit_01)
BS.Sumfit_02 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=3, tau=1, predictmethod="loo")
summary(BS.Sumfit_02)

#Fall SEAMAP
#Read Data

BS.Falldat<-read.csv("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp\\ChenhanDat\\FallCPUEbs_zonesCJFAS_remake_forLew.csv",
                 header=T)#;BS.Falldat
BS.Falldat=BS.Falldat[,-dim(BS.Falldat)[2]]
#Falldat[,-1]=scale(dat[,-1],center=T)

inpdat <- data.frame(Time=rep(1:nrow(BS.Falldat),9), Population=zones, Abundance=as.numeric(as.matrix(BS.Falldat)))
inpdat$Population <- as.factor(inpdat$Population)
BS.Fallfit_01 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=3, tau=1) # this E is ED-1 described above
summary(BS.Fallfit_01)
BS.Fallfit_02 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=3, tau=1, predictmethod="loo")
summary(BS.Fallfit_02)

#zonemean SEAMAP
#Read Data

BS.zonemeandat<-read.csv("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp\\ChenhanDat\\zonemeanCPUEbs_zonesCJFAS_remake_forLew.csv",
                  header=T)#;BS.zonemeandat
BS.zonemeandat=BS.zonemeandat[,-dim(BS.zonemeandat)[2]]
#zonemeandat[,-1]=scale(dat[,-1],center=T)

inpdat <- data.frame(Time=rep(1:nrow(BS.zonemeandat),9), Population=zones, Abundance=as.numeric(as.matrix(BS.zonemeandat)))
inpdat$Population <- as.factor(inpdat$Population)
BS.zonemeanfit_01 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=3, tau=1) # this E is ED-1 described above
summary(BS.zonemeanfit_01)
BS.zonemeanfit_02 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=3, tau=1, predictmethod="loo")
summary(BS.zonemeanfit_02)


#GP-EDM models of White Shrimp Data using GPEDM

rm(list = ls())

setwd("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp")

library("readxl");library("dplyr");library(GPEDM);library(ggplot2);library(tidyr);library(rEDM)

#Summer SEAMAP
#Read Data

WS.Sumdat<-read.csv("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp\\ChenhanDat\\SummerCPUEws_zonesCJFAS_remake_forLew.csv",
                    header=T)#;WS.Sumdat
WS.Sumdat=WS.Sumdat[,-dim(WS.Sumdat)[2]]
#WS.Sumdat[,-1]=scale(dat[,-1],center=T)
zones=c("zone11","zone14","zone15","zone16","zone17","zone18","zone19","zone20","zone21")

inpdat <- data.frame(Time=rep(1:nrow(WS.Sumdat),9), Population=zones, Abundance=as.numeric(as.matrix(WS.Sumdat)))
inpdat$Population <- as.factor(inpdat$Population)
WS.Sumfit_01 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=3, tau=1) # this E is ED-1 described above
summary(WS.Sumfit_01)
WS.Sumfit_02 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=3, tau=1, predictmethod="loo")
summary(WS.Sumfit_02)

#Fall SEAMAP
#Read Data

WS.Falldat<-read.csv("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp\\ChenhanDat\\FallCPUEws_zonesCJFAS_remake_forLew.csv",
                     header=T)#;WS.Falldat
WS.Falldat=WS.Falldat[,-dim(WS.Falldat)[2]]
#Falldat[,-1]=scale(dat[,-1],center=T)

inpdat <- data.frame(Time=rep(1:nrow(WS.Falldat),9), Population=zones, Abundance=as.numeric(as.matrix(WS.Falldat)))
inpdat$Population <- as.factor(inpdat$Population)
WS.Fallfit_01 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=3, tau=1) # this E is ED-1 described above
summary(WS.Fallfit_01)
WS.Fallfit_02 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=3, tau=1, predictmethod="loo")
summary(WS.Fallfit_02)

#zonemean SEAMAP
#Read Data

WS.zonemeandat<-read.csv("C:\\Users\\lewis.coggins\\Documents\\GitHub\\EDMShrimp\\ChenhanDat\\zonemeanCPUEws_zonesCJFAS_remake_forLew.csv",
                         header=T)#;WS.zonemeandat
WS.zonemeandat=WS.zonemeandat[,-dim(WS.zonemeandat)[2]]
#zonemeandat[,-1]=scale(dat[,-1],center=T)

inpdat <- data.frame(Time=rep(1:nrow(WS.zonemeandat),9), Population=zones, Abundance=as.numeric(as.matrix(WS.zonemeandat)))
inpdat$Population <- as.factor(inpdat$Population)
WS.zonemeanfit_01 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=3, tau=1) # this E is ED-1 described above
summary(WS.zonemeanfit_01)
WS.zonemeanfit_02 <-  fitGP(data=inpdat, y="Abundance", pop="Population", E=3, tau=1, predictmethod="loo")
summary(WS.zonemeanfit_02)




