#This is the program to calculate the age effec and the relationships between carbon and oxygen istopes in Tree-rings#.
#Author: Guobao Xu
## Date: 2019-3-2, updated in 2019-8
## xgb234@lzb.ac.cn


## part 1. intial load the packages------------
library(treeclim)
library(RColorBrewer)
library(ggplot2)
library(readxl)
library(grid) 
library(corrplot)
library(reshape)
library(dplyr)
library(plyr)
library(Hmisc)## for the rcorr function
library(gtools)
library(psych)
library(boot)##  for the bootstap analysis
library(pgirmess) ## for the multple comparison, base on non-parameter
library(multcomp)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(dplR)
library(modifiedmk)
library(TTR)
library(trend)

windowsFonts(TN = windowsFont("Times New Roman"))
print ("============================================")

#-------part 2 load the carbon and oxygen isotopoe data and one-way ANOVA and plot----------

# 2.1 load and interpolate the climate data -----------------------------------------------
 # if you run it in your computer, please change it to your workspace
setwd("E:/Rwork/WSS-AGE")## if you run it in your computer, please change it to your workspace
rm(data.frame=ls())

MTemp<-read_excel("Zhaosu_station.xls",sheet="mean temp")
MaxTemp<-read_excel("Zhaosu_station.xls",sheet="mean max temp")
MinTemp<-read_excel("Zhaosu_station.xls",sheet="mean nin temp")
Pre<-read_excel("Zhaosu_station.xls",sheet="precipitaion")
Rh<-read_excel("Zhaosu_station.xls",sheet="RH")
# herein, we used the nearest data to interpolate the NA value, using the rpackages  DMwR
library(DMwR)
MTemp<-knnImputation(MTemp, k = 10, scale = T, meth = "weighAvg", distData = NULL)
MaxTemp<-knnImputation(MaxTemp, k = 10, scale = T, meth = "weighAvg", distData = NULL)
MinTemp<-knnImputation(MinTemp, k = 10, scale = T, meth = "weighAvg", distData = NULL)
Pre<-knnImputation(Pre, k = 10, scale = T, meth = "weighAvg", distData = NULL)
Rh<-knnImputation(Rh, k = 10, scale = T, meth = "weighAvg", distData = NULL)

### calculate seasonal climate data-
MTemp5.9<-apply(MTemp[c(6:10)],1,mean)/10
Pre5.9<-apply(Pre[c(6:10)],1,sum)/10
Rh5.9<-apply(Rh[c(6:10)],1,mean)

MTemp5.8<-apply(MTemp[c(6:9)],1,mean)/10
MaxTemp5.8<-apply(MaxTemp[c(6:9)],1,mean)/10
MinTemp5.8<-apply(MinTemp[c(6:9)],1,mean)/10

Pre5.8<-apply(Pre[c(6:9)],1,sum)/10
Rh5.8<-apply(Rh[c(6:9)],1,mean)

climate5.8<-data.frame(year=c(1954:2012),MTemp5.8,MaxTemp5.8,MinTemp5.8,
                       Pre5.8,Rh5.8)
## 2.1.1 MKK test for the trend------------------
## 
MTemp5.9.ts<-ts(MTemp5.9,start=1954,end=2012,frequency = 1)

mk.test(MTemp5.9, alternative = c("two.sided", "greater", "less"), continuity = TRUE)
bbsmk(MTemp5.9,ci = 0.95,nsim = 2000,eta = 1)## this is the Modified MK test accounting the autocorrelation in the tiem sereies.
bbsmk(Pre5.9,ci = 0.95,nsim = 2000,eta = 1)


plot(SMA(MTemp5.9.ts,n=20))

## 2.1.2 plot the climate variability data -----------------
library(plotrix)
tiff("./plot/Figure 1 color. Climate in the sampling site.tiff",width=18,height = 21,units = "cm",compression = "lzw",res = 400)
windowsFonts(TN = windowsFont("Times New Roman"))  
par(mfcol=c(2,1),mgp=c(2.0,0.5,0),family="TN",ps=13)
par(mar=c(3, 3, 0.5, 7) + 0.1)
month<-1:12

mp<-barplot(apply(Pre[-1],2,mean)/10,width=0.9,space=0.1,beside=TRUE,xlim=c(0,12),
            ylim=c(0,100),xaxs="i",xlab = "Month",ylab = "Precipitation (mm)")
par(new=TRUE)
plot(mp[,1],apply(MTemp[-1],2,mean)/10,xlim=c(0.5,11.5),type="b",col="red",axes=FALSE,ylim=c(-15,20),ann=FALSE)

# We use the "pretty" function go generate nice axes
axis(at = pretty(c(-15:20)), side = 4,col=2)
mtext(4,text="Temperature (¡æ)",col="red",line=2)
par(new=TRUE)
plot(mp[,1],apply(Rh[-1],2,mean),xlim=c(0.5,11.5),type="b",pch=2,col="blue",axes=FALSE,ylim=c(20,90),ann=FALSE)
axis(at = pretty(c(20:90)), side = 4,col="blue",line =3.5)
mtext(4,text="Relative humidity (%)",col="blue",line=5)
box()
text(0.5,88,labels = "a",cex = 2.5,font= 2)

par(mar=c(3, 3, 0, 7) + 0.1,xaxs="i")

#mp1<-barplot(Pre5.9,width=0.05,space=0.001,beside=TRUE,col="gray90",#names.arg = pretty(MTemp[1]),
#            ylim=c(0,600),xaxs="i",xlab = "Year",ylab = "May-Sep Precipitation (mm)")
#axis.break(2,200,style="gap")
plot(t(MTemp[1]),Pre5.9,xlim = c(1954,2010),"l",
        ylim=c(200,600),xlab = "Year",ylab = "May-Sep precipitation (mm)")
box()
par(new=TRUE)
plot(t(MTemp[1]),MTemp5.9,xlim=c(1954,2010),"l",ylim=c(9,17),axes=F,ann=FALSE,col=2,
     xlab = "Year",ylab = "Temperature (¡æ)")
axis(at = pretty(c(9:17)), side = 4,col=2)
abline(lm(MTemp5.9[c(1:57)]~c(1954:2010)),lty=2,col=2)
# text(1980,9.5,labels = expression(paste(italic("slope"), "= 0.034, ", italic("R")^"2", "= 0.44, ", italic(" p"), " < 0.001")),col=2)
text(1980,9.5,labels = expression(paste(italic("Sen's slope"), "= 0.035, ",  italic(" p"), " < 0.001")),col=2)
mtext(4,text="Temperature (¡æ)",col="red",line=2)

par(new=TRUE)
plot(t(MTemp[1]),Rh5.9,xlim=c(1954,2010),type="l",col="blue",axes=FALSE,ylim=c(50,90),ann=FALSE)
#abline(lm(Rh5.9~c(1954:2012)),lty=2,col=4)
# We use the "pretty" function go generate nice axes
axis(at = pretty(c(20:90)), side = 4,col="blue",line =3.5)
mtext(4,text="Relative humidity (%)",col="blue",line=5)
#box()
text(1956,88,labels = "b",cex =2.5,font = 2)

dev.off()

## 2.1.3 load the carbon and oxygen isotopoe data one-way ANOVA and plot----
## 
## a. load the carbon isotope data------
## 
## here, ccor: the carbon isotope data were correected "Suess effect",
## pin13c: "pin" corrected tree-ring stable isotope.
##  "Pin" correction was completed in MATLAB using the code from  McCarroll, 2010; McCarroll et al., 2009.

ccor<-read_excel("C13.xlsx",sheet="cor")
pin13c<-read_excel("C13.xlsx",sheet="pin",skip = 1)

describe(ccor)
describe(pin13c)

## b. load the stable oxygen isotope data-----
oxy<-read_xlsx("oxy.xlsx",sheet="oxy",skip = 1)
oxy1<-oxy

oxyall<- melt(as.data.frame(oxy1),id.vars ="Year",
              obsoxy=c("Young", "Middle", "Old"))

ccorall<- melt(as.data.frame(ccor),id.vars = "Year",
               obsccor=c("Young", "Middle", "Old"))


pin1 <-as.data.frame(pin13c)

pinall<- melt(pin1,id.vars ="Year",
              obspin=c("Young", "Middle", "Old"))

# Compute the analysis of variance
res.aov <- aov(value ~ variable, data = oxyall)
# Summary of the analysis
summary(res.aov)

#Tukey multiple pairwise-comparisons

TukeyHSD(res.aov)# takes the fitted ANOVA as an argument.
glht(res.aov)
summary(glht(res.aov, linfct = mcp(variable = "Tukey")))

pin.aov <- aov(value ~ variable, data = pinall)
summary(pin.aov)

#Tukey multiple pairwise-comparisons
TukeyHSD(pin.aov)# takes the fitted ANOVA as an argument.
glht(pin.aov)
summary(glht(pin.aov, linfct = mcp(variable = "Tukey")))

ccor.aov <- aov(value ~ variable, data = ccorall)
summary(ccor.aov)

#Tukey multiple pairwise-comparisons
TukeyHSD(ccor.aov)# takes the fitted ANOVA as an argument.
glht(ccor.aov)
summary(glht(ccor.aov, linfct = mcp(variable = "Tukey")))

## 2.1.3.1 ANOVA plot-----
 cormeanplot<-
   ggplot(ccorall, aes(x = variable, y =value,
                       fill = variable)) + 
  geom_violin(trim = TRUE)+
  geom_boxplot(width=0.1, fill="white")+
   stat_summary(fun.y=mean, colour="orange", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  labs(title=" ",x="Age group", 
       y = expression(paste(delta ^13,"C (¡ë)")))+
   theme(legend.position = "none",
         plot.margin = unit(c(-0.4,0.10,-0.5,-0.1),"cm"))+
   mythemeplot()+
   stat_compare_means(method = "anova",
                      label.x="Middle", label.y = -23.5)+ 
   annotate("text", x =c("Young", "Middle", "Old"), 
            #y = c(-21.5,-21.1,-21.3)+0.5, 
            y = c(-20.6),
            label =  c("A","B","B"),
  parse = TRUE)
 
 pinmeanplot<-
   ggplot(pinall, aes(x = variable, y =value,fill = variable)) + 
  geom_violin(trim = TRUE)+
  geom_boxplot(width=0.1, fill="white")+
   stat_summary(fun.y=mean, colour="orange", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
  labs(title=" ",x="Age group", 
       y = expression(paste(delta ^13,"C (¡ë)")))+
   theme(legend.position = "none",
         plot.margin = unit(c(-0.2,0.1,-0.5,-0.1),"cm"))+
   mythemeplot()+
   stat_compare_means(method = "anova",
                      label.x="Middle", label.y = -23)+  # Add global annova p-value
   annotate("text", x =c("Young", "Middle", "Old"), 
            #y = c(-21.4,-21,-20.8)+0.5, 
            y=c(-20.3),
            label =  c("A","B","C"),hjust=c(1.2,0.5,0.5),
  parse = TRUE)
 
 oxymeanplot<-
   ggplot(oxyall, aes(x = variable, y =value,fill = variable)) + 
  geom_violin(trim = TRUE)+
  geom_boxplot(width=0.1, fill="white")+
   stat_summary(fun.y=mean, colour="orange", geom="point", 
               shape=18, size=3,show.legend = FALSE) + 
   
  labs(title=" ",x="Age group", 
       y = expression(paste(delta ^18,"O (¡ë)")))+
   theme(legend.position = "none",
         plot.margin = unit(c(-0.1,0.10,0.0,-0.1),"cm"))+
   mythemeplot()+
   stat_compare_means(method = "anova",
                       label.y = 30)+  # Add global annova p-value
   annotate("text", x =c("Young", "Middle", "Old"), 
            y = c(34.8), 
            label =  c("A","A","B"),
  parse = TRUE)

   
tiff("./plot/Figure 3.1 ANOVA.tiff", width = 21,height = 10,
     units = "cm",compression = "lzw",res = 300,
     bg = "white",family = "serif")
ggarrange( cormeanplot,pinmeanplot,oxymeanplot,
          nrow=1,ncol=3,
          label.x=0.20,label.y=0.93,
          labels=c("a","b","c"),
          font.label = list(size=16,family="serif"))
dev.off()

## 2.1.4 MKK trend test for the proxies setries-----

 oxy.y <- as.vector(unlist(oxy[2]))
 oxy.m <- as.vector(unlist(oxy[3]))
 oxy.o <- as.vector(unlist(oxy[4]))

mk.test(oxy.y, alternative = c("two.sided", "greater", "less"), continuity = TRUE)
bbsmk(oxy.y,ci = 0.95,nsim = 2000,eta = 1)

mk.test(oxy.m, alternative = c("two.sided", "greater", "less"), continuity = TRUE)
bbsmk(oxy.m,ci = 0.95,nsim = 2000,eta = 1)

mk.test(oxy.o, alternative = c("two.sided", "greater", "less"), continuity = TRUE)
bbsmk(oxy.o,ci = 0.95,nsim = 2000,eta = 1)
 
 ccor.y <- as.vector(unlist(ccor[2]))
 ccor.m <- as.vector(unlist(ccor[3]))
 ccor.o <- as.vector(unlist(ccor[4]))
 
 mk.test(ccor.y, alternative = c("two.sided", "greater", "less"), continuity = TRUE)
bbsmk(ccor.y,ci = 0.95,nsim = 2000,eta = 1)

mk.test(ccor.m, alternative = c("two.sided", "greater", "less"), continuity = TRUE)
bbsmk(ccor.m,ci = 0.95,nsim = 2000,eta = 1)

mk.test(ccor.o, alternative = c("two.sided", "greater", "less"), continuity = TRUE)
bbsmk(ccor.o,ci = 0.95,nsim = 2000,eta = 1)
 
 pin.y <- as.vector(unlist(pin[2]))
 pin.m <- as.vector(unlist(pin[3]))
 pin.o <- as.vector(unlist(pin[4]))
 
 mk.test(pin.y, alternative = c("two.sided", "greater", "less"), continuity = TRUE)
bbsmk(pin.y,ci = 0.95,nsim = 2000,eta = 1)

mk.test(pin.m, alternative = c("two.sided", "greater", "less"), continuity = TRUE)
bbsmk(pin.m,ci = 0.95,nsim = 2000,eta = 1)

mk.test(pin.o, alternative = c("two.sided", "greater", "less"), continuity = TRUE)
bbsmk(pin.o,ci = 0.95,nsim = 2000,eta = 1)

##2.1.5 plot the time series of the proxies-----
## 
ccor.p<- ggplot(data=ccorall,aes(x=year,y=value,color=variable))+
                 geom_line()+
       labs(title=" ",x=" ", 
       y = expression(paste(delta ^13,"C (¡ë)")))+
       geom_smooth(method = "lm", se = FALSE,linetype=2)+
       scale_x_continuous(expand=c(0.01,0.02))+
       annotate("text", x = 1940, y=c(-20.9,-20.6,-20.3), 
             label = c(expression(paste(italic("Sen's slope"), "= -0.0072, ",  italic(" p"), " < 0.001")),
                       expression(paste(italic("Sen's slope"), "= -0.0076, ",  italic(" p"), " < 0.001")),
                       expression(paste(italic("Sen's slope"), "= -0.0107, ",  italic(" p"), " < 0.001"))),
             col=c("#F8766D","#00BA38","#619CFF"))+## old, middle and young
  mythemeplot()+
  theme(legend.title = element_blank(),
        legend.position=c(0.6,0.1),
        plot.margin = unit(c(-0.4,0.22,-0.5,0.1),"cm"))+
  guides(col=guide_legend(ncol = 3,order=1))

pin.p<- ggplot(data=pinall,aes(x=year,y=value,color=variable))+
                 geom_line()+
       labs(title=" ",x=" ", 
       y = expression(paste(delta ^13,"C (¡ë)")))+
       geom_smooth(method = "lm", se = FALSE,linetype=2)+
       scale_x_continuous(expand=c(0.01,0.02))+
       annotate("text", x = 1940, y=c(-20.9,-20.6,-20.3), 
             label = c(expression(paste(italic("Sen's slope"), "= 0.0037, ",  italic(" p"), " = 0.019")),
                       expression(paste(italic("Sen's slope"), "= 0.0039, ",  italic(" p"), " = 0.004")),
                       expression(paste(italic("Sen's slope"), "= 0.0043, ",  italic(" p"), " < 0.001"))),
             col=c("#F8766D","#00BA38","#619CFF"))+## old, middle and young
  mythemeplot()+
  theme(legend.title = element_blank(), 
        legend.position=c(0.6,0.1),
        plot.margin = unit(c(-0.2,0.22,-0.5,0.1),"cm"))+
  guides(col=guide_legend(ncol = 3,order=1))
  
  oxy.p<- ggplot(data=oxyall,aes(x=Year,y=value,color=variable))+
                 geom_line()+
       labs(title=" ",x="Year", 
       y = expression(paste(delta ^18,"O (¡ë)")))+
       scale_x_continuous(expand=c(0.01,0.02))+
  mythemeplot()+
  theme(legend.title = element_blank(), 
        legend.position=c(0.6,0.1),
        plot.margin = unit(c(-0.1,0.22,0,0.1),"cm")) +
    guides(col=guide_legend(ncol = 3))

  ##2.1.6 output the proxies series and ANOVA results------

  tiff("./plot/Figure 3.2 time series and ANOVA.tiff", width = 21,height = 18,
     units = "cm",compression = "lzw",res = 300,
     bg = "white",family = "serif")
ggarrange( ccor.p,cormeanplot,
           pin.p,pinmeanplot,
           oxy.p,oxymeanplot,
          nrow=3,ncol=2,widths = c(2,1),
          heights = c(1,1,1.1),
          align = "v",
          label.x=c(0.1,0.20,0.1,0.2,0.1,0.2),
          label.y=0.93,
          labels=c("a","d","b","e","c","f"),
          common.legend = TRUE,legend="bottom",
          font.label = list(size=20,family="serif"))
dev.off()
# 2.2 climatic response ---------------------------------------------------

# 2.2.1 Chronology transformat pin13c-----------------------------

pin <- as.data.frame(pin13c[2:4])
rownames(pin)<-c(pin13c$Year)# rename the row names.

group<-c(rep("Young",15),rep("Medium",15),rep("Old",15))


# 2.2.2 Calculate pin13C cc and climate response--------------------------------
pin.1954<-subset(pin,rownames(pin)>1953)

pinA1<-subset.data.frame(pin.1954,select=1)
pinA2<-subset.data.frame(pin.1954,select=2)
pinA3<-subset.data.frame(pin.1954,select=3)

pinA1.mtresp <- dcc(pinA1,data.frame(MTemp), method = "corr",selection = .range(-8:9)+.mean(5:8))
pinA2.mtresp <- dcc(pinA2,data.frame(MTemp), method = "corr",selection = .range(-8:9)+.mean(5:8))
pinA3.mtresp <- dcc(pinA3,data.frame(MTemp), method = "corr",selection =.range(-8:9)+.mean(5:8))

pinA1.maxtresp <- dcc(pinA1,data.frame(MaxTemp), method = "corr",selection =.range(-8:9)+.mean(5:8))
pinA2.maxtresp <- dcc(pinA2,data.frame(MaxTemp), method = "corr",selection =.range(-8:9)+.mean(5:8))
pinA3.maxtresp <- dcc(pinA3,data.frame(MaxTemp), method = "corr",selection = .range(-8:9)+.mean(5:8))

pinA1.mintresp <- dcc(pinA1,data.frame(MinTemp), method = "corr",selection = .range(-8:9)+.mean(5:8))
pinA2.mintresp <- dcc(pinA2,data.frame(MinTemp), method = "corr",selection = .range(-8:9)+.mean(5:8))
pinA3.mintresp <- dcc(pinA3,data.frame(MinTemp), method = "corr",selection = .range(-8:9)+.mean(5:8))

pinA1.preresp <- dcc(pinA1,data.frame(Pre), method = "corr",selection = .range(-8:9)+.sum(5:8))
pinA2.preresp <- dcc(pinA2,data.frame(Pre), method = "corr",selection = .range(-8:9)+.sum(5:8))
pinA3.preresp <- dcc(pinA3,data.frame(Pre), method = "corr",selection = .range(-8:9)+.sum(5:8))

pinA1.rhresp <- dcc(pinA1,data.frame(Rh), method = "corr",selection = .range(-8:9)+.mean(5:8))
pinA2.rhresp <- dcc(pinA2,data.frame(Rh), method = "corr",selection = .range(-8:9)+.mean(5:8))
pinA3.rhresp <- dcc(pinA3,data.frame(Rh), method = "corr",selection = .range(-8:9)+.mean(5:8))

pin.MTemp<-rbind(pinA1.mtresp$coef,pinA2.mtresp$coef,pinA3.mtresp$coef)
pin.MaxTemp<-rbind(pinA1.maxtresp$coef,pinA2.maxtresp$coef,pinA3.maxtresp$coef)
pin.MinTemp<-rbind(pinA1.mintresp$coef,pinA2.mintresp$coef,pinA3.mintresp$coef)
pin.PreTemp<-rbind(pinA1.preresp$coef,pinA2.preresp$coef,pinA3.preresp$coef)
pin.RhTemp<-rbind(pinA1.rhresp$coef,pinA2.rhresp$coef,pinA3.rhresp$coef)

pin.MTemp$group<-group;
pin.MaxTemp$group<-group
pin.MinTemp$group<-group
pin.PreTemp$group<-group
pin.RhTemp$group<-group

pin.res<-cbind(pin.MTemp,pin.MaxTemp,pin.MinTemp,pin.PreTemp,pin.RhTemp)

##2.2.3 read the oxygen data and calculate the climate response------
oxy2<-as.data.frame(oxy[2:4])
row.names(oxy2)<-c(1910:2010)
describe(oxy2)

oxyA1<-subset.data.frame(oxy2,select=1)
oxyA2<-subset.data.frame(oxy2,select=2)
oxyA3<-subset.data.frame(oxy2,select=3)

####bbsmk(oxy,ci = 0.95,nsim = 2000,eta = 1)---

oxyA1.mtresp <- dcc(oxyA1,data.frame(MTemp), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxyA2.mtresp <- dcc(oxyA2,data.frame(MTemp), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxyA3.mtresp <- dcc(oxyA3,data.frame(MTemp), method = "corr",selection = .range(-8:9)+.mean(5:8))

oxyA1.maxtresp <- dcc(oxyA1,data.frame(MaxTemp), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxyA2.maxtresp <- dcc(oxyA2,data.frame(MaxTemp), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxyA3.maxtresp <- dcc(oxyA3,data.frame(MaxTemp), method = "corr",selection =.range(-8:9)+.mean(5:8))

oxyA1.mintresp <- dcc(oxyA1,data.frame(MinTemp), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxyA2.mintresp <- dcc(oxyA2,data.frame(MinTemp), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxyA3.mintresp <- dcc(oxyA3,data.frame(MinTemp), method = "corr",selection = .range(-8:9)+.mean(5:8))

oxyA1.preresp <- dcc(oxyA1,data.frame(Pre), method = "corr",selection =.range(-8:9)+.sum(5:8))
oxyA2.preresp <- dcc(oxyA2,data.frame(Pre), method = "corr",selection = .range(-8:9)+.sum(5:8))
oxyA3.preresp <- dcc(oxyA3,data.frame(Pre), method = "corr",selection =.range(-8:9)+.sum(5:8))

oxyA1.rhresp <- dcc(oxyA1,data.frame(Rh), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxyA2.rhresp <- dcc(oxyA2,data.frame(Rh), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxyA3.rhresp <- dcc(oxyA3,data.frame(Rh), method = "corr",selection = .range(-8:9)+.mean(5:8))

oxyA1.mtresp.moving <- dcc(oxyA1,data.frame(MTemp), method = "corr",selection = .range(2:9)+.mean(5:8), moving=TRUE)
oxyA2.mtresp.moving <- dcc(oxyA2,data.frame(MTemp), method = "corr",selection = .range(2:9)+.mean(5:8), moving=TRUE)
oxyA3.mtresp.moving <- dcc(oxyA3,data.frame(MTemp), method = "corr",selection = .range(2:9)+.mean(5:8), moving=TRUE)


oxyA1.maxtresp.moving <- dcc(oxyA1,data.frame(MaxTemp), method = "corr",selection = .range(2:9)+.mean(5:8), moving=TRUE)
oxyA2.maxtresp.moving <- dcc(oxyA2,data.frame(MaxTemp), method = "corr",selection = .range(2:9)+.mean(5:8), moving=TRUE)
oxyA3.maxtresp.moving <- dcc(oxyA3,data.frame(MaxTemp), method = "corr",selection =.range(2:9)+.mean(5:8), moving=TRUE)

oxyA1.mintresp.moving <- dcc(oxyA1,data.frame(MinTemp), method = "corr",selection = .range(2:9)+.mean(5:8), moving=TRUE)
oxyA2.mintresp.moving <- dcc(oxyA2,data.frame(MinTemp), method = "corr",selection = .range(2:9)+.mean(5:8), moving=TRUE)
oxyA3.mintresp.moving <- dcc(oxyA3,data.frame(MinTemp), method = "corr",selection = .range(2:9)+.mean(5:8), moving=TRUE)

oxyA1.preresp.moving <- dcc(oxyA1,data.frame(Pre), method = "corr",selection =.range(2:9)+.sum(5:8), moving=TRUE)
oxyA2.preresp.moving <- dcc(oxyA2,data.frame(Pre), method = "corr",selection = .range(2:9)+.sum(5:8), moving=TRUE)
oxyA3.preresp.moving <- dcc(oxyA3,data.frame(Pre), method = "corr",selection =.range(2:9)+.sum(5:8), moving=TRUE)

oxyA1.rhresp.moving <- dcc(oxyA1,data.frame(Rh), method = "corr",selection = .range(2:9)+.mean(5:8), moving=TRUE)
oxyA2.rhresp.moving <- dcc(oxyA2,data.frame(Rh), method = "corr",selection = .range(2:9)+.mean(5:8), moving=TRUE)
oxyA3.rhresp.moving <- dcc(oxyA3,data.frame(Rh), method = "corr",selection = .range(2:9)+.mean(5:8), moving=TRUE)
plot(oxyA3.rhresp.moving)

oxyA3.rhresp.season <- seascorr(oxyA1,climate=list(data.frame(Rh),data.frame(MTemp)),var_names = c("RH","Tem"))
plot(oxyA3.rhresp.season)

oxy.MTemp<-rbind(oxyA1.mtresp$coef,oxyA2.mtresp$coef,oxyA3.mtresp$coef)
oxy.MaxTemp<-rbind(oxyA1.maxtresp$coef,oxyA2.maxtresp$coef,oxyA3.maxtresp$coef)
oxy.MinTemp<-rbind(oxyA1.mintresp$coef,oxyA2.mintresp$coef,oxyA3.mintresp$coef)
oxy.PreTemp<-rbind(oxyA1.preresp$coef,oxyA2.preresp$coef,oxyA3.preresp$coef)
oxy.RhTemp<-rbind(oxyA1.rhresp$coef,oxyA2.rhresp$coef,oxyA3.rhresp$coef)

oxy.MTemp$group<-group
oxy.MaxTemp$group<-group
oxy.MinTemp$group<-group
oxy.PreTemp$group<-group
oxy.RhTemp$group<-group

oxy.res<-cbind(oxy.MTemp,oxy.MaxTemp,oxy.MinTemp,oxy.PreTemp,oxy.RhTemp)

##2.3 calculate the bootstrap corelation between chronlogy and climate and test the sensitivity between groups-------

## here, we need to call my function "climres_com"
## in the climres_com function, we used the bootstrap method to compare the difference in the correlations by Kruskal-Wallis test
## 

RH.oxy.climres.com <- climres_com(sty=1954, edy=2010, oxy=oxy2, Rh=Rh, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
RH.oxy.climres.season <- climres_com(sty=1954, edy=2010, oxy=oxy2, Rh=climate5.8, Brep = 1000,stm = 6,edm = 6,stmp = 6,edmp = 6)
RH.oxy.clim.com<-rbind(RH.oxy.climres.com[[2]],RH.oxy.climres.com[[4]],RH.oxy.climres.season[[4]])
RH.oxy.clim.com.sig <- data.frame(matrix(RH.oxy.clim.com$difference, ncol =3,byrow = TRUE))

MTemp.oxy.climres.com <- climres_com(sty=1954, edy=2010, oxy=oxy2, Rh=MTemp, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
MTemp.oxy.climres.season<-climres_com(sty=1954, edy=2010, oxy=oxy2, Rh=climate5.8, Brep = 1000,stm = 2,edm = 2,stmp = 2,edmp = 2)
MTemp.oxy.clim.com<-rbind(MTemp.oxy.climres.com[[2]],MTemp.oxy.climres.com[[4]],MTemp.oxy.climres.season[[4]])
MTemp.oxy.clim.com.sig <- data.frame(matrix(MTemp.oxy.clim.com$difference, ncol =3,byrow = TRUE))

MaxTemp.oxy.climres.com<-climres_com(sty=1954, edy=2010, oxy=oxy2, Rh=MaxTemp, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
MaxTemp.oxy.climres.season<-climres_com(sty=1954, edy=2010, oxy=oxy2, Rh=climate5.8, Brep = 1000,stm = 3,edm = 3,stmp = 3,edmp = 3)
MaxTemp.oxy.clim.com<-rbind(MaxTemp.oxy.climres.com[[2]],MaxTemp.oxy.climres.com[[4]],MaxTemp.oxy.climres.season[[4]])
MaxTemp.oxy.clim.com.sig<-data.frame(matrix(MaxTemp.oxy.clim.com$difference, ncol =3,byrow = TRUE))

MinTemp.oxy.climres.com<-climres_com(sty=1954, edy=2010, oxy=oxy2, Rh=MinTemp, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
MinTemp.oxy.climres.season<-climres_com(sty=1954, edy=2010, oxy=oxy2, Rh=climate5.8, Brep = 1000,stm = 4,edm = 4,stmp = 4,edmp = 4)
MinTemp.oxy.clim.com<-rbind(MinTemp.oxy.climres.com[[2]],MinTemp.oxy.climres.com[[4]],MinTemp.oxy.climres.season[[4]])
MinTemp.oxy.clim.com.sig<-data.frame(matrix(MinTemp.oxy.clim.com$difference, ncol =3,byrow = TRUE))

Pre.oxy.climres.com<-climres_com(sty=1954, edy=2010, oxy=oxy2, Rh=Pre, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
Pre.oxy.climres.season<-climres_com(sty=1954, edy=2010, oxy=oxy2, Rh=climate5.8, Brep = 1000,stm = 5,edm = 5,stmp =5,edmp = 5)
Pre.oxy.clim.com<-rbind(Pre.oxy.climres.com[[2]],Pre.oxy.climres.com[[4]],Pre.oxy.climres.season[[4]])
Pre.oxy.clim.com.sig<-data.frame(matrix(Pre.oxy.clim.com$difference, ncol =3,byrow = TRUE))

## for pin chronologies

RH.pin.climres.com<-climres_com(sty=1954, edy=2010, oxy=pin, Rh=Rh, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
RH.pin.climres.season<-climres_com(sty=1954, edy=2010, oxy=pin, Rh=climate5.8, Brep = 1000,stm = 6,edm = 6,stmp = 6,edmp = 6)
RH.pin.clim.com<-rbind(RH.pin.climres.com[[2]],RH.pin.climres.com[[4]],RH.pin.climres.season[[4]])
RH.pin.clim.com.sig<-data.frame(matrix(RH.pin.clim.com$difference, ncol =3,byrow = TRUE))

MTemp.pin.climres.com<-climres_com(sty=1954, edy=2010, oxy=pin, Rh=MTemp, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
MTemp.pin.climres.season<-climres_com(sty=1954, edy=2010, oxy=pin, Rh=climate5.8, Brep = 1000,stm = 2,edm = 2,stmp = 2,edmp = 2)
MTemp.pin.clim.com<-rbind(MTemp.pin.climres.com[[2]],MTemp.pin.climres.com[[4]],MTemp.pin.climres.season[[4]])
MTemp.pin.clim.com.sig<-data.frame(matrix(MTemp.pin.clim.com$difference, ncol =3,byrow = TRUE))

MaxTemp.pin.climres.com<-climres_com(sty=1954, edy=2010, oxy=pin, Rh=MaxTemp, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
MaxTemp.pin.climres.season<-climres_com(sty=1954, edy=2010, oxy=pin, Rh=climate5.8, Brep = 1000,stm = 3,edm = 3,stmp = 3,edmp = 3)
MaxTemp.pin.clim.com<-rbind(MaxTemp.pin.climres.com[[2]],MaxTemp.pin.climres.com[[4]],MaxTemp.pin.climres.season[[4]])
MaxTemp.pin.clim.com.sig<-data.frame(matrix(MaxTemp.pin.clim.com$difference, ncol =3,byrow = TRUE))

MinTemp.pin.climres.com<-climres_com(sty=1954, edy=2010, oxy=pin, Rh=MinTemp, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
MinTemp.pin.climres.season<-climres_com(sty=1954, edy=2010, oxy=pin, Rh=climate5.8, Brep = 1000,stm = 4,edm = 4,stmp = 4,edmp = 4)
MinTemp.pin.clim.com<-rbind(MinTemp.pin.climres.com[[2]],MinTemp.pin.climres.com[[4]],MinTemp.pin.climres.season[[4]])
MinTemp.pin.clim.com.sig<-data.frame(matrix(MinTemp.pin.clim.com$difference, ncol =3,byrow = TRUE))

Pre.pin.climres.com<-climres_com(sty=1954, edy=2010, oxy=pin, Rh=Pre, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
Pre.pin.climres.season<-climres_com(sty=1954, edy=2010, oxy=pin, Rh=climate5.8, Brep = 1000,stm = 5,edm = 5,stmp =5,edmp = 5)
Pre.pin.clim.com<-rbind(Pre.pin.climres.com[[2]],Pre.pin.climres.com[[4]],Pre.pin.climres.season[[4]])
Pre.pin.clim.com.sig<-data.frame(matrix(Pre.pin.clim.com$difference, ncol =3,byrow = TRUE))

## compare group and data
compare_group<-c(rep("Young-Medium",15),rep("Young-Old",15),rep("Medium-Old",15))
y.position<-c(rep(0.05,15),rep(0.1,15),rep(0.15,15))

MTemp.oxy.clim.com.sig.long<-melt(MTemp.oxy.clim.com.sig)
MaxTemp.oxy.clim.com.sig.long<-melt(MaxTemp.oxy.clim.com.sig)
MinTemp.oxy.clim.com.sig.long<-melt(MinTemp.oxy.clim.com.sig)
Pre.oxy.clim.com.sig.long<-melt(Pre.oxy.clim.com.sig)
Rh.oxy.clim.com.sig.long<-melt(RH.oxy.clim.com.sig)

oxy.MTemp$compare_group<-compare_group
oxy.MaxTemp$compare_group<-compare_group
oxy.MinTemp$compare_group<-compare_group
oxy.PreTemp$compare_group<-compare_group
oxy.RhTemp$compare_group<-compare_group

oxy.MTemp$sig<-MTemp.oxy.clim.com.sig.long$value
oxy.MaxTemp$sig<-MaxTemp.oxy.clim.com.sig.long$value
oxy.MinTemp$sig<-MinTemp.oxy.clim.com.sig.long$value
oxy.PreTemp$sig<-Pre.oxy.clim.com.sig.long$value
oxy.RhTemp$sig<-Rh.oxy.clim.com.sig.long$value

oxy.MTemp$y.position<-y.position
oxy.MaxTemp$y.position<-y.position
oxy.MinTemp$y.position<-y.position
oxy.PreTemp$y.position<-y.position
oxy.RhTemp$y.position<-y.position

MTemp.pin.clim.com.sig.long<-melt(MTemp.pin.clim.com.sig)
MaxTemp.pin.clim.com.sig.long<-melt(MaxTemp.pin.clim.com.sig)
MinTemp.pin.clim.com.sig.long<-melt(MinTemp.pin.clim.com.sig)
Pre.pin.clim.com.sig.long<-melt(Pre.pin.clim.com.sig)
Rh.pin.clim.com.sig.long<-melt(RH.pin.clim.com.sig)

pin.MTemp$compare_group<-compare_group
pin.MaxTemp$compare_group<-compare_group
pin.MinTemp$compare_group<-compare_group
pin.PreTemp$compare_group<-compare_group
pin.RhTemp$compare_group<-compare_group

pin.MTemp$sig<-MTemp.pin.clim.com.sig.long$value
pin.MaxTemp$sig<-MaxTemp.pin.clim.com.sig.long$value
pin.MinTemp$sig<-MinTemp.pin.clim.com.sig.long$value
pin.PreTemp$sig<-Pre.pin.clim.com.sig.long$value
pin.RhTemp$sig<-Rh.pin.clim.com.sig.long$value

pin.MTemp$y.position<-y.position
pin.MaxTemp$y.position<-y.position
pin.MinTemp$y.position<-y.position
pin.PreTemp$y.position<-y.position
pin.RhTemp$y.position<-y.position

### 2.4 output the climate response -------------

## for the pin series
## 
pinmT.res<-ggplot(data=pin.MTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(pin.MTemp,sig==TRUE & id %in% c(subset(pin.MTemp,significant==TRUE & abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,#(y.position+max(abs(oxy.MTemp$coef))),
                 #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  # geom_bar(position = "dodge",stat = "identity")+
   geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  # scale_color_manual(values=c("Young"="red",
  #                             "Medium"="darkgreen",
  #                             "Old"="royalblue",
  #                             "Young-Medium"="red",
  #                             "Medium-Old"="darkgreen",
  #                             "Young-Old"="royalblue" ))+
  scale_color_discrete(name="Age group")+
  #scale_fill_discrete(breaks=c("Young","Medium","Old"))+ ## change the order of the legend
  geom_errorbar(data=pin.MTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x=" ",y="Correlation coefficient")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position = c(0.6,0.15),legend.margin = margin(-0.2,0,0,0, unit="cm"))+
  guides(col=guide_legend(ncol = 3,order=1),
         shape=guide_legend(ncol = 3,order = 2))+
    #labs(fill="Age group",shape="Comapre group")+
  # geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',
             linetype=2,size=1,alpha=I(1/2))+
  annotate("text", x=3, y=-0.90, label="a T-mean", colour="black",family="TN",size=9,fontface=2)+
  theme(plot.margin = unit(c(0.2,0,-0.5,0),"lines"))
  #annotate("text", x=3, y=0.75, label=" ", colour="black",family="TN",size=9)


pinmaxT.res<-ggplot(data=pin.MaxTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(pin.MaxTemp,sig==TRUE & id %in% c(subset(pin.MaxTemp,significant==TRUE & abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,#(y.position+max(abs(pin.MaxTemp$coef))),
                 #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  # geom_bar(position = "dodge",stat = "identity")+
   geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  # scale_color_manual(values=c("Young"="red",
  #                             "Medium"="darkgreen",
  #                             "Old"="royalblue",
  #                             "Young-Medium"="red",
  #                             "Medium-Old"="darkgreen",
  #                             "Young-Old"="royalblue" ))+
  scale_fill_discrete(name="Age group")+
  geom_errorbar(data=pin.MaxTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x=" ",y="Correlation coefficient")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
   theme(legend.position = "none")+
  # guides(fill=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
  # geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  annotate("text", x=2.5, y=-0.90, label="b T-max", colour="black",family="TN",size=9,fontface=2)+
  theme(plot.margin = unit(c(-0.5,0,-0.5,0),"lines"))

 pinminT.res<-ggplot(data=pin.MinTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(pin.MinTemp,sig==TRUE & id %in% c(subset(pin.MinTemp,significant==TRUE & abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,#(y.position+max(abs(pin.MinTemp$coef))),
                # col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  # geom_bar(position = "dodge",stat = "identity")+
   geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  # scale_color_manual(values=c("Young"="red",
  #                             "Medium"="darkgreen",
  #                             "Old"="royalblue",
  #                             "Young-Medium"="red",
  #                             "Medium-Old"="darkgreen",
  #                             "Young-Old"="royalblue" ))+
  scale_color_discrete(name="Age group")+
  geom_errorbar(data=pin.MinTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x="Month",y="Correlation coefficient")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position = c(0.3,0.15),legend.margin = margin(-0.2,0,0,0, unit="cm"))+
  guides(col=guide_legend(ncol = 3,order=1),
         shape=guide_legend(ncol = 3,order=2))+
  # geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  annotate("text", x=1, y=0.90, label="a", colour="black",family="TN",size=9,fontface=2)+
  theme(plot.margin = unit(c(0.2,0,0.4,0.2),"lines"))
  #annotate("text", x=3.5, y=0.75, label="h T-min", colour="black",family="TN",size=12)

pinpre.res<-ggplot(data=pin.PreTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(pin.PreTemp,sig==TRUE & id %in% c(subset(pin.PreTemp,significant==TRUE & abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,#-1*(y.position+max(abs(pin.PreTemp$coef))),
                 #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  # geom_bar(position = "dodge",stat = "identity")+
   geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  #scale_fill_discrete(name="Age group")+
  # scale_color_manual(values=c("Young"="red",
  #                             "Medium"="darkgreen",
  #                             "Old"="royalblue",
  #                             "Young-Medium"="red",
  #                             "Medium-Old"="darkgreen",
  #                             "Young-Old"="royalblue" ))+
  scale_fill_discrete(name="Age group")+
  geom_errorbar(data=pin.PreTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
   labs(x=" ",y="Correlation coefficient")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position = "none")+
  guides(fill=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
  # geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  annotate("text", x=2, y=-0.90, label="c Pre", colour="black",family="TN",size=9,fontface=2)+
  theme(plot.margin = unit(c(-0.5,0,-0.5,0),"lines"))
  #annotate("text", x=3, y=0.75, label="i Pre", colour="black",family="TN",size=12)

pinRh.res<-ggplot(data=pin.RhTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(pin.RhTemp,sig==TRUE & id %in% c(subset(pin.RhTemp,significant==TRUE& abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,#-1*(y.position+max(abs(pin.RhTemp$coef))),
                 #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  # geom_bar(position = "dodge",stat = "identity")+
   geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  # scale_color_manual(values=c("Young"="red",
  #                             "Medium"="darkgreen",
  #                             "Old"="royalblue",
  #                             "Young-Medium"="red",
  #                             "Medium-Old"="darkgreen",
  #                             "Young-Old"="royalblue" ))+
  scale_fill_discrete(name="Age group")+
  geom_errorbar(data=pin.RhTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
   labs(x="Month",y="Correlation coefficient")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position = "none")+
  guides(fill=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
  #geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  annotate("text", x=2, y=-0.90, label="d RH", colour="black",family="TN",size=9,fontface=2)+
  theme(plot.margin = unit(c(-0.5,0,0.2,0),"lines"))

### for the oxy -- 

oxymT.res<-ggplot(data=oxy.MTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(oxy.MTemp,sig==TRUE & id %in% c(subset(oxy.MTemp,significant==TRUE& abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,#y.position+max(abs(oxy.MTemp$coef)),
             #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
             shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE)+
  theme(legend.title= element_blank())+
  scale_fill_discrete(guide='none')+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+## change the shape by manual
  #scale_fill_discrete(name="Compare group")+## change the label of the legend
  #geom_bar(position = "dodge",stat = "identity")+
  #scale_fill_discrete(name="Age group")+
  #scale_fill_discrete(breaks=c("Young","Medium","Old"))+ ## change the order of the legend
  geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  
  # scale_color_manual(values=c("Young"="red",
  #                             "Medium"="darkgreen",
  #                             "Old"="royalblue",
  #                             "Young-Medium"="red",
  #                             "Medium-Old"="darkgreen",
  #                             "Young-Old"="royalblue" ))+
  scale_color_discrete(name="Age group")+
  geom_errorbar(data=oxy.MTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x="  ",y=" ")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position = "none")+
  #theme(legend.position = c(0.3,0.15),legend.margin = margin(-0.2,0,0,0, unit="cm"))+
  #guides(col=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
  
  #+labs(fill="Age group",shape="Comapre group")+
  #geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  #annotate("text", x=3, y=0.90, label="a T-mean", colour="black",family="TN",size=9)+
  theme(plot.margin = unit(c(0.2,0,-0.5,-0.2),"lines"))


oxymaxT.res<-ggplot(data=oxy.MaxTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(oxy.MaxTemp,sig==TRUE & id %in% c(subset(oxy.MaxTemp,significant==TRUE& abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,#y.position+max(abs(oxy.MaxTemp$coef)),
                 #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
                 position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  #geom_bar(position = "dodge",stat = "identity")+
  geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  scale_fill_discrete(name="Age group")+
  # scale_color_manual(values=c("Young"="red",
  #                             "Medium"="darkgreen",
  #                             "Old"="royalblue",
  #                             "Young-Medium"="red",
  #                             "Medium-Old"="darkgreen",
  #                             "Young-Old"="royalblue" ))+
  geom_errorbar(data=oxy.MaxTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x="  ",y=" ")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position = "none")+
  guides(fill=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
  #geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  #annotate("text", x=3, y=0.90, label="b T-max", colour="black",family="TN",size=9)+
  theme(plot.margin = unit(c(-0.5,0,-0.5,-0.2),"lines"))

oxyminT.res<-ggplot(data=oxy.MinTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(oxy.MinTemp,sig==TRUE & id %in% c(subset(oxy.MinTemp,significant==TRUE & abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,
                 #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  # geom_bar(position = "dodge",stat = "identity")+
  geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  scale_fill_discrete(name="Age group")+
  # scale_color_manual(values=c("Young"="red",
  #                             "Medium"="darkgreen",
  #                             "Old"="royalblue",
  #                             "Young-Medium"="red",
  #                             "Medium-Old"="darkgreen",
  #                             "Young-Old"="royalblue" ))+
  geom_errorbar(data=oxy.MinTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x="Month",y="Correlation coefficient")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position ="none")+
  guides(fill=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
 # geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  annotate("text", x=1, y=0.90, label="b", colour="black",family="TN",size=9,fontface=2)+
  theme(plot.margin = unit(c(0.2,0,0.4,0.2),"lines"))

oxypre.res<-ggplot(data=oxy.PreTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(oxy.PreTemp,sig==TRUE & id %in% c(subset(oxy.PreTemp,significant==TRUE & abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,
                 #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  #geom_bar(position = "dodge",stat = "identity")+
  #scale_fill_discrete(name="Age group")+
  geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  scale_fill_discrete(name="Age group")+
  # scale_color_manual(values=c("Young"="red",
  #                             "Medium"="darkgreen",
  #                             "Old"="royalblue",
  #                             "Young-Medium"="red",
  #                             "Medium-Old"="darkgreen",
  #                             "Young-Old"="royalblue" ))+
  geom_errorbar(data=oxy.PreTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x=" ",y=" ")+
  mythemeplot1(size.num = 12,legend.size = 0.8,a=0,b=0)+
  theme(legend.position = "none")+
  guides(fill=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
  # geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  
  theme(plot.margin = unit(c(-0.5,0,-0.5,-0.2),"lines"))

oxyRh.res<-ggplot(data=oxy.RhTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef,col=group))+ 
  geom_point(data=subset(oxy.RhTemp,sig==TRUE & id %in% c(subset(oxy.RhTemp,significant==TRUE& abs(coef)>0.259)$id)),
             #aes(x=id,y=-1*(y.position-0.02+max(abs(oxy.RhTemp$coef))),
             aes(x=id,y=0.8,
             #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
             shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),col=1,
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_color_manual(name="Compare group")+
  
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  # geom_bar(position = "dodge",stat = "identity")+
  # geom_point(shape = 15, 
  #              size  = 4, 
  #            position =  position_dodge(0.5),
  #            stat = "identity")+
 
  geom_errorbar(data=oxy.RhTemp,aes(ymin=ci_lower, ymax=ci_upper,
                                    col = significant), 
                width=0.4,
                position=position_dodge(0.9),col="gray50")+
  geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  scale_fill_discrete(name="Age group")+
  # scale_color_manual(values=c("Young"="red",
  #                             "Medium"="darkgreen",
  #                             "Old"="royalblue",
  #                             "Young-Medium"="red",
  #                             "Medium-Old"="darkgreen",
  #                             "Young-Old"="royalblue" ))+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x="Month ",y=" ")+
  mythemeplot1(size.num = 12,legend.size = 0.8,a=0,b=0)+
  theme(legend.position = "none")+
  
  #geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
 # annotate("text", x=2.5, y=0.90, label="e RH", colour="black",family="TN",size=9)+
  theme(plot.margin = unit(c(-0.5,0,0.2,-0.2),"lines")) 

library(ggpubr)
tiff("./plot/Figure 5. climate response for the carbon and oxyge isotpe-merge.tiff",width = 30, height = 40,family = "serif",
     units = "cm",compression = "lzw",bg="white",res = 300)
ggarrange(pinmT.res,oxymT.res,  
          pinmaxT.res,oxymaxT.res,
          #pinminT.res,oxyminT.res,
          pinpre.res,oxypre.res,
          pinRh.res,oxyRh.res,
          nrow = 4,ncol = 2)

dev.off()

tiff("./plot/Figure S3 climate response for the carbon and oxyge isotpe-merge.tiff",width = 15, height = 20,family = "serif",
     units = "cm",compression = "lzw",bg="white",res = 300)
ggarrange(pinminT.res,oxyminT.res,
          nrow = 2,ncol = 1)
dev.off()

### 3. correlations and ACF to detect the potential legacy effects in different trees, as well as the response after AR models-------

## 3.1 claculate the correlation between tree-ring series----

oxy.sector <- as.matrix(cbind(oxy2,pin))
colnames(oxy.sector)<-c("oxy.Y","oxy.M","oxy.O","Pin.Y","Pin.M","Pin.O")
cc.proxy<-rcorr(oxy.sector, type="pearson")
corrplot(corr = cc.proxy$r[c(1:6),c(1:6)],order="AOE",type="upper",tl.pos="d")
corrplot(corr=cc.proxy$r[c(1:6),c(1:6)],add=TRUE,type="lower",method = "number", order = "AOE",diag=FALSE,tl.pos="n", cl.pos="n")

cc.proxy1<-rcorr(diff(oxy.sector), type="pearson")
corrplot(corr = cc.proxy1$r[c(1:6),c(1:6)],order="AOE",type="upper",tl.pos="d")
corrplot(corr=cc.proxy1$r[c(1:6),c(1:6)],add=TRUE,type="lower",method = "number", order = "AOE",diag=FALSE,tl.pos="n", cl.pos="n")


## 3.2 calculate ACF and plot for the acf results------------
##3.2.1 ACF for the climate data------
temacf1<-acf(MTemp5.9,type="correlation",plot=TRUE,
             na.action = na.pass,lag.max = 10,ci=c(0.95,0.99))

climacf<-NULL
for(i in 2:6){
precf1<-acf(climate5.8[i],type="correlation",plot=TRUE,
             na.action = na.pass,lag.max = 20,ci=c(0.95,0.99))
climacf<- rbind(climacf, precf1$acf)
}

rownames(climacf)<-c("T-mean", "T-max","T-min", "Pre", "RH")
colnames(climacf) <- c(seq(0:21))


tiff("./plot/Figure S5. clim-acf-ci0.95.tiff",width=11,height = 8,units = "cm",
     compression = "lzw",bg="white",res = 300,family = "serif")
par(mgp=c(2,0.5,0),mar=c(3, 3, 0.2, 0.5) + 0.1)
##col=brewer.pal(5,"Set3") you can use the colorful
barplot(climacf[,c(1:11)],col=c("black","gray30","gray50","gray70","gray90"),border="white",font.axis=1, 
        beside=T, xlab="Lag (years)", axis.lty = 0.05,
        ylab="Autocorrelation coefficient",ylim=c(-0.5,1.05),font.lab=1,names.arg=c(0:10))
abline(h=0,lty=1)
abline(h=c(0.25,-0.25),lty=2)
legend(20,0.95,legend=rownames(climacf),col=c("black","gray30","gray50","gray70","gray90"),
       pch = 15,bty = "n",ncol = 3)
#text(42,0.95,labels = "c",cex=2.5)
box()
#axis(side = 1, at = seq_along(count) - 0.5, tick = FALSE)
dev.off()

###3.2.2 ACF fot the tree-ring data------
oxyacf1<-acf(oxy2[1],type="correlation",plot=TRUE,
    na.action = na.pass,lag.max = 10,ci=0.99)
oxyacf2<-acf(oxy2[2],type="correlation",plot=TRUE,
            na.action = na.pass,lag.max = 10,ci=0.99)
oxyacf3<-acf(oxy2[3],type="correlation",plot=TRUE,
            na.action = na.pass,lag.max = 10,ci=0.99)

pinacf1<-acf(pin[1],type="correlation",plot=TRUE,
    na.action = na.pass,lag.max = 10,ci=0.99)
pinacf2<-acf(pin[2],type="correlation",plot=TRUE,
             na.action = na.pass,lag.max = 10,ci=0.99)
pinacf3<-acf(pin[3],type="correlation",plot=TRUE,
             na.action = na.pass,lag.max = 10,ci=0.99)

oxyacf<-rbind(oxyacf1$acf,oxyacf2$acf,oxyacf3$acf)
pinacf<-rbind(pinacf1$acf,pinacf2$acf,pinacf3$acf)
colnames(oxyacf)=c(oxyacf1$lag)
rownames(oxyacf)=c("Young","Middle","Old")

colnames(pinacf)=c(pinacf1$lag)
rownames(pinacf)=c("Young","Middle","Old")

# 3.2.3 Grouped barplot and correlation output-------

tiff("./plot/Figure 4 carbon oxygen acf-ci0.95.tiff",width=21,height = 18,units = "cm",
     compression = "lzw",bg="white",res = 300)
windowsFonts(TN = windowsFont("Times New Roman"))  
par(mfrow=c(2,2),mgp=c(2.0,0.5,0),family="TN",ps=13)
par(mar=c(0, 0, 0.2, 0) + 0.1)

corrplot(corr = cc.proxy$r[c(1:6),c(1:6)],type="upper",
         col=brewer.pal(n=10, name="PuOr"),cl.lim = c(0, 1),
         tl.pos="d",tl.col = 1,tl.cex=1.2,
         p.mat = cc.proxy$P, sig.level = 0.05,insig ="pch",pch.cex = 3,pch.col = rgb(255, 0, 0,100, maxColorValue=255))
# corrplot(corr=cc.proxy$r[c(1:6),c(1:6)],add=TRUE,type="lower",method = "number", order = "AOE",number.cex = 1,number.font=2,col=1,
#          diag=FALSE,tl.pos="n", cl.pos="n",p.mat = cc.proxy$P, sig.level = 0.05,insig ="pch",pch.cex = 3,pch.col = rgb(255, 0, 0, 100, maxColorValue=255))
mtext("a", side = 3, line =-3.0, outer = FALSE,at=-0.05,
      cex = 2.5, col = 1,font = 2)

corrplot(corr = cc.proxy1$r[c(1:6),c(1:6)],type="upper",
         col=brewer.pal(n=10, name="PuOr"),cl.lim = c(0, 1),
         tl.pos="d",tl.col = 1,tl.cex=1.2,
         p.mat = cc.proxy1$P, sig.level = 0.05,insig ="pch",pch.cex = 3,
         pch.col = rgb(255, 0, 0,100, maxColorValue=255))
# corrplot(corr=cc.proxy1$r[c(1:6),c(1:6)],add=TRUE,type="lower",method = "number", order = "AOE",number.cex = 1,number.font=2,col=1,
#          diag=FALSE,tl.pos="n", cl.pos="n",p.mat = cc.proxy1$P, sig.level = 0.05,insig ="pch",pch.cex = 3,pch.col = rgb(255, 0, 0, 100, maxColorValue=255))

mtext("b", side = 3, line =-3.2, outer = FALSE,at=-0.05,
     cex = 2.5, col = 1,font = 2,ps=11)

par(mar=c(3, 3, 0, 0.5) + 0.1)

barplot(pinacf, col=c("gray40","gray80",1) , border="white",font.axis=1, 
        beside=T, xlab="Lag (years)", 
        ylab="Autocorrelation coefficient",ylim=c(-0.3,1.05),font.lab=1)
abline(h=0,lty=1)
abline(h=c(0.20,-0.20),lty=2)
legend(20,0.9,legend=rownames(oxyacf),pch = 15,col =c("gray40","gray80",1),bty = "n")
text(42,0.95,labels = "c",cex=3.1,font=2)
box()

barplot(oxyacf, col=c("gray40","gray80",1) , border="white",font.axis=1, 
        beside=T, xlab="Lag (years)", 
        ylab="Autocorrelation coefficient",ylim=c(-0.3,1.05),font.lab=1)
abline(h=0,lty=1)
abline(h=c(0.20,-0.20),lty=2)
legend(20,0.9,legend=rownames(oxyacf),pch = 15,col =c("gray40","gray80",1),bty = "n")
text(42,0.95,labels = "d",cex=2.9,font=2)
box()
dev.off()

##3.3 remove the autocorrealtion using the AR(10) model------
##3.3.1 for the stable isotopes series----
## for the carbon pin series
pin.AR<-NULL
for(i in seq(1:3)){
pinA1.AR<-arima(pin[i],order = c(10,0,0),method="CSS")## using the order 3 of AR in the arima moder
pinA1.AR.dat<-pinA1.AR$res
pin.AR <- cbind(pin.AR,pinA1.AR.dat)
}

plot.ts(pin.AR[,1])
lines(pin.AR[,2],col=2)
lines(pin.AR[,3],col=3)

mode(pin.AR)
pin.AR<-data.frame(x = pin.AR) 
colnames(pin.AR)<- c("S1","S2","S4")
rownames(pin.AR)<- c(1910:2010)

oxy.AR<-NULL
for(i in seq(1:3)){
  oxyA1.AR<-arima(oxy2[i],order = c(10,0,0),method="CSS")## using the order 3 of AR in the arima moder
  oxyA1.AR.dat<-oxyA1.AR$res
  oxy.AR <- cbind(oxy.AR,oxyA1.AR.dat)
}
oxy.AR<-data.frame(x = oxy.AR) 
colnames(oxy.AR)<- c("S1","S2","S4")
rownames(oxy.AR) <- c(1910:2010)

## 3.3.2 remove the autocorrelation in the climate data----
Rh.AR<-NULL
for(i in seq(2:13)){
  RhA1.AR<-arima(Rh[i],order = c(10,0,0), method="CSS")## using the order 3 of AR in the arima moder
  RhA1.AR.dat<-RhA1.AR$res
  Rh.AR <- cbind(Rh.AR,RhA1.AR.dat)
}
colnames(Rh.AR)<- c(1:12)
Rh.AR<-data.frame(cbind(year =c(1954:2012),Rh.AR))

MTemp.AR<-NULL
for(i in seq(2:13)){
  MTempA1.AR<-arima(MTemp[i],order = c(10,0,0),method="CSS")
  ## using the order 3 of AR in the arima moder
  MTempA1.AR.dat<-MTempA1.AR$res
  MTemp.AR <- cbind(MTemp.AR,MTempA1.AR.dat)
}
colnames(MTemp.AR)<- c(1:12)
MTemp.AR<-data.frame(cbind(year =c(1954:2012),MTemp.AR))

MaxTemp.AR<-NULL
for(i in seq(2:13)){
  MaxTempA1.AR<-arima(MaxTemp[i],order = c(10,0,0),method="CSS")
  ## using the order 3 of AR in the arima moder
  MaxTempA1.AR.dat<-MaxTempA1.AR$res
  MaxTemp.AR <- cbind(MaxTemp.AR,MaxTempA1.AR.dat)
}
colnames(MaxTemp.AR)<- c(1:12)
MaxTemp.AR<-data.frame(cbind(year =c(1954:2012),MaxTemp.AR))

MinTemp.AR<-NULL
for(i in seq(2:13)){
  MinTempA1.AR<-arima(MinTemp[i],order = c(10,0,0),method="CSS")
  ## using the order 3 of AR in the arima moder
  MinTempA1.AR.dat<-MinTempA1.AR$res
  MinTemp.AR <- cbind(MinTemp.AR,MinTempA1.AR.dat)
}
colnames(MinTemp.AR)<- c(1:12)
MinTemp.AR<-data.frame(cbind(year =c(1954:2012),MinTemp.AR))

Pre.AR<-NULL
for(i in seq(2:13)){
  PreA1.AR<-arima(Pre[i],order = c(10,0,0),method="CSS")
  ## using the order 3 of AR in the arima moder
  PreA1.AR.dat<-PreA1.AR$res
  Pre.AR <- cbind(Pre.AR,PreA1.AR.dat)
}
colnames(Pre.AR)<- c(1:12)
Pre.AR<-data.frame(cbind(year =c(1954:2012),Pre.AR))

head(climate5.8)
climate5.8.AR<-NULL
for(i in seq(2:6)){
  climate5.8A1.AR<-arima(climate5.8[i],order = c(10,0,0),method="CSS")
  ## using the order 3 of AR in the arima moder
  climate5.8A1.AR.dat<-climate5.8A1.AR$res
  climate5.8.AR <- cbind(climate5.8.AR,climate5.8A1.AR.dat)
}
colnames(climate5.8.AR)<- colnames(climate5.8)[c(-1)]
climate5.8.AR<-data.frame(cbind(year =c(1954:2012),climate5.8.AR))

## 3.3.3 climate response after AR----
## 3.3.3.1 Climate response analysis after AR remove----

pin.AR.1954<-subset(pin.AR,rownames(pin.AR)>1953)

pin.ARA1<-subset.data.frame(pin.AR.1954,select=1)
pin.ARA2<-subset.data.frame(pin.AR.1954,select=2)
pin.ARA3<-subset.data.frame(pin.AR.1954,select=3)

pin.ARA1.mtresp <- dcc(pin.ARA1,data.frame(MTemp.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))
pin.ARA2.mtresp <- dcc(pin.ARA2,data.frame(MTemp.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))
pin.ARA3.mtresp <- dcc(pin.ARA3,data.frame(MTemp.AR), method = "corr",selection =.range(-8:9)+.mean(5:8))

pin.ARA1.maxtresp <- dcc(pin.ARA1,data.frame(MaxTemp.AR), method = "corr",selection =.range(-8:9)+.mean(5:8))
pin.ARA2.maxtresp <- dcc(pin.ARA2,data.frame(MaxTemp.AR), method = "corr",selection =.range(-8:9)+.mean(5:8))
pin.ARA3.maxtresp <- dcc(pin.ARA3,data.frame(MaxTemp.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))

pin.ARA1.mintresp <- dcc(pin.ARA1,data.frame(MinTemp.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))
pin.ARA2.mintresp <- dcc(pin.ARA2,data.frame(MinTemp.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))
pin.ARA3.mintresp <- dcc(pin.ARA3,data.frame(MinTemp.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))

pin.ARA1.preresp <- dcc(pin.ARA1,data.frame(Pre.AR), method = "corr",selection = .range(-8:9)+.sum(5:8))
pin.ARA2.preresp <- dcc(pin.ARA2,data.frame(Pre.AR), method = "corr",selection = .range(-8:9)+.sum(5:8))
pin.ARA3.preresp <- dcc(pin.ARA3,data.frame(Pre.AR), method = "corr",selection = .range(-8:9)+.sum(5:8))

pin.ARA1.rhresp <- dcc(pin.ARA1,data.frame(Rh.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))
pin.ARA2.rhresp <- dcc(pin.ARA2,data.frame(Rh.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))
pin.ARA3.rhresp <- dcc(pin.ARA3,data.frame(Rh.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))

pin.AR.MTemp<-rbind(pin.ARA1.mtresp$coef,pin.ARA2.mtresp$coef,pin.ARA3.mtresp$coef)
pin.AR.MaxTemp<-rbind(pin.ARA1.maxtresp$coef,pin.ARA2.maxtresp$coef,pin.ARA3.maxtresp$coef)
pin.AR.MinTemp<-rbind(pin.ARA1.mintresp$coef,pin.ARA2.mintresp$coef,pin.ARA3.mintresp$coef)
pin.AR.PreTemp<-rbind(pin.ARA1.preresp$coef,pin.ARA2.preresp$coef,pin.ARA3.preresp$coef)
pin.AR.RhTemp<-rbind(pin.ARA1.rhresp$coef,pin.ARA2.rhresp$coef,pin.ARA3.rhresp$coef)

pin.AR.MTemp$group<-group;
pin.AR.MaxTemp$group<-group
pin.AR.MinTemp$group<-group
pin.AR.PreTemp$group<-group
pin.AR.RhTemp$group<-group

## for the oxygen data
oxy.ARA1<-subset.data.frame(oxy.AR,select=1)
oxy.ARA2<-subset.data.frame(oxy.AR,select=2)
oxy.ARA3<-subset.data.frame(oxy.AR,select=3)

oxy.ARA1.mtresp <- dcc(oxy.ARA1,data.frame(MTemp.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxy.ARA2.mtresp <- dcc(oxy.ARA2,data.frame(MTemp.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxy.ARA3.mtresp <- dcc(oxy.ARA3,data.frame(MTemp.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))

oxy.ARA1.maxtresp <- dcc(oxy.ARA1,data.frame(MaxTemp.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxy.ARA2.maxtresp <- dcc(oxy.ARA2,data.frame(MaxTemp.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxy.ARA3.maxtresp <- dcc(oxy.ARA3,data.frame(MaxTemp.AR), method = "corr",selection =.range(-8:9)+.mean(5:8))

oxy.ARA1.mintresp <- dcc(oxy.ARA1,data.frame(MinTemp.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxy.ARA2.mintresp <- dcc(oxy.ARA2,data.frame(MinTemp.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxy.ARA3.mintresp <- dcc(oxy.ARA3,data.frame(MinTemp.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))

oxy.ARA1.preresp <- dcc(oxy.ARA1,data.frame(Pre.AR), method = "corr",selection =.range(-8:9)+.sum(5:8))
oxy.ARA2.preresp <- dcc(oxy.ARA2,data.frame(Pre.AR), method = "corr",selection = .range(-8:9)+.sum(5:8))
oxy.ARA3.preresp <- dcc(oxy.ARA3,data.frame(Pre.AR), method = "corr",selection =.range(-8:9)+.sum(5:8))

oxy.ARA1.rhresp <- dcc(oxy.ARA1,data.frame(Rh.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxy.ARA2.rhresp <- dcc(oxy.ARA2,data.frame(Rh.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))
oxy.ARA3.rhresp <- dcc(oxy.ARA3,data.frame(Rh.AR), method = "corr",selection = .range(-8:9)+.mean(5:8))

oxy.AR.MTemp<-rbind(oxy.ARA1.mtresp$coef,oxy.ARA2.mtresp$coef,oxy.ARA3.mtresp$coef)
oxy.AR.MaxTemp<-rbind(oxy.ARA1.maxtresp$coef,oxy.ARA2.maxtresp$coef,oxy.ARA3.maxtresp$coef)
oxy.AR.MinTemp<-rbind(oxy.ARA1.mintresp$coef,oxy.ARA2.mintresp$coef,oxy.ARA3.mintresp$coef)
oxy.AR.PreTemp<-rbind(oxy.ARA1.preresp$coef,oxy.ARA2.preresp$coef,oxy.ARA3.preresp$coef)
oxy.AR.RhTemp<-rbind(oxy.ARA1.rhresp$coef,oxy.ARA2.rhresp$coef,oxy.ARA3.rhresp$coef)

oxy.AR.MTemp$group<-group
oxy.AR.MaxTemp$group<-group
oxy.AR.MinTemp$group<-group
oxy.AR.PreTemp$group<-group
oxy.AR.RhTemp$group<-group

oxy.AR.res<-cbind(oxy.AR.MTemp,oxy.AR.MaxTemp,oxy.AR.MinTemp,oxy.AR.PreTemp,oxy.AR.RhTemp)

##3.3.3.2 calculate the bootstrap corelation between chronlogy and climate and test the sensitivity between groups----

RH.oxy.ARclimres.com<-climres_com(sty=1954, edy=2010, oxy=oxy.AR, Rh=Rh.AR, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
RH.oxy.ARclimres.season<-climres_com(sty=1954, edy=2010, oxy=oxy.AR, Rh=climate5.8.AR, Brep = 1000,stm = 6,edm = 6,stmp = 6,edmp = 6)
RH.oxy.ARclim.com<-rbind(RH.oxy.ARclimres.com[[2]],RH.oxy.ARclimres.com[[4]],RH.oxy.ARclimres.season[[4]])
RH.oxy.ARclim.com.sig<-data.frame(matrix(RH.oxy.ARclim.com$difference, ncol =3,byrow = TRUE))

MTemp.oxy.ARclimres.com<-climres_com(sty=1954, edy=2010, oxy=oxy.AR, Rh=MTemp.AR, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
MTemp.oxy.ARclimres.season<-climres_com(sty=1954, edy=2010, oxy=oxy.AR, Rh=climate5.8.AR, Brep = 1000,stm = 2,edm = 2,stmp = 2,edmp = 2)
MTemp.oxy.ARclim.com<-rbind(MTemp.oxy.ARclimres.com[[2]],MTemp.oxy.ARclimres.com[[4]],MTemp.oxy.ARclimres.season[[4]])
MTemp.oxy.ARclim.com.sig<-data.frame(matrix(MTemp.oxy.ARclim.com$difference, ncol =3,byrow = TRUE))

MaxTemp.oxy.ARclimres.com<-climres_com(sty=1954, edy=2010, oxy=oxy.AR, Rh=MaxTemp.AR, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
MaxTemp.oxy.ARclimres.season<-climres_com(sty=1954, edy=2010, oxy=oxy.AR, Rh=climate5.8.AR, Brep = 1000,stm = 3,edm = 3,stmp = 3,edmp = 3)
MaxTemp.oxy.ARclim.com<-rbind(MaxTemp.oxy.ARclimres.com[[2]],MaxTemp.oxy.ARclimres.com[[4]],MaxTemp.oxy.ARclimres.season[[4]])
MaxTemp.oxy.ARclim.com.sig<-data.frame(matrix(MaxTemp.oxy.ARclim.com$difference, ncol =3,byrow = TRUE))

MinTemp.oxy.ARclimres.com<-climres_com(sty=1954, edy=2010, oxy=oxy.AR, Rh=MinTemp.AR, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
MinTemp.oxy.ARclimres.season<-climres_com(sty=1954, edy=2010, oxy=oxy.AR, Rh=climate5.8.AR, Brep = 1000,stm = 4,edm = 4,stmp = 4,edmp = 4)
MinTemp.oxy.ARclim.com<-rbind(MinTemp.oxy.ARclimres.com[[2]],MinTemp.oxy.ARclimres.com[[4]],MinTemp.oxy.ARclimres.season[[4]])
MinTemp.oxy.ARclim.com.sig<-data.frame(matrix(MinTemp.oxy.ARclim.com$difference, ncol =3,byrow = TRUE))

Pre.oxy.ARclimres.com<-climres_com(sty=1954, edy=2010, oxy=oxy.AR, Rh=Pre.AR, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
Pre.oxy.ARclimres.season<-climres_com(sty=1954, edy=2010, oxy=oxy.AR, Rh=climate5.8.AR, Brep = 1000,stm = 5,edm = 5,stmp =5,edmp = 5)
Pre.oxy.ARclim.com<-rbind(Pre.oxy.ARclimres.com[[2]],Pre.oxy.ARclimres.com[[4]],Pre.oxy.ARclimres.season[[4]])
Pre.oxy.ARclim.com.sig<-data.frame(matrix(Pre.oxy.ARclim.com$difference, ncol =3,byrow = TRUE))

## for pin.AR chronologies

Rh.pin.ARclimres.com<-climres_com(sty=1954, edy=2010, oxy=pin.AR, Rh=Rh.AR, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
Rh.pin.ARclimres.season<-climres_com(sty=1954, edy=2010, oxy=pin.AR, Rh=climate5.8.AR, Brep = 1000,stm = 6,edm = 6,stmp = 6,edmp = 6)
Rh.pin.ARclim.com<-rbind(Rh.pin.ARclimres.com[[2]],Rh.pin.ARclimres.com[[4]],Rh.pin.ARclimres.season[[4]])
Rh.pin.ARclim.com.sig<-data.frame(matrix(Rh.pin.ARclim.com$difference, ncol =3,byrow = TRUE))

MTemp.pin.ARclimres.com<-climres_com(sty=1954, edy=2010, oxy=pin.AR, Rh=MTemp, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
MTemp.pin.ARclimres.season<-climres_com(sty=1954, edy=2010, oxy=pin.AR, Rh=climate5.8.AR, Brep = 1000,stm = 2,edm = 2,stmp = 2,edmp = 2)
MTemp.pin.ARclim.com<-rbind(MTemp.pin.ARclimres.com[[2]],MTemp.pin.ARclimres.com[[4]],MTemp.pin.ARclimres.season[[4]])
MTemp.pin.ARclim.com.sig<-data.frame(matrix(MTemp.pin.ARclim.com$difference, ncol =3,byrow = TRUE))

MaxTemp.pin.ARclimres.com<-climres_com(sty=1954, edy=2010, oxy=pin.AR, Rh=MaxTemp, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
MaxTemp.pin.ARclimres.season<-climres_com(sty=1954, edy=2010, oxy=pin.AR, Rh=climate5.8.AR, Brep = 1000,stm = 3,edm = 3,stmp = 3,edmp = 3)
MaxTemp.pin.ARclim.com<-rbind(MaxTemp.pin.ARclimres.com[[2]],MaxTemp.pin.ARclimres.com[[4]],MaxTemp.pin.ARclimres.season[[4]])
MaxTemp.pin.ARclim.com.sig<-data.frame(matrix(MaxTemp.pin.ARclim.com$difference, ncol =3,byrow = TRUE))

MinTemp.pin.ARclimres.com <- climres_com(sty=1954, edy=2010, oxy=pin.AR, Rh=MinTemp, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
MinTemp.pin.ARclimres.season<-climres_com(sty=1954, edy=2010, oxy=pin.AR, Rh=climate5.8.AR, Brep = 1000,stm = 4,edm = 4,stmp = 4,edmp = 4)
MinTemp.pin.ARclim.com<-rbind(MinTemp.pin.ARclimres.com[[2]],MinTemp.pin.ARclimres.com[[4]],MinTemp.pin.ARclimres.season[[4]])
MinTemp.pin.ARclim.com.sig<-data.frame(matrix(MinTemp.pin.ARclim.com$difference, ncol =3,byrow = TRUE))

Pre.pin.ARclimres.com<-climres_com(sty=1954, edy=2010, oxy=pin.AR, Rh=Pre, Brep = 1000,stm = 2,edm = 10,stmp = 9,edmp = 13)
Pre.pin.ARclimres.season<-climres_com(sty=1954, edy=2010, oxy=pin.AR, Rh=climate5.8.AR, Brep = 1000,stm = 5,edm = 5,stmp =5,edmp = 5)
Pre.pin.ARclim.com<-rbind(Pre.pin.ARclimres.com[[2]],Pre.pin.ARclimres.com[[4]],Pre.pin.ARclimres.season[[4]])
Pre.pin.ARclim.com.sig<-data.frame(matrix(Pre.pin.ARclim.com$difference, ncol =3,byrow = TRUE))

## 3.3.4 output the correlations------
##  3.3.4.1 prepare the data for each subfigure-----
## ## compare group and data

y.position<-c(rep(0.05,15),rep(0.1,15),rep(0.15,15))

MTemp.oxy.ARclim.com.sig.long<-melt(MTemp.oxy.ARclim.com.sig)
MaxTemp.oxy.ARclim.com.sig.long<-melt(MaxTemp.oxy.ARclim.com.sig)
MinTemp.oxy.ARclim.com.sig.long<-melt(MinTemp.oxy.ARclim.com.sig)
Pre.oxy.ARclim.com.sig.long<-melt(Pre.oxy.ARclim.com.sig)
Rh.oxy.ARclim.com.sig.long<-melt(RH.oxy.ARclim.com.sig)

oxy.AR.MTemp$compare_group<-compare_group
oxy.AR.MaxTemp$compare_group<-compare_group
oxy.AR.MinTemp$compare_group<-compare_group
oxy.AR.PreTemp$compare_group<-compare_group
oxy.AR.RhTemp$compare_group<-compare_group

oxy.AR.MTemp$sig<-MTemp.oxy.ARclim.com.sig.long$value
oxy.AR.MaxTemp$sig<-MaxTemp.oxy.ARclim.com.sig.long$value
oxy.AR.MinTemp$sig<-MinTemp.oxy.ARclim.com.sig.long$value
oxy.AR.PreTemp$sig<-Pre.oxy.ARclim.com.sig.long$value
oxy.AR.RhTemp$sig<-Rh.oxy.ARclim.com.sig.long$value

oxy.AR.MTemp$y.position<-y.position
oxy.AR.MaxTemp$y.position<-y.position
oxy.AR.MinTemp$y.position<-y.position
oxy.AR.PreTemp$y.position<-y.position
oxy.AR.RhTemp$y.position<-y.position

MTemp.pin.ARclim.com.sig.long<-melt(MTemp.pin.ARclim.com.sig)
MaxTemp.pin.ARclim.com.sig.long<-melt(MaxTemp.pin.ARclim.com.sig)
MinTemp.pin.ARclim.com.sig.long<-melt(MinTemp.pin.ARclim.com.sig)
Pre.pin.ARclim.com.sig.long<-melt(Pre.pin.ARclim.com.sig)
Rh.pin.ARclim.com.sig.long<-melt(Rh.pin.ARclim.com.sig)

pin.AR.MTemp$compare_group<-compare_group
pin.AR.MaxTemp$compare_group<-compare_group
pin.AR.MinTemp$compare_group<-compare_group
pin.AR.PreTemp$compare_group<-compare_group
pin.AR.RhTemp$compare_group<-compare_group

pin.AR.MTemp$sig<-MTemp.pin.ARclim.com.sig.long$value
pin.AR.MaxTemp$sig<-MaxTemp.pin.ARclim.com.sig.long$value
pin.AR.MinTemp$sig<-MinTemp.pin.ARclim.com.sig.long$value
pin.AR.PreTemp$sig<-Pre.pin.ARclim.com.sig.long$value
pin.AR.RhTemp$sig<-Rh.pin.ARclim.com.sig.long$value

pin.AR.MTemp$y.position<-y.position
pin.AR.MaxTemp$y.position<-y.position
pin.AR.MinTemp$y.position<-y.position
pin.AR.PreTemp$y.position<-y.position
pin.AR.RhTemp$y.position<-y.position

### 3.3.4.2 output the climate response after AR removing-------------
## for pin----
pinmT.ARres<-ggplot(data=pin.AR.MTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(pin.AR.MTemp,sig==TRUE & id %in% c(subset(pin.AR.MTemp,significant==TRUE & abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,#(y.position+max(abs(pin.AR.MTemp$coef))),
                 #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  #geom_bar(position = "dodge",stat = "identity")+
 geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  scale_color_discrete(name="Age group")+
  #scale_fill_discrete(breaks=c("Young","Medium","Old"))+ ## change the order of the legend
  geom_errorbar(data=pin.AR.MTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x=" ",y="Correlation coefficient")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
   theme(legend.position = c(0.6,0.15),legend.margin = margin(-0.2,0,0,0, unit="cm"))+
  #scale_fill_discrete(labels = c("Young", "Medium", "Old"))+
  guides(col=guide_legend(ncol = 3,order=1),
         shape=guide_legend(ncol = 3,order=2))+
  #labs(fill="Age group",shape="Comapre group")+
  #geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  annotate("text", x=2.5, y=-0.85, label="a T-mean", colour="black",family="TN",size=9,fontface=2)+
  theme(plot.margin = unit(c(0.2,0,-0.5,0),"lines"))

pinmaxT.ARres<-ggplot(data=pin.AR.MaxTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(pin.AR.MaxTemp,sig==TRUE & id %in% c(subset(pin.AR.MaxTemp,significant==TRUE& abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,#(y.position+max(abs(pin.AR.MaxTemp$coef))),
                 #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  #geom_bar(position = "dodge",stat = "identity")+
  geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  #scale_color_discrete(name="Age group")+
  geom_errorbar(data=pin.AR.MaxTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x="  ",y="Correlation coefficient ")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position = "none")+
  # guides(fill=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
  # geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  theme(plot.margin = unit(c(-0.5,0,-0.5,0),"lines"))+
  annotate("text", x=2.5, y=-0.85, label="b T-max", colour="black",family="TN",size=9,fontface=2)

pinminT.ARres<-ggplot(data=pin.AR.MinTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(pin.AR.MinTemp,sig==TRUE & id %in% c(subset(pin.AR.MinTemp,significant==TRUE & abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,
                 #(y.position+max(abs(pin.AR.MinTemp$coef))),
                 #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  # geom_bar(position = "dodge",stat = "identity")+
  geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  scale_fill_discrete(name="Age group")+
  geom_errorbar(data=pin.AR.MinTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x="  ",y="Correlation coefficient")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position = "none")+
  guides(fill=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
  # geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  theme(plot.margin = unit(c(-0.5,0,-0.5,0),"lines"))+
  annotate("text", x=2.5, y=-0.85, label="c T-min", colour="black",family="TN",size=9,fontface=2)

pinpre.ARres<-ggplot(data=pin.AR.PreTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(pin.AR.PreTemp,sig==TRUE & id %in% c(subset(pin.AR.PreTemp,significant==TRUE & abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,#-1*(y.position+max(abs(pin.AR.PreTemp$coef))),
                 # col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  #geom_bar(position = "dodge",stat = "identity")+
  geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  scale_color_discrete(name="Age group")+
  geom_errorbar(data=pin.AR.PreTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x="  ",y="Correlation coefficient")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position = "none")+
  guides(fill=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
  #geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  theme(plot.margin = unit(c(-0.5,0,-0.5,0),"lines"))+
 annotate("text", x=2, y=-0.85, label="d Pre", colour="black",family="TN",size=9,fontface=2)

pinRh.ARres<-ggplot(data=pin.AR.RhTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(pin.AR.RhTemp,sig==TRUE & id %in% c(subset(pin.AR.RhTemp,abs(coef)>=0.259)$id)),
             aes(x=id,y=0.8,#-1*(y.position+max(abs(pin.AR.RhTemp$coef))),
                 #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  #geom_bar(position = "dodge",stat = "identity")+
  geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  scale_color_discrete(name="Age group")+
  geom_errorbar(data=pin.AR.RhTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x="Month ",y="Correlation coefficient")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position = "none")+
  guides(fill=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
  #geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  annotate("text", x=2, y=-0.85, label="e RH", colour="black",family="TN",size=9,fontface=2)+
  theme(plot.margin = unit(c(-0.5,0,0.2,0),"lines"))


### for the oxy ---- 
oxymT.ARres<-ggplot(data=oxy.AR.MTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(oxy.AR.MTemp,sig==TRUE & id %in% c(subset(oxy.AR.MTemp,significant==TRUE & abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,#y.position+max(abs(oxy.AR.MTemp$coef)),
                 #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+## change the shape by manual
  #scale_fill_discrete(name="Compare group")+## change the label of the legend
  # geom_bar(position = "dodge",stat = "identity")+
  geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  scale_fill_discrete(name="Age group")+
  #scale_fill_discrete(breaks=c("Young","Medium","Old"))+ ## change the order of the legend
  geom_errorbar(data=oxy.AR.MTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x="  ",y=" ")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position = "none")+
  # theme(legend.position = c(0.3,0.15),legend.margin = margin(-0.2,0,0,0, unit="cm"))+
  # guides(fill=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
  #theme(legend.spacing.y = unit(0, "cm"))+
  
  #+labs(fill="Age group",shape="Comapre group")+
  #geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  # annotate("text", x=3, y=0.90, label="a T-mean", colour="black",family="TN",size=9)+
  theme(plot.margin = unit(c(0.2,0,-0.5,-0.2),"lines"))


oxymaxT.ARres<-ggplot(data=oxy.AR.MaxTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(oxy.AR.MaxTemp,sig==TRUE & id %in% c(subset(oxy.AR.MaxTemp,significant==TRUE & abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,#y.position+max(abs(oxy.AR.MaxTemp$coef)),
                 #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  #geom_bar(position = "dodge",stat = "identity")+
  geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  scale_fill_discrete(name="Age group")+
  geom_errorbar(data=oxy.AR.MaxTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x="  ",y=" ")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position = "none")+
  #geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  #annotate("text", x=3, y=0.90, label="b T-max", colour="black",family="TN",size=9)+
  theme(plot.margin = unit(c(-0.5,0,-0.5,-0.2),"lines"))

oxyminT.ARres<-ggplot(data=oxy.AR.MinTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  # geom_point(data=subset(oxy.AR.MinTemp,sig==TRUE & id %in% c(subset(oxy.AR.MinTemp, abs(coef)>0.259)$id)),
  #            aes(x=id,y=0.8,#y.position+max(abs(oxy.MinTemp$coef)),
  #                # col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
  #                shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
  #            position=position_dodge(0.9),size=4)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  #geom_bar(position = "dodge",stat = "identity")+
  geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  scale_fill_discrete(name="Age group")+
  geom_errorbar(data=oxy.AR.MinTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x="  ",y=" ")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position ="none")+
  guides(fill=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
 # geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  # annotate("text", x=3, y=0.90, label="c T-min", colour="black",family="TN",size=9)+
  theme(plot.margin = unit(c(-0.5,0,-0.5,-0.2),"lines"))

oxypre.ARres<-ggplot(data=oxy.AR.PreTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(oxy.AR.PreTemp,sig==TRUE & id %in% c(subset(oxy.AR.PreTemp,significant==TRUE & abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,#-1*(y.position+max(abs(oxy.AR.PreTemp$coef))),
                # col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  #geom_bar(position = "dodge",stat = "identity")+
  geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  scale_fill_discrete(name="Age group")+
  geom_errorbar(data=oxy.AR.PreTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x=" ",y=" ")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position = "none")+
  guides(fill=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
  #geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  #annotate("text", x=3, y=0.90, label="d Pre", colour="black",family="TN",size=9)+
  theme(plot.margin = unit(c(-0.5,0,-0.5,-0.2),"lines"))

oxyRh.ARres<-ggplot(data=oxy.AR.RhTemp, aes(fill=factor(group,levels = c("Young","Medium", "Old" )),x=id, y=coef))+ 
  geom_point(data=subset(oxy.AR.RhTemp,sig==TRUE & id %in% c(subset(oxy.AR.RhTemp,significant==TRUE & abs(coef)>0.259)$id)),
             aes(x=id,y=0.8,#-1*(y.position-0.02+max(abs(oxy.AR.RhTemp$coef))),
                 #col=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" )),
                 shape=factor(compare_group,levels = c("Young-Medium","Medium-Old", "Young-Old" ))),
             position=position_dodge(0.9),size=3)+
  guides(color=FALSE,fill="none")+
  theme(legend.title= element_blank())+
  scale_shape_manual(name="Compare group",values = c(1,2,8))+
  #scale_shape_discrete(name="Compare group")+## change the label of thew legend
  #geom_bar(position = "dodge",stat = "identity")+
  geom_point(aes(col=factor(group,levels = c("Young","Medium", "Old")),
                   x=id, y=coef),
               shape = 15, 
               size  = 5, 
             position = position_dodge(0.9),
           )+ 
  scale_fill_discrete(name="Age group")+
  geom_errorbar(data=oxy.AR.RhTemp,aes(ymin=ci_lower, ymax=ci_upper,col = significant), width=.4,position=position_dodge(0.9),col="gray50")+
  scale_x_continuous(breaks=c(1:15),labels=c("a","s","o","n","d","J","F","M","A","M","J","J","A","S","M-A"),expand = c(0.005,0.005))+
  scale_y_continuous(limits = c(-0.92,0.92))+
  labs(x="Month ",y="")+
  mythemeplot1(size.num = 12,legend.size = 0.8,0,0)+
  theme(legend.position = "none")+
  guides(fill=guide_legend(ncol = 3),shape=guide_legend(ncol = 3))+
  #geom_hline(yintercept=c(0),colour='gray',linetype=1,size=1)+
  geom_hline(yintercept=c(-0.26,0.26),colour='black',linetype=2,size=1,alpha=I(1/2))+
  # annotate("text", x=2.5, y=0.90, label="e RH", colour="black",family="TN",size=9)+
  theme(plot.margin = unit(c(-0.5,0,0.2,-0.2),"lines")) 

## output
tiff("./plot/Figure S3 climate response for the carbon and oxyge isotpe-AR model.tiff",width = 35, height = 46,
     units = "cm",compression = "lzw",bg="white",
     res = 300,family = "serif")
ggarrange(pinmT.ARres,oxymT.ARres,
          pinmaxT.ARres,oxymaxT.ARres,
          pinminT.ARres,oxyminT.ARres,
          pinpre.ARres,oxypre.ARres,
          pinRh.ARres,oxyRh.ARres,
          nrow = 5,ncol = 2)
dev.off()

### 4. The reconstructions from different combinations of tree-ring proxies----------
###4.1 For the climate response and climate reconstructions-----
## Step 1 z-scale all of the tree-ring stable isotope series 
## Step 2 get the climate response, here we focused on relative humidity--
## step 3 get the RE, CE and R2 to assess the reconstructions
## 
z.oxy<-data.frame(scale(oxy2))
z.pin<-data.frame(scale(pin))

##  All d13c from different trees ages
z.pin.ts<-data.frame(apply(z.pin,1,FUN = mean))

pinall.rhresp <- dcc(z.pin.ts,data.frame(Rh), method = "corr",selection = .range(1:9)+.mean(5:8))
sk.pinall<-skills(pinall.rhresp, target = .mean(5:8), model = "ols", calibration = "-50%",
       timespan = NULL)
plot(sk.pinall)
sk.pinall$full.model$rsquare

sk.pinall1<-skills(pinall.rhresp, target = .mean(5:8), model = "ols", calibration = "49%")
plot(sk.pinall1)
sk.pinall1$full.model$rsquare
sk.pinall1$RE;sk.pinall1$CE
sk.pinall1$r.cal

## d13C for the young trees
pinA1.rhresp <- dcc(z.pin[1],data.frame(Rh), method = "corr",selection = .range(1:9)+.mean(5:8))
sk.pinA1<-skills(pinA1.rhresp, target = .mean(5:8), model = "ols", calibration = "-50%",
                  timespan = NULL)
plot(sk.pinA1)

sk.pinA1.1<-skills(pinA1.rhresp, target = .mean(5:8), model = "ols", calibration = "49%",
                  timespan = NULL)
plot(sk.pinA1.1)
sk.pinA1.1$RE
sk.pinA1.1$CE
sk.pinA1.1$cal.model
sk.pinA1.1$DW

## d13C for the old trees: pinA3
pinA3.rhresp <- dcc(z.pin[3],data.frame(Rh), method = "corr",selection = .range(1:9)+.mean(5:8))
sk.pinA3<-skills(pinA3.rhresp, target = .mean(5:8), model = "ols", calibration = "-50%",
                  timespan = NULL)
plot(sk.pinA3)
sk.pinA3$RE;sk.pinA3$CE

sk.pinA3.1<-skills(pinA3.rhresp, target = .mean(5:8), model = "ols", calibration = "49%",
                  timespan = NULL)
plot(sk.pinA3.1)
sk.pinA3.1$RE;sk.pinA3.1$CE
sk.pinA3.1$full.model

## All d18O from all tree ages
z.oxy.ts<-data.frame(apply(z.oxy,1,FUN = mean))

oxyall.rhresp <- dcc(z.oxy.ts,data.frame(Rh), method = "corr",selection = .range(1:9)+.mean(5:8))
plot(oxyall.rhresp)
sk.oxyall<-skills(oxyall.rhresp, target = .mean(5:8), model = "ols", calibration = "-50%",
                  timespan = NULL)
plot(sk.oxyall)
sk.oxyall$full.model$rsquare
sk.oxyall1<-skills(oxyall.rhresp, target = .mean(5:8), model = "ols", calibration = "49%",
                  timespan = NULL)
plot(sk.oxyall1)
 sk.oxyall1$RE;sk.oxyall1$CE
sk.oxyall1$cal.model

## d18O from old trees:oxyA3
oxyA3.rhresp <- dcc(z.oxy[3],data.frame(Rh), method = "corr",selection = .range(1:9)+.mean(5:8))
sk.oxyA3<-skills(oxyA3.rhresp, target = .mean(5:8), model = "ols", calibration = "-50%",
                 timespan = NULL)
plot(sk.oxyA3)

sk.oxyA3.1<-skills(oxyA3.rhresp, target = .mean(5:8), model = "ols", calibration = "49%",
                 timespan = NULL)
plot(sk.oxyA3.1)
sk.oxyA3.1$RE;sk.oxyA3.1$CE
sk.oxyA3.1$cal.model

## d18O from young trees
oxyA1.rhresp <- dcc(z.oxy[1],data.frame(Rh), method = "corr",selection = .range(1:9)+.mean(5:8))
sk.oxyA1<-skills(oxyA1.rhresp, target = .mean(5:8), model = "ols", calibration = "-50%",
                 timespan = NULL)
plot(sk.oxyA1)
sk.oxyA1$RE;sk.oxyA1$CE
sk.oxyA1$cal.model
sk.oxyA1$full.model

sk.oxyA1.1<-skills(oxyA1.rhresp, target = .mean(5:8), model = "ols", calibration = "49%",
                 timespan = NULL)
plot(sk.oxyA1.1)
sk.oxyA1.1$RE;sk.oxyA1.1$CE
sk.oxyA1.1$cal.model

## d13cPin + d18O from young trees
z.pin.oxy.A1<-data.frame(apply(cbind(z.pin[1],z.oxy[1]),1,FUN = mean))

pin.oxy.A1.rhresp <- dcc(z.pin.oxy.A1,data.frame(Rh), method = "corr",selection = .range(1:9)+.mean(5:8))
plot(pin.oxy.A1.rhresp)
sk.pin.oxy.A1<-skills(pin.oxy.A1.rhresp, target = .mean(5:8), model = "ols", calibration = "-50%",
                  timespan = NULL)
plot(sk.pin.oxy.A1)
sk.pin.oxy.A1$full.model$rsquare

sk.pin.oxy.A1.1<-skills(pin.oxy.A1.rhresp, target = .mean(5:8), model = "ols", calibration = "49%",
                  timespan = NULL)
plot(sk.pin.oxy.A1.1)
sk.pin.oxy.A1.1$RE;sk.pin.oxy.A1.1$CE
sk.pin.oxy.A1.1$cal.model
sk.pin.oxy.A1.1$full.model
sk.pin.oxy.A1$full.model$rsquare


## d13cPin + d18O from old trees
z.pin.oxy.A3<-data.frame(apply(cbind(z.pin[3],z.oxy[3]),1,FUN = mean))

pin.oxy.A3.rhresp <- dcc(z.pin.oxy.A3,data.frame(Rh), method = "corr",selection = .range(1:9)+.mean(5:8))
plot(pin.oxy.A3.rhresp)
sk.pin.oxy.A3<-skills(pin.oxy.A3.rhresp, target = .mean(5:8), model = "ols", calibration = "-50%",
                      timespan = NULL)
plot(sk.pin.oxy.A3)
sk.pin.oxy.A3$full.model$rsquare

sk.pin.oxy.A3.1<-skills(pin.oxy.A3.rhresp, target = .mean(5:8), model = "ols", calibration = "49%",
                      timespan = NULL)
plot(sk.pin.oxy.A3.1)
sk.pin.oxy.A3.1$RE;sk.pin.oxy.A3.1$CE
sk.pin.oxy.A3.1$cal.model
sk.pin.oxy.A3.1$full.model$rsquare

##d13cPin from young trees + d18O from old trees
z.pinA1.oxy.A3<-data.frame(apply(cbind(z.pin[1],z.oxy[3]),1,FUN = mean))
pinA1.oxy.A3.rhresp <- dcc(z.pinA1.oxy.A3,data.frame(Rh), method = "corr",selection = .range(1:9)+.mean(5:8))
plot(pinA1.oxy.A3.rhresp)
sk.pinA1.oxy.A3<-skills(pinA1.oxy.A3.rhresp, target = .mean(5:8), model = "ols", calibration = "-50%",
                      timespan = NULL)
plot(sk.pinA1.oxy.A3)
sk.pinA1.oxy.A3$full.model$rsquare

sk.pinA1.oxy.A3.1<-skills(pinA1.oxy.A3.rhresp, target = .mean(5:8), model = "ols", calibration = "49%",
                      timespan = NULL)
plot(sk.pinA1.oxy.A3.1)
sk.pinA1.oxy.A3.1$RE;sk.pinA1.oxy.A3.1$CE
sk.pinA1.oxy.A3.1$cal.model
sk.pinA1.oxy.A3$full.model$rsquare

## d13cPin + d18O from All trees
z.pin.oxy.all<-data.frame(apply(cbind(z.pin,z.oxy),1,FUN = mean))

pin.oxy.all.rhresp <- dcc(z.pin.oxy.all,data.frame(Rh), method = "corr",selection = .range(1:9)+.mean(5:8))
plot(pin.oxy.all.rhresp)
sk.pin.oxy.all<-skills(pin.oxy.all.rhresp, target = .mean(5:8), model = "ols", calibration = "-50%",
                      timespan = NULL)
plot(sk.pin.oxy.all)
sk.pin.oxy.all$full.model$rsquare

sk.pin.oxy.all.1<-skills(pin.oxy.all.rhresp, target = .mean(5:8), model = "ols", calibration = "49%",
                      timespan = NULL)
plot(sk.pin.oxy.all.1)
sk.pin.oxy.all.1$RE;sk.pin.oxy.all.1$CE
sk.pin.oxy.all.1$cal.model
sk.pin.oxy.all.1$full.model$rsquare

##5.1.2 output the figure-------
##
#### pin reconstruction
#### Here, we use our "getseries" function to assign the reconstruction series
sk.pinallrecon<-ts(getseries(z.pin.ts,sk.pinall),
                        start = 1910,end=2010)
sk.pinA1recon<-ts(getseries(z.pin[1],sk.pinA1),
                        start = 1910,end=2010)
sk.pinA3recon<-ts(getseries(z.pin[3],sk.pinA3),
                        start = 1910,end=2010)

#### oxy reconstruction
sk.oxyallrecon<-ts(getseries(z.oxy.ts,sk.oxyall),
                        start = 1910,end=2010)
sk.oxyA1recon<-ts(getseries(z.oxy[1],sk.oxyA1),
                        start = 1910,end=2010)
sk.oxyA3recon<-ts(getseries(z.oxy[3],sk.oxyA3),
                        start = 1910,end=2010)

#### combine reconstruction
sk.pin.oxy.allrecon<-ts(getseries(z.pin.oxy.all,sk.pin.oxy.all),
                        start = 1910,end=2010)
sk.pin.oxy.A1recon<-ts(getseries(z.pin.oxy.A1,sk.pin.oxy.A1),
                        start = 1910,end=2010)
sk.pin.oxy.A3recon<-ts(getseries(z.pin.oxy.A3,sk.pin.oxy.A3),
                        start = 1910,end=2010)
sk.pinA1.oxy.A3recon<-ts(getseries(z.pinA1.oxy.A3,sk.pinA1.oxy.A3),
                        start = 1910,end=2010)

tiff("./plot/Figure 6. climate reconstruction comparisons-4 spline20.tiff",width=24, height=23,units="cm",res = 300)
par(family="serif",ps=13)
#par(mfrow=c(3,1),mgp=c(2.0,0.5,0),mar=c(3, 3, 0.5, 1) + 0.1,cex=1)
par(mfrow=c(3,1),mgp=c(2.0,0.5,0),mar=c(2, 3, 0, 0.5) + 0.1,cex=1)
plot(sk.pinall$years,sk.pinall$full$x,"l",xlim = c(1910,2010),
     ylim=c(48,86),xlab = " ", ylab = " ",las=1)

lines(window(sk.pinallrecon[,2],end=1954),col=2,lty=3)
lines(sk.pinall$cal.years,sk.pinall$pred.cal,"l",col=2)
lines(sk.pinall$ver.years,sk.pinall$pred.ver,"l",col=2,lty=2)
lines(splinesmoother2((sk.pinallrecon[,2]),20),
      col=2,lwd=2)

lines(window(sk.pinA1recon[,2],end=1954),col=4,lty=3)
lines(sk.pinA1$cal.years,sk.pinA1$pred.cal,"l",col=4)
lines(sk.pinA1$ver.years,sk.pinA1$pred.ver,"l",col=4,lty=2)
lines(splinesmoother2((sk.pinA1recon[,2]),20),
      col=4,lwd=2)

lines(window(sk.pinA3recon[,2],end=1954),col=6,lty=3)
lines(sk.pinA3$cal.years,sk.pinA3$pred.cal,"l",col=6)
lines(sk.pinA3$ver.years,sk.pinA3$pred.ver,"l",col=6,lty=2)
lines(splinesmoother2((sk.pinA3recon[,2]),20),
      col=6,lwd=2)


legend(1940,88,legend = c(expression(paste(delta ^13,"C-young")),
                          expression(paste(delta ^13,"C-old")),
                          expression(paste(delta ^13,"C-all"))
                          ),
                       col=c(4,6,2),lty=0,
                       text.col = c(4,6,2),bty = "n")

legend(1960,88,legend = c("Observation","Calibration (color)","Verification (color)","Recon before instrument"),
                       col=c(1,1,1,1),lty=c(1,1,2,3),
                       bty = "n")
text(1910,85,labels = "a",cex=2.5,font = 2)

plot(sk.oxyall$years,sk.oxyall$full$x,"l",
     xlim = c(1910,2010),ylim=c(48,86),
     xlab = " ", ylab = "May-August relative humidity (%)",las=1)
lines(window(sk.oxyallrecon[,2],end=1954),col=2,lty=3)
lines(sk.oxyall$cal.years,sk.oxyall$pred.cal,"l",col=2)
lines(sk.oxyall$ver.years,sk.oxyall$pred.ver,"l",col=2,lty=2)
lines(splinesmoother2((sk.oxyallrecon[,2]),20),
      col=2,lwd=2)

lines(window(sk.oxyA1recon[,2],end=1954),col=4,lty=3)
lines(sk.oxyA1$cal.years,sk.oxyA1$pred.cal,"l",col=4)
lines(sk.oxyA1$ver.years,sk.oxyA1$pred.ver,"l",col=4,lty=2)
lines(splinesmoother2((sk.oxyA1recon[,2]),20),
      col=4,lwd=2)

lines(window(sk.oxyA3recon[,2],end=1954),col=6,lty=3)
lines(sk.oxyA3$cal.years,sk.oxyA3$pred.cal,"l",col=6)
lines(sk.oxyA3$ver.years,sk.oxyA3$pred.ver,"l",col=6,lty=2)
lines(splinesmoother2((sk.oxyA3recon[,2]),20),
      col=6,lwd=2)

legend(1940,88,
       legend = c(expression(paste(delta ^18,"O-young")),
                  expression(paste(delta ^18,"O-old")),
                  expression(paste(delta ^18,"O-all"))                          ),
       col=c(1,1,1,1),lty=0,
       text.col = c(4,6,2),bty = "n")
legend(1960,88,
       legend = c("Observation","Calibration (color)","Verification (color)","Recon before instrument"),
       col=c(1,1,1,1),lty=c(1,1,2,3),
       text.col = c(1,1,1),bty = "n")
text(1910,84,labels = "b",cex=2.5,font=2)


par(mar=c(3, 3, 0, 0.5) + 0.1,cex=1)
plot(sk.pin.oxy.all$years,sk.pin.oxy.all$full$x,"l",
     xlim = c(1910,2010),ylim=c(48,86),
     xlab = "Year", ylab = " ",las=1)

lines(window(sk.pin.oxy.allrecon[,2],end=1954),col=2,lty=3)
lines(sk.pin.oxy.all$cal.years,sk.pin.oxy.all$pred.cal,"l",col=2)
lines(sk.pin.oxy.all$ver.years,sk.pin.oxy.all$pred.ver,"l",col=2,lty=2)
lines(splinesmoother2((sk.pin.oxy.allrecon[,2]),20),
      col=2,lwd=2)
# boxplot(window(sk.pin.oxy.allrecon[,2],end=1954),col=2)

lines(window(sk.pin.oxy.A1recon[,2],end=1954),col=4,lty=3)
lines(sk.pin.oxy.A1$cal.years,sk.pin.oxy.A1$pred.cal,"l",col=4)
lines(sk.pin.oxy.A1$ver.years,sk.pin.oxy.A1$pred.ver,"l",col=4,lty=2)
lines(splinesmoother2((sk.pin.oxy.A1recon[,2]),20),
      col=4,lwd=2)

lines(window(sk.pin.oxy.A3recon[,2],end=1954),col=6,lty=3)
lines(sk.pin.oxy.A3$cal.years,sk.pin.oxy.A3$pred.cal,"l",col=6)
lines(sk.pin.oxy.A3$ver.years,sk.pin.oxy.A3$pred.ver,"l",col=6,lty=2)
lines(splinesmoother2((sk.pin.oxy.A3recon[,2]),20),
      col=6,lwd=2)

lines(window(sk.pinA1.oxy.A3recon[,2],end=1954),col=3,lty=3)
lines(sk.pinA1.oxy.A3$cal.years,sk.pinA1.oxy.A3$pred.cal,"l",col=3)
lines(sk.pinA1.oxy.A3$ver.years,sk.pinA1.oxy.A3$pred.ver,"l",col=3,lty=2)
lines(splinesmoother2((sk.pinA1.oxy.A3recon[,2]),20),
      col=3,lwd=2)

legend(1930,88,legend = c(expression(paste(delta ^13,"C-young + ", 
                                           delta ^18,"O-young")),
                          expression(paste(delta ^13,"C-old + ", 
                                           delta ^18,"O-old")),
                          expression(paste(delta ^13,"C-young + ", 
                                           delta ^18,"O-old")),
                          expression(paste(delta ^13,"C-all + ", 
                                           delta ^18,"O-all"))),
       col=c(4,6,3,2),lty=0,text.col = c(4,6,3,2),bty = "n")
legend(1960,88,legend = c("Observation","Calibration (color)","Verification (color)","Recon before instrument"), 
       col=1,lty=c(1,1,2,3),text.col = 1,bty = "n")
text(1910,85,labels = "c",cex=2.5,font=2)
dev.off()

##5.1.3 scatter plot and penalty analysis------
scatt_data.c<- rbind(sk.pinA1$full,sk.pinA3$full,sk.pinall$full)
scatt_data.o<- rbind(sk.oxyA1$full,sk.oxyA3$full,sk.oxyall$full)
scatt_data.co<- rbind(sk.pin.oxy.A1$full,sk.pin.oxy.A3$full,
                      sk.pin.oxy.all$full,sk.pinA1.oxy.A3$full)

summary(sk.oxyA1$full.model)

group1<-c(rep("young",57),
              rep("old",57),
              rep("all",57))
group2<-c(group1,rep("C-young & O-old",57))

scatt_data.c$group<-group1
scatt_data.o$group<-group1
scatt_data.co$group<-group2

summary(glm(x~poly(y,2),data=subset(scatt_data.c,group=="old")))
## here, 1 - (Residual Deviance/Null Deviance) will give the R2.
summary(lm(x~poly(y,2),data=subset(scatt_data.c,group=="all")))

summary(glm(x~poly(y,2),data=subset(scatt_data.c,group=="old")))
## here, 1 - (Residual Deviance/Null Deviance) will give the R2.
summary(lm(x~poly(y,2),data=subset(scatt_data.o,group=="old")))

summary(glm(x~poly(y,2),data=subset(scatt_data.c,group=="old")))
## here, 1 - (Residual Deviance/Null Deviance) will give the R2.
summary(lm(x~poly(y,2),data=subset(scatt_data.co,
                                   group=="C-young & O-old")))

summary(glm(x~y,data=subset(scatt_data.c,group=="all")))
## here, 1 - (Residual Deviance/Null Deviance) will give the R2.
fit<-lm(x~y,data=subset(scatt_data.co,group=="C-young & O-old"))
summary(fit)
BIC(fit)
AIC(fit)
sqrt(mean(fit$residuals^2))# calculate RMSE
RMSE(fit)

scatplot1<-ggplot(data=scatt_data.c,aes(x=y,y=x,col=group,shape=group))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE, fullrange=TRUE)+
  scale_color_manual(values=c(2,6,4),labels = c("All","Old","Young"))+
  guides(shape=FALSE)+
  mythemeplot()+
  ylab("May-August relative humidity (%)")+
  xlab("Z-score")+
  theme(legend.position = c(0.7,0.85),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
    labs(title = c(expression(paste(delta ^"13","C "))))+
    annotate("text", x = -1, y=rev(c(50,52,54)), 
             label = c(expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.59")),
                       expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.19")),
                       expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.60"))),
             col=c(2,6,4))

scatplot2<-ggplot(data=scatt_data.o,aes(x=y,y=x,col=group,shape=group))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE, fullrange=TRUE)+
  scale_color_manual(values=c(2,6,4),labels = c("All","Old","Young"))+
  guides(shape=FALSE)+
  mythemeplot()+
  ylab("May-August relative humidity (%)")+
  xlab("Z-score")+
  theme(legend.position = c(0.7,0.85),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
    labs(title = c(expression(paste(delta ^"18"," O"))))+
  annotate("text", x = -1.5, y=rev(c(50,52,54)), 
             label = c(expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.37")),
                       expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.34")),
                       expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.35"))),
             col=c(2,6,4))

scatplot3<- ggplot(data=scatt_data.co,aes(x=y,y=x,col=group,shape=group))+
  geom_point()+
  geom_smooth(method=lm,se=FALSE, fullrange=TRUE)+
  scale_color_manual(values=c(2,3,6,4),
                     labels = c("All","C-young & O-old","Old","Young"))+
  guides(shape=FALSE)+
  mythemeplot()+
  ylab("May-August relative humidity (%)")+
  xlab("Z-score")+
  theme(legend.position = c(0.7,0.85),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
    labs(title = c(expression(paste(delta ^"13","C + ",delta ^"18"," O"))))+
  annotate("text", x = -1, y=rev(c(50,52,54,56)), 
             label = c(expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.60")),
                       expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.60")),
                       expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.41")),
                       expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.60"))),
             col=c(2,3,6,4))

tiff("./plot/Figure S4 scattle plot for the regression.tiff", width = 25,height = 10,
     units = "cm",compression = "lzw",res = 300,
     bg = "white",family = "serif")
ggarrange( scatplot1,scatplot2,scatplot3,
          nrow=1,ncol=3,
          label.x=0.90,label.y=0.92,
          labels=c("a","b","c"),
          font.label = list(size=22,family="serif"))
dev.off()

scatplotc1<-ggplot(data=scatt_data.c,aes(x=y,y=x,col=group,shape=group))+
  geom_point()+
  #geom_smooth(method=lm,se=FALSE, fullrange=TRUE)+
  geom_smooth(method="gam",formula = y ~ poly(x,2),
              se=FALSE)+
  scale_color_manual(values=c(2,6,4),labels = c("All","Old","Young"))+
  guides(shape=FALSE)+
  mythemeplot()+
  ylab("May-August relative humidity (%)")+
  xlab("Z-score")+
  theme(legend.position = c(0.7,0.85),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
    labs(title = c(expression(paste(delta ^"13","C "))))+
    annotate("text", x = -1, y=rev(c(50,52,54)), 
             label = c(expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.59")),
                       expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.18")),
                       expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.64"))),
             col=c(2,6,4))

scatplotc2<-ggplot(data=scatt_data.o,aes(x=y,y=x,col=group,shape=group))+
  geom_point()+
  #geom_smooth(method=lm,se=FALSE, fullrange=TRUE)+
  geom_smooth(method="lm",formula = y ~ poly(x,2),
              se=FALSE)+
  scale_color_manual(values=c(2,6,4),
                     labels = c("All","Old","Young"))+
  guides(shape=FALSE)+
  mythemeplot()+
  ylab("May-August relative humidity (%)")+
  xlab("Z-score")+
  theme(legend.position = c(0.7,0.85),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
    labs(title = c(expression(paste(delta ^"18"," O"))))+
  annotate("text", x = -1.5, y=rev(c(50,52,54)), 
             label = c(expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.35")),
                       expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.33")),
                       expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.37"))),
             col=c(2,6,4))

scatplotc3<- ggplot(data=scatt_data.co,aes(x=y,y=x,col=group,shape=group))+
  geom_point()+
  #geom_smooth(method=lm,se=FALSE, fullrange=TRUE)+
  geom_smooth(method="lm",formula = y ~ poly(x,2),
              se=FALSE)+
  scale_color_manual(values=c(2,3,6,4),
                     labels = c("All","C-young & O-old","Old","Young"))+
  guides(shape=FALSE)+
  mythemeplot()+
  ylab("May-August relative humidity (%)")+
  xlab("Z-score")+
  theme(legend.position = c(0.7,0.85),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))+
    labs(title = c(expression(paste(delta ^"13","C + ",delta ^"18"," O"))))+
  annotate("text", x = -1, y=rev(c(50,52,54,56)), 
             label = c(expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.59")),
                       expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.61")),
                       expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.40")),
                       expression(paste(italic(R)^"2",
                                        ""["adj"],"= 0.61"))),
             col=c(2,3,6,4))

## We used non-linear regression and we get the similar R2.
tiff("./plot/Figure S5 scattle plot for the regression.tiff", width = 25,height = 10,
     units = "cm",compression = "lzw",res = 300,
     bg = "white",family = "serif")
ggarrange( scatplotc1,scatplotc2,scatplotc3,
          nrow=1,ncol=3,
          label.x=0.90,label.y=0.92,
          labels=c("a","b","c"),
          font.label = list(size=22,family="serif"))
dev.off()


##5.2 make the table for the statistics-------
r.calibration<-c(sk.pinA1$r.cal, sk.pinall$r.cal,
                 sk.oxyA3$r.cal, sk.oxyall$r.cal, 
                 sk.pin.oxy.A1$r.cal,sk.pin.oxy.A3$r.cal,sk.pin.oxy.all$r.cal)
R2_calibration<-c(sk.pinA1$cal.model$rsquare, sk.pinall$cal.model$rsquare,
                  sk.oxyA3$cal.model$rsquare, sk.oxyall$cal.model$rsquare, 
                  sk.pin.oxy.A1$cal.model$rsquare,sk.pin.oxy.A3$cal.model$rsquare,sk.pin.oxy.all$cal.model$rsquare)
r.full<-c(sk.pinA1$r.ful, sk.pinall$r.ful,
          sk.oxyA3$r.ful, sk.oxyall$r.ful, 
          sk.pin.oxy.A1$r.ful,sk.pin.oxy.A3$r.ful,sk.pin.oxy.all$r.ful)
RE<-c(sk.pinA1$RE, sk.pinall$RE,
      sk.oxyA3$RE, sk.oxyall$RE, 
      sk.pin.oxy.A1$RE,sk.pin.oxy.A3$RE,sk.pin.oxy.all$RE)
CE<-c(sk.pinA1$CE, sk.pinall$CE,
      sk.oxyA3$CE, sk.oxyall$CE, 
      sk.pin.oxy.A1$CE,sk.pin.oxy.A3$CE,sk.pin.oxy.all$CE)
R2<-c(sk.pinA1$full.model$rsquare, sk.pinall$full.model$rsquare,
      sk.oxyA3$full.model$rsquare, sk.oxyall$full.model$rsquare, 
      sk.pin.oxy.A1$full.model$rsquare,sk.pin.oxy.A3$full.model$rsquare,sk.pin.oxy.all$full.model$rsquare)

statis<-data.frame(r.calibration*-1,R2_calibration,RE,CE,r.full*-1,R2)
statis<-round(statis,digits = 2)


## 6. Calculate the Tree-ring width---------
AWL<-read.rwl("./width/AWL.RWL")
AWL1<-read.rwl("./width/AWL1.RWL")
AWL2<-read.rwl("./width/AWL2.RWL")
AWL3<-read.rwl("./width/AWL3.RWL")
AWLall<-read.rwl("./width/AWL-all.RWL")
colnames(AWL1)
TRW.all<-rowMeans(AWLall,na.rm=TRUE)
plot(AWLall, plot.type="spag")

TRW.allrwi<-detrend(rwl = AWLall, method = "ModNegExp")
## bulid a mean rwi chronology
TRW.allcrn <- chron(TRW.allrwi, prefix = "CAM")

young.rwi <- subset(TRW.allrwi,
                    select=c('AWL01','AWL05','AWL06',
                             'AWL3-1A','AWL3-2B'))
young.rwimean <-apply(young.rwi,1,FUN=mean,na.rm=TRUE)

mid.rwi <- subset(TRW.allrwi,
              select=c('AWL1-13B','AWL1-17','AWL1-19','AWL1-22A'))#AWL1-13B,AWL1-17£¬19, 22A

mid.rwimean <-apply(mid.rwi,1,FUN=mean,na.rm=TRUE)

old.rwi <- subset(TRW.allrwi,
              select=c('AWL1-1A','AWL1-1B','AWL1-2A',
                       'AWL1-4A','AWL1-6A',
                            'AWL1-14A','AWL1-15A',
                       'AWL1-25A','AWL1-24A','AWL1-27B'))
old.rwimean<-apply(old.rwi,1,FUN=mean,na.rm=TRUE)
plot(c(1540:2010),old.rwimean,"l")

TRW.ids <- read.ids(AWLall, stc = c(4, 2, 3))

foo <- rwi.stats.running(TRW.allrwi, TRW.ids,
 window.length = 30)

TRW.allsss <- sss(TRW.allrwi)
yrs <- time(TRW.allcrn)
plot(yrs, TRW.allsss , type = "l", xlab = "", ylab = "",
 axes = FALSE, col = "blue")


plot(subset(TRW.allcrn,samp.depth>10),
     add.spline=T,nyrs=30)

tiff("./plot/Figure S1-2 tree width for the carbon and oxyge isotpe analysis.tiff",width = 16, height = 10,
     units = "cm",compression = "lzw",bg="white",res = 300,
     family = "serif")
par(mar=c(2.5, 2.5, 0.5, 3) + 0.1,mgp=c(1.5,0.5,0))
plot(yrs,TRW.allcrn[,1], xlim=c(1540,2010), 
     ylim = c(0.3,1.72),
     type = "l",
     xlab = "Year", ylab = "RWI (Ring width index)")

 lines(c(1540:2010),old.rwimean,"l",col="#619CFF",lwd=0.7)
 lines(c(1540:2010),mid.rwimean,"l",col="#00BA38",lwd=0.7)
 lines(c(1540:2010),young.rwimean,"l",col="#F8766D",lwd=0.7)
 lines(yrs, ffcsaps(TRW.allcrn[,1], nyrs = 30), col = "red", lwd = 2)
#axis(1); axis(2); axis(3);

legend(1550,0.5,
       legend = c("Master chronology","30-year spline",
                  "Young","Middle","Old","SSS"),
       col=c(1,"red","#F8766D","#00BA38", "#619CFF","#E69F00"),
       lty=1,bty = "n",bg="transparent",
       ncol=3)
# par(new=TRUE)
# plot(foo$end.year, foo$eps,xlim=c(1540,2010),type="b",col="blue",axes=FALSE,ylim=c(0.2,1),ann=FALSE)
par(new=TRUE)
plot(yrs, TRW.allsss,xlim=c(1540,2010),type="l",col="#E69F00",axes=FALSE,ylim=c(0.2,1),ann=FALSE)
abline(h=0.85,col="#E69F00",lty="dashed")
axis(4, at = pretty(TRW.allsss),
     col="#E69F00",
     col.axis="#E69F00")
mtext(4,text="SSS (subsample signal strength)",col="#E69F00",line=2)
dev.off()



young <- subset(AWL,select=c('AWL01','AWL05','AWL06'))
young2 <- subset(AWL3,select=c('AWL3-1A','AWL3-2B'))
young.trw <- cbind(subset(young,rownames(young)>1909),subset(young2,rownames(young2)>1909))
summary(apply(young.trw,1,FUN=mean,na.rm=TRUE))

mid <- subset(AWL1,select=c('AWL1-13B','AWL1-17','AWL1-19','AWL1-22A'))#AWL1-13B,AWL1-17£¬19, 22A
mid.trw <- subset(mid,rownames(mid)>1909)
summary(apply(mid.trw,1,FUN=mean,na.rm=TRUE))

old <- subset(AWL1,select=c('AWL1-1A','AWL1-1B','AWL1-2A','AWL1-4A','AWL1-6A',
                            'AWL1-14A','AWL1-15A','AWL1-25','AWL1-24','AWL1-27B'))
old.trw <- subset(old,rownames(old)>1909)
summary(apply(old.trw,1,FUN=mean,na.rm=TRUE))

core <- read.xlsx("./width/cores.xlsx")

age.plot<-ggplot(core)+
  geom_segment(aes(y=x,x=start,yend=x,xend=end,col=Group))+
  scale_x_continuous(expand = c(0.005, 0.005))+
  scale_y_continuous(expand = c(0.005, 0.005))+
  #scale_color_manual(values=c("darkgreen", "#E69F00", "#56B4E9"))+
  scale_color_manual(values=c("Middle"="#00BA38",
                              "Old"= "#619CFF",
                              "Young"="#F8766D"))+
  geom_hline(yintercept = c(2,3),col=1,lty=2)+
  annotate(geom = "rect",xmin = 1910,xmax = 2010,ymin = 1.7,ymax = 4,alpha=0.3)+
  annotate(geom = "rect",xmin = 1920,xmax = 2010,ymin = 1.3,ymax = 1.7,alpha=0.3)+
  annotate(geom = "rect",xmin = 1910,xmax = 2010,ymin = 0.8,ymax = 1.3,alpha=0.3)+
  mythemeplot()+ theme(legend.position = c(0.1,0.15))+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  xlab("Year")+ylab("")

## output
tiff("./plot/Figure 1 tree age for the carbon and oxyge isotpe analysis.tiff",width = 16, height = 8,
     units = "cm",compression = "lzw",bg="white",res = 300)
age.plot
dev.off()


