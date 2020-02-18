 ##2.3 calculate the bootstrap corelation between chronlogy and climate and test the sensitivity between groups
  library(pgirmess)## this package is used multi-comparison
  require(boot)
  
  climres_com<-function(sty,edy,oxy,Rh,Brep,stm,edm,stmp,edmp){
    ## sty is the start year, edy is the end of the year for the climate data, oxy is chronology with year as rownames
    ## Rh is the climatic data and the first conlumn is year
    ## Brep is the times of the bootstrap
    ## stm,edm,stmp,edpm is the start month and end month for the current year and the previous year!!
    ## here, we cannot estimate the seasonal climate data in the function
  
  sty;edy; stm; edm;Brep
  RP<-NULL
  KS<-NULL
  for (odmonth in c(stm:edm)) {
    chronclim.data<-cbind(subset(oxy,rownames(oxy)>sty-1)[1],subset(Rh,Rh[1]< edy+1)[odmonth])
    chronclim.data1<-cbind(subset(oxy,rownames(oxy)>sty-1)[2],subset(Rh,Rh[1]< edy+1)[odmonth])
    chronclim.data2<-cbind(subset(oxy,rownames(oxy)>sty-1)[3],subset(Rh,Rh[1]< edy+1)[odmonth])
    
    # Bootstrap the Pearson correlation coefficient-
    n = length(sty:edy)
    #Brep = 1000
    pearson <- function(d,i=c(1:n)){
      return(cor(d[i,1],d[i,2]))
    }
    bootcorr <- boot(data=chronclim.data,statistic=pearson,R=Brep)
    bootcorr1<- boot(data=chronclim.data1,statistic=pearson,R=Brep)
    bootcorr2<- boot(data=chronclim.data2,statistic=pearson,R=Brep)
    #boot.ci(bootcorr,conf=.95)
    ### check the distribution of the boostrap correlations
    ##par(mfrow=c(2,1))
    ##hist(bootcorr$t,main="Bootstrap Pearson Sample Correlation Coefficients")
    ###plot(ecdf(bootcorr$t),main="ECDF of Bootstrap Correlation Coefficients")
    bootcoorall<-as.data.frame(rbind(bootcorr$t,bootcorr1$t,bootcorr2$t))
    bootcoorall$group<-c(rep("A1",1000),rep("A2",1000),rep("A3",1000))
    #head(bootcoorall)
    
    ## test the significant difference of the correlations----
    ## this is based on the Mann-Whitney-Wilcoxon Test, pairs comparison
    #wilcox.test(bootcoorall$V1~bootcoorall$group)
    PT<-pairwise.wilcox.test(bootcoorall$V1, bootcoorall$group)
    PT
    PT$p.value
    RP<-rbind(RP,PT$p.value)
    
    
    ## Multiple comparison test after Kruskal-Wallis
    
    #boxplot(V1 ~ group, data=bootcoorall)
    KT<-kruskalmc(bootcoorall$V1, bootcoorall$group)
    KT99<-kruskalmc(bootcoorall$V1, bootcoorall$group, probs=0.01)
    #kruskalmc(bootcoorall$V1, bootcoorall$group, cont="one-tailed")
    #kruskalmc(bootcoorall$V1, bootcoorall$group, cont="two-tailed")
    KS<-rbind(KS,KT99$dif.com)
  }
  
  ## for previous year
  RP.p<-NULL
  KS.p<-NULL
  for (odmonth in c(stmp:edmp)) {
    chronclim.data<-cbind(subset(oxy,rownames(oxy)>sty)[1],
                          subset(Rh,Rh[1]< edy)[odmonth])
    chronclim.data1<-cbind(subset(oxy,rownames(oxy)>sty)[2],
                           subset(Rh,Rh[1]< edy)[odmonth])
    chronclim.data2<-cbind(subset(oxy,rownames(oxy)>sty)[3],
                           subset(Rh,Rh[1]< edy)[odmonth])
    
    # Bootstrap the Pearson correlation coefficient-
    #library(boot)
    n = length(sty:edy)
    pearson <- function(d,i=c(1:n)){
      return(cor(d[i,1],d[i,2]))
    }
    bootcorr <- boot(data=chronclim.data,statistic=pearson,R=Brep)
    bootcorr1<- boot(data=chronclim.data1,statistic=pearson,R=Brep)
    bootcorr2<- boot(data=chronclim.data2,statistic=pearson,R=Brep)
    #boot.ci(bootcorr,conf=.95)
    ### check the distribution of the boostrap correlations
    ##par(mfrow=c(2,1))
    ##hist(bootcorr$t,main="Bootstrap Pearson Sample Correlation Coefficients")
    ###plot(ecdf(bootcorr$t),main="ECDF of Bootstrap Correlation Coefficients")
    bootcoorall<-as.data.frame(rbind(bootcorr$t,bootcorr1$t,bootcorr2$t))
    bootcoorall$group<-c(rep("A1",1000),rep("A2",1000),rep("A3",1000))
    #head(bootcoorall)
    
    ## test the significant difference of the correlations----
    ## this is based on the Mann-Whitney-Wilcoxon Test, pairs comparison
    #wilcox.test(bootcoorall$V1~bootcoorall$group)
    PT.p<-pairwise.wilcox.test(bootcoorall$V1, bootcoorall$group)
    PT.p
    PT.p$p.value
    RP.p<-rbind(RP.p,PT.p$p.value)
    
    
    ## Multiple comparison test after Kruskal-Wallis
    KT.p<-kruskalmc(bootcoorall$V1, bootcoorall$group)
    KT99.p<-kruskalmc(bootcoorall$V1, bootcoorall$group, probs=0.01)
    #kruskalmc(bootcoorall$V1, bootcoorall$group, cont="one-tailed")
    #kruskalmc(bootcoorall$V1, bootcoorall$group, cont="two-tailed")
    KS.p<-rbind(KS.p,KT99.p$dif.com)
  }
  return(list(RP.p,KS.p,RP,KS))
  }


####get the reconstruction sereies
getseries<-function(x,modelcoef){
  #here x should be [year, proxy]??modelcoef is the skill based on treeclim.
  series<-x[,1]*modelcoef$coef.full[2]+modelcoef$coef.full[1]
  x$series<-series
  return(x)
}


##########################################################################################
#function splinesmoother2
#
#as for splinesmoother, but deals with singlecolumns
##########################################################################################

splinesmoother2 <- function(x,smoothing) {
  
  spline.p.for.smooth.Pspline<-function(year)  {
    p <- .5/((6*(cos(2*pi/year)-1)^2/(cos(2*pi/year)+2)))
    return(p)
  }
  
  
  ifelse ((is.null(ncol(x))),x2 <- ts.union(x,x), x2 <- x)
  
  fyarray <- start(x2)[1]
  lyarray <- end(x2)[1]
  begin <- fy(x2)
  end <- ly(x2)
  
  smoothedarray <- x2
  
  spline.p <- spline.p.for.smooth.Pspline(smoothing)
  
  
  for (i in 1:ncol(x2)) {
    temp <- smooth.Pspline((begin[i]:end[i]),x2[(begin[i]-fyarray+1):(end[i]-fyarray+1),i],spar=spline.p,method=1)
    smoothedarray[(begin[i]-fyarray+1):(end[i]-fyarray+1),i]  <- temp$ysmth
  }
  
  if (is.null(ncol(x))) {smoothedarray <- smoothedarray[,-1]}
  return(smoothedarray)
}




