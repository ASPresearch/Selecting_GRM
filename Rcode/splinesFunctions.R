Qnn<-function(x, q) # q: quantil en % (1-100)
{
  x<-as.numeric(x)
  Q<-quantile(x, q/100, na.rm=TRUE)
  return(Q)
}

#QInf i QSup determinen quin quartil volem calcular, metInf i metSup quin  el llindar dels quartils #que volem utilitzar (QInf<metInf, QSup>metSup)

Initialselection<-function(allData, metilData, corrData, 
                           QInf=25, metInf, QSup=75, metSup, 
                           Adjust=FALSE, pAdj)
{
  QSupDades<-apply(metilData,1,Qnn, QSup) ###de les metilacions
  QInfDades<-apply(metilData,1,Qnn, QInf)  ###de les metilacions
  pv<-corrData[,5]
  if (Adjust)   pv<-corrData[,7]
  selCond  <- corrData[ ,2] < 0 & 
    pv < pAdj & 
    QSupDades > metSup & 
    QInfDades < metInf
  Cors <- corrData[selCond,]
  DadesSel <- allData[selCond, ]  
  return(DadesSel)
}


calculaSplines<-function(mat)
{
  Qmet <- quantile(mat$met, probs=c(0.25,0.5,0.75))
  reg<-lm(expr~bs(met,knots=Qmet, degree=2), data=mat)
  summ<-summary(reg)
  return(summ$coef)
}


#Scatterplot amb la corva spline dibuixada

plotWithSplines<-function(mat,titleText)
{
  x<-mat[,1]
  y<-mat[,2]
  maxy<-max(y)
  miny<-min(y)
  plot(x,y,  xlim=c(0,1),ylim=c(miny, maxy),main=titleText)
  Qmet<-quantile(x, probs=c(0.25,0.5,0.75))
  reg<-lm(expr~bs(met,knots=Qmet, degree=2), data=mat)
  #reg<-lm(expr~bs(met, df=5,intercept=FALSE), data=mat)
  minval<-min(x)
  maxval<-max(x)
  u<-seq(minval,maxval,by=.1)
  B<-data.frame(met=u)
  Y<-predict(reg,newdata=B)
  lines(u,Y,lwd=2,col="red")
}


#Grafic dels elements d'un cluster

PlotAGroup<-function(nameFile, Members, ListGenes)
{
  if (nameFile !="") pdf(nameFile)
  for (IdGene in Members ) 
  {
    myTitle<-paste("Gene",IdGene, sep=" ")
    mat<-ListGenes[[IdGene]]
    plotWithSplines(mat,myTitle)
  }
  if (nameFile !="") dev.off()
}

