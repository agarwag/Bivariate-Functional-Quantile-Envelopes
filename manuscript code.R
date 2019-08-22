## Code for Multivariate Functional Quantile Envelopes with Application to Radiosonde Wind Data ##

##Reading radiosonde wind data
library(quantreg)
setwd("/Users/agarwag/Desktop/KAUST/project 2/radiosonde wind")
load(file="Denver.station.RData",verbose=T)

str(Denver.station)
Denver.station=as.data.frame(Denver.station)
head(Denver.station);tail(Denver.station);
dim(Denver.station)
attach(Denver.station)


#important functions

eli=function(qln) 
{
  dirs=dim(qln)[1]
  clnid=1
  cln=1
  stop=FALSE
  filn=c(1:dirs)
  ss=0
  sln=0
  while (!stop)
  {
    abc=dertln(filn,clnid)
    a=qln[abc[1],]
    b=qln[abc[2],]
    c=qln[abc[3],]
    tt=actln(a,b,c)
    if (ss==0&tt==1) sln=filn[clnid]
    if (tt*ss==1&cln==sln) stop=TRUE
    ss=tt
    m=length(filn)
    if (tt==1)
    {
      if (clnid==m) {cln=filn[1] 
      clnid=1} else 
      {cln=filn[clnid+1] 
      clnid=clnid+1}
    } else { 
      clnidt=clnid
      if (clnid==1) {cln=filn[m] 
      clnid=m-1} else 
      {cln=filn[clnid-1] 
      clnid=clnid-1}
      filn=filn[-clnidt]
    }
  } 
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts[i,]=inln(actlns[i,],actlns[j,])
  }
  z=list(actlns=actlns,intpts=intpts)
  z                 
}

inln <- function(a,b) ## compute intersection of two lines
{
  inln=-c(b[2]*a[3]-a[2]*b[3],-b[1]*a[3]+a[1]*b[3])/(a[1]*b[2]-a[2]*b[1])
  inln
}

actln <- function(a,b,c) ## tell if a line is active
{
  aa=inln(a,b)
  cc=inln(a,c)
  if (sum((aa-cc)*c[1:2])>0) actln=1 else actln=0
  actln
}

dertln <- function(x,y)  ## determine lines
{
  k=length(x)
  
  if (y==1) bln=x[k] else bln=x[y-1]
  
  if (y==k) aln=x[1] else aln=x[y+1]
  
  cln=x[y]
  
  c(bln,cln,aln)
}

mqli <- function(x, prob, dirs=100) ## generate directions
{
  ang = seq(-pi,pi,length=dirs+1)
  dir = cbind(cos(ang[1:dirs]),sin(ang[1:dirs]))
  cbind(dir,-apply(dir,1,function(u) quantile(x %*% u, prob, type=1)))
}

abcline <- function(a,b,c,...) ## plot lines in terms of coefficients
{
  if (b == 0) abline(v=-c/a,...) else  abline(-c/b,-a/b,...)
}

qdir <- function(dirs, x=NULL)
{
  ang = seq(0,2*pi,length=dirs+1)
  ang = ang[1:dirs]-ang[1]/2-ang[dirs]/2
  cbind(cos(ang),sin(ang))
}


linecolor = 'black'
cntcolor = 'grey30'
pointcolor = 'grey60'
cntpbs=c(0.0002,0.0001,0.00001,0.000001)
dirs=1000
clrs=rep(linecolor,length(cntpbs))
xlab='Horizontal wind speed'
ylab='Vertical wind speed'


complete=complete.cases(u70,v70,u100,v100,u200,v200,u250,v250,u300,v300,u400,v400,u500,v500,u700,v700)

horz70=u70; horz100=u100; horz200=u200; horz250=u250; horz300=u300;horz400=u400; horz500=u500; horz700=u700; horz850=u850; horz1000=u1000
vert70=v70; vert100=v100; vert200=v200; vert250=v250; vert300=v300;vert400=v400; vert500=v500; vert700=v700; vert850=v850; vert1000=v1000 

xx70=horz70[complete]; xx100=horz100[complete]; xx200=horz200[complete]; xx250=horz250[complete]; xx300=horz300[complete]; xx400=horz400[complete]; xx500=horz500[complete];xx700=horz700[complete]
yy70=vert70[complete]; yy100=vert100[complete]; yy200=vert200[complete]; yy250=vert250[complete]; yy300=vert300[complete]; yy400=vert400[complete]; yy500=vert500[complete];yy700=vert700[complete]


######################## Figure 1 #######################################
##extreme quantile envelopes

library(evir)

qd = qdir(dirs) 

nextremes = round(sum(complete)-sum(complete)*0.95) 
## fit gpd on 5% extremes (573)
## riskmeasures gives estimate of cntpbs^th directional quantile
qln70 = cbind(qd,-apply(qd,1,
                        function(u)
                          -riskmeasures(gpd(cbind(xx70,yy70) %*% -u, nextremes=573),1-cntpbs[2])[2] ))

qln100 = cbind(qd,-apply(qd,1,
                         function(u)
                           -riskmeasures(gpd(cbind(xx100,yy100) %*% -u, nextremes=573),1-cntpbs[2])[2] ))
qln200 = cbind(qd,-apply(qd,1,
                         function(u)
                           -riskmeasures(gpd(cbind(xx200,yy200) %*% -u, nextremes=573),1-cntpbs[2])[2] ))
qln250 = cbind(qd,-apply(qd,1,
                         function(u)
                           -riskmeasures(gpd(cbind(xx250,yy250) %*% -u, nextremes=573),1-cntpbs[2])[2] ))
qln300 = cbind(qd,-apply(qd,1,
                         function(u)
                           -riskmeasures(gpd(cbind(xx300,yy300) %*% -u, nextremes=573),1-cntpbs[2])[2] ))
qln400 = cbind(qd,-apply(qd,1,
                         function(u)
                           -riskmeasures(gpd(cbind(xx400,yy400) %*% -u, nextremes=573),1-cntpbs[2])[2] ))
qln500 = cbind(qd,-apply(qd,1,
                         function(u)
                           -riskmeasures(gpd(cbind(xx500,yy500) %*% -u, nextremes=573),1-cntpbs[2])[2] ))
qln700 = cbind(qd,-apply(qd,1,
                         function(u)
                           -riskmeasures(gpd(cbind(xx700,yy700) %*% -u, nextremes=573),1-cntpbs[2])[2] ))

cnt70 = eli(qln70); cnt100 = eli(qln100); cnt200 = eli(qln200); cnt250 = eli(qln250); cnt300 = eli(qln300); cnt400 = eli(qln400); cnt500 = eli(qln500); cnt700 = eli(qln700)

library(spatstat)
df70=data.frame(xx70,yy70,year[complete],julian[complete],hour[complete]); df100=data.frame(xx100,yy100,year[complete],julian[complete],hour[complete]);df200=data.frame(xx200,yy200,year[complete],julian[complete],hour[complete]);df250=data.frame(xx250,yy250,year[complete],julian[complete],hour[complete]);df250=data.frame(xx250,yy250,year[complete],julian[complete],hour[complete]);df300=data.frame(xx300,yy300,year[complete],julian[complete],hour[complete]);df400=data.frame(xx400,yy400,year[complete],julian[complete],hour[complete]); df500=data.frame(xx500,yy500,year[complete],julian[complete],hour[complete]);df700=data.frame(xx700,yy700,year[complete],julian[complete],hour[complete])
bound70 <- owin(poly=data.frame(x=cnt70$intpts[,1], y=cnt70$intpts[,2])); bound100 <- owin(poly=data.frame(x=cnt100$intpts[,1], y=cnt100$intpts[,2])); bound200 <- owin(poly=data.frame(x=cnt200$intpts[,1], y=cnt200$intpts[,2])); bound250 <- owin(poly=data.frame(x=cnt250$intpts[,1], y=cnt250$intpts[,2])); bound300 <- owin(poly=data.frame(x=cnt300$intpts[,1], y=cnt300$intpts[,2])); bound400 <- owin(poly=data.frame(x=cnt400$intpts[,1], y=cnt400$intpts[,2])); bound500 <- owin(poly=data.frame(x=cnt500$intpts[,1], y=cnt500$intpts[,2])); bound700 <- owin(poly=data.frame(x=cnt700$intpts[,1], y=cnt700$intpts[,2]))
isin70<-inside.owin(x=xx70,y=yy70,w=bound70); isin100<-inside.owin(x=xx100,y=yy100,w=bound100); isin200<-inside.owin(x=xx200,y=yy200,w=bound200); isin250<-inside.owin(x=xx250,y=yy250,w=bound250); isin300<-inside.owin(x=xx300,y=yy300,w=bound300); isin400<-inside.owin(x=xx400,y=yy400,w=bound400); isin500<-inside.owin(x=xx500,y=yy500,w=bound500); isin700<-inside.owin(x=xx700,y=yy700,w=bound700)
point_out70 <- df70[!isin70,]; point_out100 <- df100[!isin100,]; point_out200 <- df200[!isin200,]; point_out250 <- df250[!isin250,]; point_out300 <- df300[!isin300,]; point_out400 <- df400[!isin400,]; point_out500 <- df500[!isin500,]; point_out700 <- df700[!isin700,]



##mark bivariate outliers
par(mfrow=c(2,4))
plot(xx70, yy70, xlab=xlab, ylab=ylab, main="Pressure level: 70", type='n',xlim=c(-85,120), ylim=c(-100,100));points(xx70,yy70,pch=16,col=pointcolor,cex=0.4);polygon(cnt70$intpts,border=cntcolor);points(point_out70$xx70,point_out70$yy70,pch=16,col="red",cex=0.6)
plot(xx100, yy100, xlab=xlab, ylab=ylab, main="Pressure level: 100", type='n',xlim=c(-85,120), ylim=c(-100,100));points(xx100,yy100,pch=16,col=pointcolor,cex=0.4);polygon(cnt100$intpts,border=cntcolor);points(point_out100$xx100,point_out100$yy100,pch=16,col="red",cex=0.6)
plot(xx200, yy200, xlab=xlab, ylab=ylab, main="Pressure level: 200", type='n',xlim=c(-85,120), ylim=c(-100,100));points(xx200,yy200,pch=16,col=pointcolor,cex=0.4);polygon(cnt200$intpts,border=cntcolor);points(point_out200$xx200,point_out200$yy200,pch=16,col="red",cex=0.6)
plot(xx250, yy250, xlab=xlab, ylab=ylab, main="Pressure level: 250", type='n',xlim=c(-85,120), ylim=c(-100,100));points(xx250,yy250,pch=16,col=pointcolor,cex=0.4);polygon(cnt250$intpts,border=cntcolor);points(point_out250$xx250,point_out250$yy250,pch=16,col="red",cex=0.6)
plot(xx300, yy300, xlab=xlab, ylab=ylab, main="Pressure level: 300", type='n',xlim=c(-85,120), ylim=c(-100,100));points(xx300,yy300,pch=16,col=pointcolor,cex=0.4);polygon(cnt300$intpts,border=cntcolor);points(point_out300$xx300,point_out300$yy300,pch=16,col="red",cex=0.6)
plot(xx400, yy400, xlab=xlab, ylab=ylab, main="Pressure level: 400", type='n',xlim=c(-85,120), ylim=c(-100,100));points(xx400,yy400,pch=16,col=pointcolor,cex=0.4);polygon(cnt400$intpts,border=cntcolor);points(point_out400$xx400,point_out400$yy400,pch=16,col="red",cex=0.6)
plot(xx500, yy500, xlab=xlab, ylab=ylab, main="Pressure level: 500", type='n',xlim=c(-85,120), ylim=c(-100,100));points(xx500,yy500,pch=16,col=pointcolor,cex=0.4);polygon(cnt500$intpts,border=cntcolor);points(point_out500$xx500,point_out500$yy500,pch=16,col="red",cex=0.6)
plot(xx700, yy700, xlab=xlab, ylab=ylab, main="Pressure level: 700", type='n',xlim=c(-85,120), ylim=c(-100,100));points(xx700,yy700,pch=16,col=pointcolor,cex=0.4);polygon(cnt700$intpts,border=cntcolor);points(point_out700$xx700,point_out700$yy700,pch=16,col="red",cex=0.6)




######################## Figure 2 #######################################
#Outlier Detection and Visualization

time1=year[complete]<=1986
time2=year[complete]>=1987

horizontalspeeds_t1=c(xx70[time1],xx100[time1],xx200[time1],xx250[time1],xx300[time1],xx400[time1],xx500[time1],xx700[time1])
verticalspeeds_t1=c(yy70[time1],yy100[time1],yy200[time1],yy250[time1], yy300[time1],yy400[time1],yy500[time1],yy700[time1])                       
pressure_t1=c(rep(70,sum(time1)), rep(100,sum(time1)),
              rep(200,sum(time1)), rep(250,sum(time1)),
              rep(300,sum(time1)), rep(400,sum(time1)),
              rep(500,sum(time1)), rep(700,sum(time1)))


horizontalspeeds_t2=c(xx70[time2],xx100[time2],xx200[time2],xx250[time2],xx300[time2],xx400[time2],xx500[time2],xx700[time2])
verticalspeeds_t2=c(yy70[time2],yy100[time2],yy200[time2],yy250[time2], yy300[time2],yy400[time2],yy500[time2],yy700[time2])                       
pressure_t2=c(rep(70,sum(time2)), rep(100,sum(time2)),
              rep(200,sum(time2)), rep(250,sum(time2)),
              rep(300,sum(time2)), rep(400,sum(time2)),
              rep(500,sum(time2)), rep(700,sum(time2)))

 

#### Figure 2(a)########
# For first time period


xx=horizontalspeeds_t1 
yy=verticalspeeds_t1
ttt =pressure_t1  
#
xlab='Horizontal wind speed'
ylab='Vertical wind speed'


dirs=100
days=c(70,100,200,250,300,400,500,700)
clrs=c('black','black','black')
pchs=c(3,4,8)

library(quantreg)
library(splines)

cntpbs=c(0.05)    ### for pth quantiles

qd = qdir(dirs)
df=6
predmat = matrix(0,length(days),dirs)   
xy = cbind(xx,yy)


for (kk in 1:dirs)
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=cntpbs[1]),list(ttt=days))
##predicting the tau^th directional quantile



for (kk in 1:length(days))
{
  qln=cbind(qd,-predmat[kk,])
  clnid=1
  cln=1
  stop=FALSE
  filn=c(1:dirs)
  ss=0
  sln=0
  while (!stop)
  {
    abc=dertln(filn,clnid)
    a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
    b=qln[abc[2],]
    c=qln[abc[3],]
    tt=actln(a,b,c)           ## active line:1 else 0
    if (ss==0&tt==1) sln=filn[clnid]
    if (tt*ss==1&cln==sln) stop=TRUE
    #print('ss,tt,sln,cln')
    #print(c(ss,tt,sln,cln))
    ss=tt
    m=length(filn)
    if (tt==1)
    {
      if (clnid==m) {cln=filn[1] 
      clnid=1} else 
      {cln=filn[clnid+1] 
      clnid=clnid+1}
    } else { 
      clnidt=clnid
      if (clnid==1) {cln=filn[m] 
      clnid=m-1} else 
      {cln=filn[clnid-1] 
      clnid=clnid-1}
      filn=filn[-clnidt]
    }
  } 
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts=matrix(rep(0,2*m),nrow=m)
  intpts700_0.05=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(kk==1){intpts70_0.05t1=matrix(rep(0,2*m),nrow=m);intpts70_0.05t1=intpts} else if(kk==2){intpts100_0.05t1=matrix(rep(0,2*m),nrow=m);intpts100_0.05t1=intpts} else if(kk==3){intpts200_0.05t1=matrix(rep(0,2*m),nrow=m);intpts200_0.05t1=intpts} else if(kk==4){intpts250_0.05t1=matrix(rep(0,2*m),nrow=m);intpts250_0.05t1=intpts} else if(kk==5){intpts300_0.05t1=matrix(rep(0,2*m),nrow=m);intpts300_0.05t1=intpts} else if(kk==6){intpts400_0.05t1=matrix(rep(0,2*m),nrow=m);intpts400_0.05t1=intpts} else if(kk==7){intpts500_0.05t1=matrix(rep(0,2*m),nrow=m);intpts500_0.05t1=intpts} else{intpts700_0.05t1=matrix(rep(0,2*m),nrow=m);intpts700_0.05t1=intpts}
  
  #polygon(intpts,border=brewer.pal(8,"Spectral")[kk],lwd=2)   #intpts contains the cordinates the vertices of polygon   
}




################### p=0.0125  ########
cntpbs=c(0.0125)
for (kk in 1:dirs)
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=cntpbs[1]),list(ttt=days))
##predicting the tau^th directional quantile


for (kk in 1:length(days))
{
  qln=cbind(qd,-predmat[kk,])
  clnid=1
  cln=1
  stop=FALSE
  filn=c(1:dirs)
  ss=0
  sln=0
  while (!stop)
  {
    abc=dertln(filn,clnid)
    a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
    b=qln[abc[2],]
    c=qln[abc[3],]
    tt=actln(a,b,c)           ## active line:1 else 0
    if (ss==0&tt==1) sln=filn[clnid]
    if (tt*ss==1&cln==sln) stop=TRUE
    #print('ss,tt,sln,cln')
    #print(c(ss,tt,sln,cln))
    ss=tt
    m=length(filn)
    if (tt==1)
    {
      if (clnid==m) {cln=filn[1] 
      clnid=1} else 
      {cln=filn[clnid+1] 
      clnid=clnid+1}
    } else { 
      clnidt=clnid
      if (clnid==1) {cln=filn[m] 
      clnid=m-1} else 
      {cln=filn[clnid-1] 
      clnid=clnid-1}
      filn=filn[-clnidt]
    }
  } 
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(kk==1){intpts70_0.0125t1=matrix(rep(0,2*m),nrow=m);intpts70_0.0125t1=intpts} else if(kk==2){intpts100_0.0125t1=matrix(rep(0,2*m),nrow=m);intpts100_0.0125t1=intpts} else if(kk==3){intpts200_0.0125t1=matrix(rep(0,2*m),nrow=m);intpts200_0.0125t1=intpts} else if(kk==4){intpts250_0.0125t1=matrix(rep(0,2*m),nrow=m);intpts250_0.0125t1=intpts} else if(kk==5){intpts300_0.0125t1=matrix(rep(0,2*m),nrow=m);intpts300_0.0125t1=intpts} else if(kk==6){intpts400_0.0125t1=matrix(rep(0,2*m),nrow=m);intpts400_0.0125t1=intpts} else if(kk==7){intpts500_0.0125t1=matrix(rep(0,2*m),nrow=m);intpts500_0.0125t1=intpts} else{intpts700_0.0125t1=matrix(rep(0,2*m),nrow=m);intpts700_0.0125t1=intpts}
  
}



### median
median_launch_t1 =matrix(0,length(days), 2)
## obtaining median launch
probs = rep(0.495, length(days))
for (kk in 1:length(days))
{
  stop=FALSE
  while(!stop){
    for (s in 1:dirs)
      predmat[,s] = predict(rq(xy %*% qd[s,] ~ bs(ttt,df=df), tau=probs[kk]),list(ttt=days))
    ##predicting the tau^th directional quantile
    
    
    qln=cbind(qd,-predmat[kk,])
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
      b=qln[abc[2],]
      c=qln[abc[3],]
      if(abs(a[1])==abs(c[1]) & abs(a[2])==abs(c[2])) {
        probs[kk]=probs[kk]-0.005
        break}
      tt=actln(a,b,c)           ## active line:1 else 0
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
    
  }   
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts=matrix(rep(0,2*m),nrow=m)
  
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts) > 1e+10 )>=2){
    intpts=intpts[-which(abs(intpts[,1])>1e+10),]}
  
  #if(kk==1){median_sim70=matrix(rep(0,2*m),nrow=m); median_sim70=intpts} else if(kk==2){median_sim100=matrix(rep(0,2*m),nrow=m);median_sim100=intpts} else if(kk==3){median_sim200=matrix(rep(0,2*m),nrow=m);median_sim200=intpts} else if(kk==4){median_sim250=matrix(rep(0,2*m),nrow=m);median_sim250=intpts} else if(kk==5){median_sim300=matrix(rep(0,2*m),nrow=m);median_sim300=intpts} else if(kk==6){median_sim400=matrix(rep(0,2*m),nrow=m);median_sim400=intpts} else if(kk==7){median_sim500=matrix(rep(0,2*m),nrow=m);median_sim500=intpts} else{median_sim700=matrix(rep(0,2*m),nrow=m);median_sim700=intpts}
  
  median_launch_t1[kk,] = c(mean(intpts[,1]),mean(intpts[,2]))
}

####### outliers for first time period
xx70t1=xx70[time1]; xx100t1=xx100[time1]; xx200t1=xx200[time1]; xx250t1=xx250[time1]; xx300t1=xx300[time1]; xx400t1=xx400[time1]; xx500t1=xx500[time1];xx700t1=xx700[time1]
yy70t1=yy70[time1]; yy100t1=yy100[time1]; yy200t1=yy200[time1]; yy250t1=yy250[time1]; yy300t1=yy300[time1]; yy400t1=yy400[time1]; yy500t1=yy500[time1];yy700t1=yy700[time1]

yearcomp=year[complete];juliancomp=julian[complete];hourcomp=hour[complete]

library(evir)

dirs=1000
qd = qdir(dirs)      ## (costheta, sintheta) nrows=dirs
cntpbs=c(0.0002,0.0001,0.00001,0.000001)

qln70t1 = cbind(qd,-apply(qd,1,
                          function(u)
                            -riskmeasures(gpd(cbind(xx70t1,yy70t1) %*% -u, nextremes=263),1-cntpbs[2])[2] ))

qln100t1 = cbind(qd,-apply(qd,1,
                           function(u)
                             -riskmeasures(gpd(cbind(xx100t1,yy100t1) %*% -u, nextremes=263),1-cntpbs[2])[2] ))
qln200t1 = cbind(qd,-apply(qd,1,
                           function(u)
                             -riskmeasures(gpd(cbind(xx200t1,yy200t1) %*% -u, nextremes=263),1-cntpbs[2])[2] ))
qln250t1 = cbind(qd,-apply(qd,1,
                           function(u)
                             -riskmeasures(gpd(cbind(xx250t1,yy250t1) %*% -u, nextremes=263),1-cntpbs[2])[2] ))
qln300t1 = cbind(qd,-apply(qd,1,
                           function(u)
                             -riskmeasures(gpd(cbind(xx300t1,yy300t1) %*% -u, nextremes=263),1-cntpbs[2])[2] ))
qln400t1 = cbind(qd,-apply(qd,1,
                           function(u)
                             -riskmeasures(gpd(cbind(xx400t1,yy400t1) %*% -u, nextremes=263),1-cntpbs[2])[2] ))
qln500t1 = cbind(qd,-apply(qd,1,
                           function(u)
                             -riskmeasures(gpd(cbind(xx500t1,yy500t1) %*% -u, nextremes=263),1-cntpbs[2])[2] ))
qln700t1 = cbind(qd,-apply(qd,1,
                           function(u)
                             -riskmeasures(gpd(cbind(xx700t1,yy700t1) %*% -u, nextremes=263),1-cntpbs[2])[2] ))

cnt70t1 = eli(qln70t1); cnt100t1 = eli(qln100t1); cnt200t1 = eli(qln200t1); cnt250t1 = eli(qln250t1); cnt300t1 = eli(qln300t1); cnt400t1 = eli(qln400t1); cnt500t1 = eli(qln500t1); cnt700t1 = eli(qln700t1)

library(spatstat)
df70t1=data.frame(xx70t1,yy70t1,yearcomp[time1],juliancomp[time1],hourcomp[time1]); df100t1=data.frame(xx100t1,yy100t1,yearcomp[time1],juliancomp[time1],hourcomp[time1]);df200t1=data.frame(xx200t1,yy200t1,yearcomp[time1],juliancomp[time1],hourcomp[time1]);df250t1=data.frame(xx250t1,yy250t1,yearcomp[time1],juliancomp[time1],hourcomp[time1]);df300t1=data.frame(xx300t1,yy300t1,yearcomp[time1],juliancomp[time1],hourcomp[time1]);df400t1=data.frame(xx400t1,yy400t1,yearcomp[time1],juliancomp[time1],hourcomp[time1]); df500t1=data.frame(xx500t1,yy500t1,yearcomp[time1],juliancomp[time1],hourcomp[time1]);df700t1=data.frame(xx700t1,yy700t1,yearcomp[time1],juliancomp[time1],hourcomp[time1])
bound70t1 <- owin(poly=data.frame(x=cnt70t1$intpts[,1], y=cnt70t1$intpts[,2])); bound100t1 <- owin(poly=data.frame(x=cnt100t1$intpts[,1], y=cnt100t1$intpts[,2])); bound200t1 <- owin(poly=data.frame(x=cnt200t1$intpts[,1], y=cnt200t1$intpts[,2])); bound250t1 <- owin(poly=data.frame(x=cnt250t1$intpts[,1], y=cnt250t1$intpts[,2])); bound300t1 <- owin(poly=data.frame(x=cnt300t1$intpts[,1], y=cnt300t1$intpts[,2])); bound400t1 <- owin(poly=data.frame(x=cnt400t1$intpts[,1], y=cnt400t1$intpts[,2])); bound500t1 <- owin(poly=data.frame(x=cnt500t1$intpts[,1], y=cnt500t1$intpts[,2])); bound700t1 <- owin(poly=data.frame(x=cnt700t1$intpts[,1], y=cnt700t1$intpts[,2]))
isin70t1<-inside.owin(x=xx70t1,y=yy70t1,w=bound70t1); isin100t1<-inside.owin(x=xx100t1,y=yy100t1,w=bound100t1); isin200t1<-inside.owin(x=xx200t1,y=yy200t1,w=bound200t1); isin250t1<-inside.owin(x=xx250t1,y=yy250t1,w=bound250t1); isin300t1<-inside.owin(x=xx300t1,y=yy300t1,w=bound300t1); isin400t1<-inside.owin(x=xx400t1,y=yy400t1,w=bound400t1); isin500t1<-inside.owin(x=xx500t1,y=yy500t1,w=bound500t1); isin700t1<-inside.owin(x=xx700t1,y=yy700t1,w=bound700t1)
point_out70t1 <- df70t1[!isin70t1,]; point_out100t1 <- df100t1[!isin100t1,]; point_out200t1 <- df200t1[!isin200t1,]; point_out250t1 <- df250t1[!isin250t1,]; point_out300t1 <- df300t1[!isin300t1,]; point_out400t1 <- df400t1[!isin400t1,]; point_out500t1 <- df500t1[!isin500t1,]; point_out700t1 <- df700t1[!isin700t1,]

names(point_out100t1)=names(point_out70t1); names(point_out200t1)=names(point_out70t1); names(point_out250t1)=names(point_out70t1); names(point_out300t1)=names(point_out70t1); names(point_out400t1)=names(point_out70t1); names(point_out500t1)=names(point_out70t1); names(point_out700t1)=names(point_out70t1)
outlierst1=unique(data.frame(rbind(point_out70t1[-c(1,2)],point_out100t1[-c(1,2)],point_out200t1[-c(1,2)],point_out250t1[-c(1,2)],point_out300t1[-c(1,2)],point_out400t1[-c(1,2)],point_out500t1[-c(1,2)],point_out700t1[-c(1,2)])))
indext1=as.numeric(rownames(outlierst1))

library(scatterplot3d)
library(plot3D)
library(RColorBrewer)
horizontalspeedt1=rbind(xx70t1[indext1],xx100t1[indext1],xx200t1[indext1],xx250t1[indext1],xx300t1[indext1],xx400t1[indext1],xx500t1[indext1],xx700t1[indext1])
verticalspeedt1=rbind(yy70t1[indext1],yy100t1[indext1],yy200t1[indext1],yy250t1[indext1],yy300t1[indext1],yy400t1[indext1],yy500t1[indext1],yy700t1[indext1])
pressurelevelt1=c(70,100,200,250,300,400,500,700)
col=brewer.pal(8,"Blues")

s=scatterplot3d(horizontalspeedt1[,1],verticalspeedt1[,1],pressurelevelt1,color="red",pch=20, 
                xlim = c(-90,140), ylim = c(-112,88),zlim = c(0,800),cex.symbols = 0.5,
                xlab = "Horizontal speed", ylab = "",zlab = "Pressure level",
                main="1962-1986",angle=10)

dims <- par("usr")
x <- dims[1]+ 0.9*diff(dims[1:2])-1.3
y <- dims[3]+ 0.08*diff(dims[3:4])-0.22
text(x,y,"Vertical speed",srt=10)

for(i in 2:length(indext1)){
  s$points3d(horizontalspeedt1[,i],verticalspeedt1[,i],pressurelevelt1,col="red",pch=20,cex=0.5)
}

### Extreme quantile envelope p=0.0001 ####
s$points3d(cnt70t1$intpts[,1],cnt70t1$intpts[,2],rep(70,nrow(cnt70t1$intpts)),type="l",lwd=2,lty=2,col=col[8])
s$points3d(cnt100t1$intpts[,1],cnt100t1$intpts[,2],rep(100,nrow(cnt100t1$intpts)),type="l",lwd=2,lty=2,col=col[7])
s$points3d(cnt200t1$intpts[,1],cnt200t1$intpts[,2],rep(200,nrow(cnt200t1$intpts)),type="l",lwd=2,lty=2,col=col[6])
s$points3d(cnt250t1$intpts[,1],cnt250t1$intpts[,2],rep(250,nrow(cnt250t1$intpts)),type="l",lwd=2,lty=2,col=col[5])
s$points3d(cnt300t1$intpts[,1],cnt300t1$intpts[,2],rep(300,nrow(cnt300t1$intpts)),type="l",lwd=2,lty=2,col=col[4])
s$points3d(cnt400t1$intpts[,1],cnt400t1$intpts[,2],rep(400,nrow(cnt400t1$intpts)),type="l",lwd=2,lty=2,col=col[3])
s$points3d(cnt500t1$intpts[,1],cnt500t1$intpts[,2],rep(500,nrow(cnt500t1$intpts)),type="l",lwd=2,lty=2,col=col[2])
s$points3d(cnt700t1$intpts[,1],cnt700t1$intpts[,2],rep(700,nrow(cnt700t1$intpts)),type="l",lwd=2,lty=2,col=col[1])

### Predicted quantile envelope using complete data for p=0.05
s$points3d(intpts70_0.05t1[,1],intpts70_0.05t1[,2],rep(70,nrow(intpts70_0.05t1)),type="l",lwd=2,col=col[8])
s$points3d(intpts100_0.05t1[,1],intpts100_0.05t1[,2],rep(100,nrow(intpts100_0.05t1)),type="l",lwd=2,col=col[7])
s$points3d(intpts200_0.05t1[,1],intpts200_0.05t1[,2],rep(200,nrow(intpts200_0.05t1)),type="l",lwd=2,col=col[6])
s$points3d(intpts250_0.05t1[,1],intpts250_0.05t1[,2],rep(250,nrow(intpts250_0.05t1)),type="l",lwd=2,col=col[5])
s$points3d(intpts300_0.05t1[,1],intpts300_0.05t1[,2],rep(300,nrow(intpts300_0.05t1)),type="l",lwd=2,col=col[4])
s$points3d(intpts400_0.05t1[,1],intpts400_0.05t1[,2],rep(400,nrow(intpts400_0.05t1)),type="l",lwd=2,col=col[3])
s$points3d(intpts500_0.05t1[,1],intpts500_0.05t1[,2],rep(500,nrow(intpts500_0.05t1)),type="l",lwd=2,col=col[2])
s$points3d(intpts700_0.05t1[,1],intpts700_0.05t1[,2],rep(700,nrow(intpts700_0.05t1)),type="l",lwd=2,col=col[1])

### Predicted quantile envelope using complete data for p=0.0125
s$points3d(intpts70_0.0125t1[,1],intpts70_0.0125t1[,2],rep(70,nrow(intpts70_0.0125t1)),type="l",lwd=2,col=col[8])
s$points3d(intpts100_0.0125t1[,1],intpts100_0.0125t1[,2],rep(100,nrow(intpts100_0.0125t1)),type="l",lwd=2,col=col[7])
s$points3d(intpts200_0.0125t1[,1],intpts200_0.0125t1[,2],rep(200,nrow(intpts200_0.0125t1)),type="l",lwd=2,col=col[6])
s$points3d(intpts250_0.0125t1[,1],intpts250_0.0125t1[,2],rep(250,nrow(intpts250_0.0125t1)),type="l",lwd=2,col=col[5])
s$points3d(intpts300_0.0125t1[,1],intpts300_0.0125t1[,2],rep(300,nrow(intpts300_0.0125t1)),type="l",lwd=2,col=col[4])
s$points3d(intpts400_0.0125t1[,1],intpts400_0.0125t1[,2],rep(400,nrow(intpts400_0.0125t1)),type="l",lwd=2,col=col[3])
s$points3d(intpts500_0.0125t1[,1],intpts500_0.0125t1[,2],rep(500,nrow(intpts500_0.0125t1)),type="l",lwd=2,col=col[2])
s$points3d(intpts700_0.0125t1[,1],intpts700_0.0125t1[,2],rep(700,nrow(intpts700_0.0125t1)),type="l",lwd=2,col=col[1])

legend("topright", c("70","100","200","250","300","400","500","700"), 
       lty=c(1,1), lwd=c(2.5,2.5),col=rev(brewer.pal(8,"Blues")),cex=0.6) 

s$points3d(median_launch_t1[,1],median_launch_t1[,2],pressurelevelt1,type="l",col='black', lwd=2)


##### Figure 2(b) #######
## For second time period

xx=horizontalspeeds_t2 
yy=verticalspeeds_t2
ttt =pressure_t2 
#
dirs=100

days=c(70,100,200,250,300,400,500,700)
clrs=c('black','black','black')
pchs=c(3,4,8)

qd = qdir(dirs)
df=6
predmat = matrix(0,length(days),dirs)   
dim(predmat)
xy = cbind(xx,yy)


##### 0.05 ################

cntpbs=c(0.05)    ### for pth quantiles


for (kk in 1:dirs)
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=cntpbs[1]),list(ttt=days))
##predicting the tau^th directional quantile


for (kk in 1:length(days))
{
  qln=cbind(qd,-predmat[kk,])
  clnid=1
  cln=1
  stop=FALSE
  filn=c(1:dirs)
  ss=0
  sln=0
  while (!stop)
  {
    abc=dertln(filn,clnid)
    a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
    b=qln[abc[2],]
    c=qln[abc[3],]
    tt=actln(a,b,c)           ## active line:1 else 0
    if (ss==0&tt==1) sln=filn[clnid]
    if (tt*ss==1&cln==sln) stop=TRUE
    #print('ss,tt,sln,cln')
    #print(c(ss,tt,sln,cln))
    ss=tt
    m=length(filn)
    if (tt==1)
    {
      if (clnid==m) {cln=filn[1] 
      clnid=1} else 
      {cln=filn[clnid+1] 
      clnid=clnid+1}
    } else { 
      clnidt=clnid
      if (clnid==1) {cln=filn[m] 
      clnid=m-1} else 
      {cln=filn[clnid-1] 
      clnid=clnid-1}
      filn=filn[-clnidt]
    }
  } 
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts=matrix(rep(0,2*m),nrow=m)
  intpts700_0.05=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(kk==1){intpts70_0.05t2=matrix(rep(0,2*m),nrow=m);intpts70_0.05t2=intpts} else if(kk==2){intpts100_0.05t2=matrix(rep(0,2*m),nrow=m);intpts100_0.05t2=intpts} else if(kk==3){intpts200_0.05t2=matrix(rep(0,2*m),nrow=m);intpts200_0.05t2=intpts} else if(kk==4){intpts250_0.05t2=matrix(rep(0,2*m),nrow=m);intpts250_0.05t2=intpts} else if(kk==5){intpts300_0.05t2=matrix(rep(0,2*m),nrow=m);intpts300_0.05t2=intpts} else if(kk==6){intpts400_0.05t2=matrix(rep(0,2*m),nrow=m);intpts400_0.05t2=intpts} else if(kk==7){intpts500_0.05t2=matrix(rep(0,2*m),nrow=m);intpts500_0.05t2=intpts} else{intpts700_0.05t2=matrix(rep(0,2*m),nrow=m);intpts700_0.05t2=intpts}
  
  #polygon(intpts,border=brewer.pal(8,"Spectral")[kk],lwd=2)   #intpts contains the cordinates the vertices of polygon   
}


##### 0.0125 ################

cntpbs=c(0.0125)    ### for pth quantiles

#for (k in 1:length(cntpbs))
#{
for (kk in 1:dirs)
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=cntpbs[1]),list(ttt=days))
##predicting the tau^th directional quantile


for (kk in 1:length(days))
{
  qln=cbind(qd,-predmat[kk,])
  clnid=1
  cln=1
  stop=FALSE
  filn=c(1:dirs)
  ss=0
  sln=0
  while (!stop)
  {
    abc=dertln(filn,clnid)
    a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
    b=qln[abc[2],]
    c=qln[abc[3],]
    tt=actln(a,b,c)           ## active line:1 else 0
    if (ss==0&tt==1) sln=filn[clnid]
    if (tt*ss==1&cln==sln) stop=TRUE
    #print('ss,tt,sln,cln')
    #print(c(ss,tt,sln,cln))
    ss=tt
    m=length(filn)
    if (tt==1)
    {
      if (clnid==m) {cln=filn[1] 
      clnid=1} else 
      {cln=filn[clnid+1] 
      clnid=clnid+1}
    } else { 
      clnidt=clnid
      if (clnid==1) {cln=filn[m] 
      clnid=m-1} else 
      {cln=filn[clnid-1] 
      clnid=clnid-1}
      filn=filn[-clnidt]
    }
  } 
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts=matrix(rep(0,2*m),nrow=m)
  intpts700_0.05=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(kk==1){intpts70_0.0125t2=matrix(rep(0,2*m),nrow=m);intpts70_0.0125t2=intpts} else if(kk==2){intpts100_0.0125t2=matrix(rep(0,2*m),nrow=m);intpts100_0.0125t2=intpts} else if(kk==3){intpts200_0.0125t2=matrix(rep(0,2*m),nrow=m);intpts200_0.0125t2=intpts} else if(kk==4){intpts250_0.0125t2=matrix(rep(0,2*m),nrow=m);intpts250_0.0125t2=intpts} else if(kk==5){intpts300_0.0125t2=matrix(rep(0,2*m),nrow=m);intpts300_0.0125t2=intpts} else if(kk==6){intpts400_0.0125t2=matrix(rep(0,2*m),nrow=m);intpts400_0.0125t2=intpts} else if(kk==7){intpts500_0.0125t2=matrix(rep(0,2*m),nrow=m);intpts500_0.0125t2=intpts} else{intpts700_0.0125t2=matrix(rep(0,2*m),nrow=m);intpts700_0.0125t2=intpts}
  
  #polygon(intpts,border=brewer.pal(8,"Spectral")[kk],lwd=2)   #intpts contains the cordinates the vertices of polygon   
}


### median
median_launch_t2 =matrix(0,length(days), 2)
## obtaining median launch
probs = rep(0.495, length(days))
for (kk in 1:length(days))
{
  stop=FALSE
  while(!stop){
    for (s in 1:dirs)
      predmat[,s] = predict(rq(xy %*% qd[s,] ~ bs(ttt,df=df), tau=probs[kk]),list(ttt=days))
    ##predicting the tau^th directional quantile
    
    
    qln=cbind(qd,-predmat[kk,])
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
      b=qln[abc[2],]
      c=qln[abc[3],]
      if(abs(a[1])==abs(c[1]) & abs(a[2])==abs(c[2])) {
        probs[kk]=probs[kk]-0.005
        break}
      tt=actln(a,b,c)           ## active line:1 else 0
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
    
  }   
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts=matrix(rep(0,2*m),nrow=m)
  
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts) > 1e+10 )>=2){
    intpts=intpts[-which(abs(intpts[,1])>1e+10),]}
  
  median_launch_t2[kk,] = c(mean(intpts[,1]),mean(intpts[,2]))
}

####### outliers for second time period
xx70t2=xx70[time2]; xx100t2=xx100[time2]; xx200t2=xx200[time2]; xx250t2=xx250[time2]; xx300t2=xx300[time2]; xx400t2=xx400[time2]; xx500t2=xx500[time2];xx700t2=xx700[time2]
yy70t2=yy70[time2]; yy100t2=yy100[time2]; yy200t2=yy200[time2]; yy250t2=yy250[time2]; yy300t2=yy300[time2]; yy400t2=yy400[time2]; yy500t2=yy500[time2];yy700t2=yy700[time2]

yearcomp=year[complete];juliancomp=julian[complete];hourcomp=hour[complete]


dirs=1000
qd = qdir(dirs)      ## (costheta, sintheta) nrows=dirs


qln70t2 = cbind(qd,-apply(qd,1,
                          function(u)
                            -riskmeasures(gpd(cbind(xx70t2,yy70t2) %*% -u, nextremes=310),1-cntpbs[2])[2] ))

qln100t2 = cbind(qd,-apply(qd,1,
                           function(u)
                             -riskmeasures(gpd(cbind(xx100t2,yy100t2) %*% -u, nextremes=310),1-cntpbs[2])[2] ))
qln200t2 = cbind(qd,-apply(qd,1,
                           function(u)
                             -riskmeasures(gpd(cbind(xx200t2,yy200t2) %*% -u, nextremes=310),1-cntpbs[2])[2] ))
qln250t2 = cbind(qd,-apply(qd,1,
                           function(u)
                             -riskmeasures(gpd(cbind(xx250t2,yy250t2) %*% -u, nextremes=310),1-cntpbs[2])[2] ))
qln300t2 = cbind(qd,-apply(qd,1,
                           function(u)
                             -riskmeasures(gpd(cbind(xx300t2,yy300t2) %*% -u, nextremes=310),1-cntpbs[2])[2] ))
qln400t2 = cbind(qd,-apply(qd,1,
                           function(u)
                             -riskmeasures(gpd(cbind(xx400t2,yy400t2) %*% -u, nextremes=310),1-cntpbs[2])[2] ))
qln500t2 = cbind(qd,-apply(qd,1,
                           function(u)
                             -riskmeasures(gpd(cbind(xx500t2,yy500t2) %*% -u, nextremes=310),1-cntpbs[2])[2] ))
qln700t2 = cbind(qd,-apply(qd,1,
                           function(u)
                             -riskmeasures(gpd(cbind(xx700t2,yy700t2) %*% -u, nextremes=310),1-cntpbs[2])[2] ))

cnt70t2 = eli(qln70t2); cnt100t2 = eli(qln100t2); cnt200t2 = eli(qln200t2); cnt250t2 = eli(qln250t2); cnt300t2 = eli(qln300t2); cnt400t2 = eli(qln400t2); cnt500t2 = eli(qln500t2); cnt700t2 = eli(qln700t2)


library(spatstat)
df70t2=data.frame(xx70t2,yy70t2,yearcomp[time2],juliancomp[time2],hourcomp[time2]); df100t2=data.frame(xx100t2,yy100t2,yearcomp[time2],juliancomp[time2],hourcomp[time2]);df200t2=data.frame(xx200t2,yy200t2,yearcomp[time2],juliancomp[time2],hourcomp[time2]);df250t2=data.frame(xx250t2,yy250t2,yearcomp[time2],juliancomp[time2],hourcomp[time2]);df300t2=data.frame(xx300t2,yy300t2,yearcomp[time2],juliancomp[time2],hourcomp[time2]);df400t2=data.frame(xx400t2,yy400t2,yearcomp[time2],juliancomp[time2],hourcomp[time2]); df500t2=data.frame(xx500t2,yy500t2,yearcomp[time2],juliancomp[time2],hourcomp[time2]);df700t2=data.frame(xx700t2,yy700t2,yearcomp[time2],juliancomp[time2],hourcomp[time2])
bound70t2 <- owin(poly=data.frame(x=cnt70t2$intpts[,1], y=cnt70t2$intpts[,2])); bound100t2 <- owin(poly=data.frame(x=cnt100t2$intpts[,1], y=cnt100t2$intpts[,2])); bound200t2 <- owin(poly=data.frame(x=cnt200t2$intpts[,1], y=cnt200t2$intpts[,2])); bound250t2 <- owin(poly=data.frame(x=cnt250t2$intpts[,1], y=cnt250t2$intpts[,2])); bound300t2 <- owin(poly=data.frame(x=cnt300t2$intpts[,1], y=cnt300t2$intpts[,2])); bound400t2 <- owin(poly=data.frame(x=cnt400t2$intpts[,1], y=cnt400t2$intpts[,2])); bound500t2 <- owin(poly=data.frame(x=cnt500t2$intpts[,1], y=cnt500t2$intpts[,2])); bound700t2 <- owin(poly=data.frame(x=cnt700t2$intpts[,1], y=cnt700t2$intpts[,2]))
isin70t2<-inside.owin(x=xx70t2,y=yy70t2,w=bound70t2); isin100t2<-inside.owin(x=xx100t2,y=yy100t2,w=bound100t2); isin200t2<-inside.owin(x=xx200t2,y=yy200t2,w=bound200t2); isin250t2<-inside.owin(x=xx250t2,y=yy250t2,w=bound250t2); isin300t2<-inside.owin(x=xx300t2,y=yy300t2,w=bound300t2); isin400t2<-inside.owin(x=xx400t2,y=yy400t2,w=bound400t2); isin500t2<-inside.owin(x=xx500t2,y=yy500t2,w=bound500t2); isin700t2<-inside.owin(x=xx700t2,y=yy700t2,w=bound700t2)
point_out70t2 <- df70t2[!isin70t2,]; point_out100t2 <- df100t2[!isin100t2,]; point_out200t2 <- df200t2[!isin200t2,]; point_out250t2 <- df250t2[!isin250t2,]; point_out300t2 <- df300t2[!isin300t2,]; point_out400t2 <- df400t2[!isin400t2,]; point_out500t2 <- df500t2[!isin500t2,]; point_out700t2 <- df700t2[!isin700t2,]


names(point_out200t2)=names(point_out70t2); names(point_out200t2)=names(point_out70t2); names(point_out250t2)=names(point_out70t2); names(point_out300t2)=names(point_out70t2); names(point_out400t2)=names(point_out70t2); names(point_out500t2)=names(point_out70t2); names(point_out700t2)=names(point_out70t2)
outlierst2=unique(data.frame(rbind(point_out70t2[-c(1,2)],point_out200t2[-c(1,2)],point_out200t2[-c(1,2)],point_out250t2[-c(1,2)],point_out300t2[-c(1,2)],point_out400t2[-c(1,2)],point_out500t2[-c(1,2)],point_out700t2[-c(1,2)])))

indext2=as.numeric(rownames(outlierst2))

horizontalspeedt2=rbind(xx70t2[indext2],xx100t2[indext2],xx200t2[indext2],xx250t2[indext2],xx300t2[indext2],xx400t2[indext2],xx500t2[indext2],xx700t2[indext2])
verticalspeedt2=rbind(yy70t2[indext2],yy100t2[indext2],yy200t2[indext2],yy250t2[indext2],yy300t2[indext2],yy400t2[indext2],yy500t2[indext2],yy700t2[indext2])
pressurelevelt2=c(70,100,200,250,300,400,500,700)
col=brewer.pal(8,"Blues")
s=scatterplot3d(horizontalspeedt2[,1],verticalspeedt2[,1],pressurelevelt2,color="red",pch=20, 
                xlim = c(-90,140), ylim = c(-112,88),zlim = c(0,800), cex.symbols = 0.5,
                xlab = "Horizontal speed", ylab = "",zlab = "Pressure level",
                main="1987-2011",angle=10)
dims <- par("usr")
x <- dims[1]+ 0.9*diff(dims[1:2])-1.3
y <- dims[3]+ 0.08*diff(dims[3:4])-0.19
text(x,y,"Vertical speed",srt=10)


for(i in 2:length(indext2)){
  s$points3d(horizontalspeedt2[,i],verticalspeedt2[,i],pressurelevelt2,col="red",pch=20,cex = 0.5)
}

### Extreme quantile envelope p=0.0001 ####
s$points3d(cnt70t2$intpts[,1],cnt70t2$intpts[,2],rep(70,nrow(cnt70t2$intpts)),type="l",lwd=2,lty=2,col=col[8])
s$points3d(cnt100t2$intpts[,1],cnt100t2$intpts[,2],rep(100,nrow(cnt100t2$intpts)),type="l",lwd=2,lty=2,col=col[7])
s$points3d(cnt200t2$intpts[,1],cnt200t2$intpts[,2],rep(200,nrow(cnt200t2$intpts)),type="l",lwd=2,lty=2,col=col[6])
s$points3d(cnt250t2$intpts[,1],cnt250t2$intpts[,2],rep(250,nrow(cnt250t2$intpts)),type="l",lwd=2,lty=2,col=col[5])
s$points3d(cnt300t2$intpts[,1],cnt300t2$intpts[,2],rep(300,nrow(cnt300t2$intpts)),type="l",lwd=2,lty=2,col=col[4])
s$points3d(cnt400t2$intpts[,1],cnt400t2$intpts[,2],rep(400,nrow(cnt400t2$intpts)),type="l",lwd=2,lty=2,col=col[3])
s$points3d(cnt500t2$intpts[,1],cnt500t2$intpts[,2],rep(500,nrow(cnt500t2$intpts)),type="l",lwd=2,lty=2,col=col[2])
s$points3d(cnt700t2$intpts[,1],cnt700t2$intpts[,2],rep(700,nrow(cnt700t2$intpts)),type="l",lwd=2,lty=2,col=col[1])

### Predicted quantile envelope using complete data for p=0.05
s$points3d(intpts70_0.05t2[,1],intpts70_0.05t2[,2],rep(70,nrow(intpts70_0.05t2)),type="l",lwd=2,col=col[8])
s$points3d(intpts100_0.05t2[,1],intpts100_0.05t2[,2],rep(100,nrow(intpts100_0.05t2)),type="l",lwd=2,col=col[7])
s$points3d(intpts200_0.05t2[,1],intpts200_0.05t2[,2],rep(200,nrow(intpts200_0.05t2)),type="l",lwd=2,col=col[6])
s$points3d(intpts250_0.05t2[,1],intpts250_0.05t2[,2],rep(250,nrow(intpts250_0.05t2)),type="l",lwd=2,col=col[5])
s$points3d(intpts300_0.05t2[,1],intpts300_0.05t2[,2],rep(300,nrow(intpts300_0.05t2)),type="l",lwd=2,col=col[4])
s$points3d(intpts400_0.05t2[,1],intpts400_0.05t2[,2],rep(400,nrow(intpts400_0.05t2)),type="l",lwd=2,col=col[3])
s$points3d(intpts500_0.05t2[,1],intpts500_0.05t2[,2],rep(500,nrow(intpts500_0.05t2)),type="l",lwd=2,col=col[2])
s$points3d(intpts700_0.05t2[,1],intpts700_0.05t2[,2],rep(700,nrow(intpts700_0.05t2)),type="l",lwd=2,col=col[1])

### Predicted quantile envelope using complete data for p=0.0125
s$points3d(intpts70_0.0125t2[,1],intpts70_0.0125t2[,2],rep(70,nrow(intpts70_0.0125t2)),type="l",lwd=2,col=col[8])
s$points3d(intpts100_0.0125t2[,1],intpts100_0.0125t2[,2],rep(100,nrow(intpts100_0.0125t2)),type="l",lwd=2,col=col[7])
s$points3d(intpts200_0.0125t2[,1],intpts200_0.0125t2[,2],rep(200,nrow(intpts200_0.0125t2)),type="l",lwd=2,col=col[6])
s$points3d(intpts250_0.0125t2[,1],intpts250_0.0125t2[,2],rep(250,nrow(intpts250_0.0125t2)),type="l",lwd=2,col=col[5])
s$points3d(intpts300_0.0125t2[,1],intpts300_0.0125t2[,2],rep(300,nrow(intpts300_0.0125t2)),type="l",lwd=2,col=col[4])
s$points3d(intpts400_0.0125t2[,1],intpts400_0.0125t2[,2],rep(400,nrow(intpts400_0.0125t2)),type="l",lwd=2,col=col[3])
s$points3d(intpts500_0.0125t2[,1],intpts500_0.0125t2[,2],rep(500,nrow(intpts500_0.0125t2)),type="l",lwd=2,col=col[2])
s$points3d(intpts700_0.0125t2[,1],intpts700_0.0125t2[,2],rep(700,nrow(intpts700_0.0125t2)),type="l",lwd=2,col=col[1])

legend("topright", c("70","100","200","250","300","400","500","700"), 
       lty=c(1,1), lwd=c(2.5,2.5),col=rev(brewer.pal(8,"Blues")),cex=0.6) 


s$points3d(median_launch_t2[,1],median_launch_t2[,2],pressurelevelt2,type="l",col='black', lwd=2)





######################## Figure 3 #######################################
# connected envelopes

complete=complete.cases(u70,v70,u100,v100,u200,v200,u250,v250,u300,v300,u400,v400,u500,v500,u700,v700)

horz70=u70; horz100=u100; horz200=u200; horz250=u250; horz300=u300;horz400=u400; horz500=u500; horz700=u700; horz850=u850; horz1000=u1000
vert70=v70; vert100=v100; vert200=v200; vert250=v250; vert300=v300;vert400=v400; vert500=v500; vert700=v700; vert850=v850; vert1000=v1000 

xx70=horz70[complete]; xx100=horz100[complete]; xx200=horz200[complete]; xx250=horz250[complete]; xx300=horz300[complete]; xx400=horz400[complete]; xx500=horz500[complete];xx700=horz700[complete]
yy70=vert70[complete]; yy100=vert100[complete]; yy200=vert200[complete]; yy250=vert250[complete]; yy300=vert300[complete]; yy400=vert400[complete]; yy500=vert500[complete];yy700=vert700[complete]
t= c(70,100,200,250,300,400,500,700)

time1=year[complete]<=1986  ## 5259 launches
time2=year[complete]>=1987  ## 6193 launches


library(DepthProc)
comb70 = data.frame(xx70[time1],yy70[time1]); comb100 = data.frame(xx100[time1],yy100[time1]); comb200 = data.frame(xx200[time1],yy200[time1]); comb250 = data.frame(xx250[time1],yy250[time1]); comb300 = data.frame(xx300[time1],yy300[time1]); comb400 = data.frame(xx400[time1],yy400[time1]); comb500 = data.frame(xx500[time1],yy500[time1]); comb700 = data.frame(xx700[time1],yy700[time1])
comb70 = data.frame(xx70[time2],yy70[time2]); comb100 = data.frame(xx100[time2],yy100[time2]); comb200 = data.frame(xx200[time2],yy200[time2]); comb250 = data.frame(xx250[time2],yy250[time2]); comb300 = data.frame(xx300[time2],yy300[time2]); comb400 = data.frame(xx400[time2],yy400[time2]); comb500 = data.frame(xx500[time2],yy500[time2]); comb700 = data.frame(xx700[time2],yy700[time2])

###### connected envelope for 0.05 
xx=list();yy=list()

xx[[1]]=comb70[,1]; xx[[2]]=comb100[,1];xx[[3]]=comb200[,1];xx[[4]]=comb250[,1];xx[[5]]=comb300[,1];xx[[6]]=comb400[,1];xx[[7]]=comb500[,1];xx[[8]]=comb700[,1]
yy[[1]]=comb70[,2]; yy[[2]]=comb100[,2];yy[[3]]=comb200[,2];yy[[4]]=comb250[,2];yy[[5]]=comb300[,2];yy[[6]]=comb400[,2];yy[[7]]=comb500[,2];yy[[8]]=comb700[,2]

int_dir_list=list()
for(p in 1:8){
  int_dir=matrix(0,nrow=dirs,ncol=2)
  
  qln = mqli(cbind(xx[[p]],yy[[p]]),prob=0.05,dirs=dirs)
  
  #### calculate active lines
  clnid=1
  cln=1
  stop=FALSE
  filn=c(1:dirs)
  ss=0
  sln=0
  while (!stop)
  {
    abc=dertln(filn,clnid)
    a=qln[abc[1],]
    b=qln[abc[2],]
    c=qln[abc[3],]
    tt=actln(a,b,c)
    if (ss==0&tt==1) sln=filn[clnid]
    if (tt*ss==1&cln==sln) stop=TRUE
   
    ss=tt
    m=length(filn)
    if (tt==1)
    {
      if (clnid==m) {cln=filn[1] 
      clnid=1} else 
      {cln=filn[clnid+1] 
      clnid=clnid+1}
    } else { 
      clnidt=clnid
      if (clnid==1) {cln=filn[m] 
      clnid=m-1} else 
      {cln=filn[clnid-1] 
      clnid=clnid-1}
      filn=filn[-clnidt]
    }
  } 
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts[i,]=inln(actlns[i,],actlns[j,])
  }
  
  activelines=rep(0,dirs)
  activelines[filn]=filn
  k=1
  ### intersection points corresponding to active directions
  for(i in 1:100){
    if(activelines[i] != 0) {int_dir[i,]=intpts[k,]; k=k+1 }
    
  }
  ### aprroximation of points corresponding to inactive directions
  zero_start=c()
  count=numeric()
  index=1
  
  i=1 
  while(i <= 100){
    count[index]=0
    if(activelines[i]==0){ 
      zero_start[index]=i
      count[index]=count[index]+1
      
      if(i+1<=100){             ### in case the last direction is inactive
        for(j in (i+1):100){
          if(activelines[j]==0){
            count[index]=count[index]+1
            i=i+1
          } else {index= index+1
          break}
        }
      }
    }
    i=i+1
  }
  
  
  for(i in 1:length(zero_start)){
    
    ##taking care of if 1st direction is inactive
    if(zero_start[i]==1){
      x1=intpts[dim(intpts)[1],1]
      y1=intpts[dim(intpts)[1],2]
    }else {
      x1=intpts[which(filn== zero_start[i]-1),1]
      y1=intpts[which(filn==zero_start[i]-1),2]
    }
    
    #taking care of end point
    if(i==length(zero_start) & count[length(count)]!=0){
      x2=intpts[1,1]
      y2=intpts[1,2]
      
    }
    else{
      x2=intpts[which(filn==zero_start[i]+count[i]),1]
      y2=intpts[which(filn==zero_start[i]+count[i]),2]
    }
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
    int_dir[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
  }
  
  int_dir_list[[p]]=int_dir
}


horizontalspeed=rbind(comb70[,1],comb100[,1],comb200[,1],comb250[,1],comb300[,1],comb400[,1],comb500[,1],comb700[,1])
verticalspeed=rbind(comb70[,2],comb100[,2],comb200[,2],comb250[,2], comb300[,2],comb400[,2],comb500[,2],comb700[,2])


pressurelevel=c(70,100,200,250,300,400,500,700)
col=rainbow(64)

library(scatterplot3d)
s=scatterplot3d(horizontalspeed[,1],verticalspeed[,1],pressurelevel,type="p" ,
                xlim = c(-50,100), ylim = c(-100,50), zlim = c(0,800), cex.symbols =0.4, pch = 20,
                xlab = "Horizontal speed", ylab = "",zlab = "Pressure level",
                angle=10, color ="grey60", main="1962-1986")
dims <- par("usr")
x <- dims[1]+ 0.9*diff(dims[1:2])-1
y <- dims[3]+ 0.08*diff(dims[3:4])-0.21
text(x,y,"Vertical speed",srt=10)

for(i in 1: sum(time1)){
  s$points3d(horizontalspeed[,i],verticalspeed[,i],pressurelevel, cex=0.4,pch=20,col="grey60")
}

horizontalintpts=rbind(int_dir_list[[1]][,1], int_dir_list[[2]][,1], int_dir_list[[3]][,1], int_dir_list[[4]][,1], int_dir_list[[5]][,1], int_dir_list[[6]][,1], int_dir_list[[7]][,1], int_dir_list[[8]][,1])
verticalintpts=rbind(int_dir_list[[1]][,2], int_dir_list[[2]][,2], int_dir_list[[3]][,2], int_dir_list[[4]][,2], int_dir_list[[5]][,2], int_dir_list[[6]][,2], int_dir_list[[7]][,2], int_dir_list[[8]][,2])
pressurelevel=c(70,100,200,250,300,400,500,700)

col=rainbow(100)

for(i in 1:100){
  s$points3d(horizontalintpts[,i],verticalintpts[,i],pressurelevel, cex=0.4,pch=20,col=col[i],type = "l")
 
}

for(i in 1:8){
  s$points3d(int_dir_list[[i]][,1],int_dir_list[[i]][,2],rep(pressurelevel[i],100),type="l",lwd=2,col="black")
}




######################## Figure 4 #######################################
# Functional outlier simulation

library(tools)
library(RandomFields)
modelsimple <- RMbiwm(nudiag = c(1, 2), nured = 1, rhored = 0.5, cdiag = c(0.2, 0.2), s = c(0.1, 0.1, 0.1))
t=seq(70, 700, length.out = 25) ## 1D locations
t_scale=(t-min(t))/(max(t)-min(t))

samplesize=100
ucomp=matrix(0,nrow = samplesize,ncol = length(t_scale))
vcomp=matrix(0,nrow = samplesize,ncol = length(t_scale))


for(i in 1:(samplesize)){
  simu <- RFsimulate(modelsimple, x=t_scale)   ##by default all gaussian random fields have zero mean
  
  ucomp[i,]=simu@data[,1]
  vcomp[i,]=simu@data[,2]
  
}
#shape
for(i in (samplesize-1):(samplesize-1)){
  simu <- matrix(RFsimulate(modelsimple, t_scale),length(t_scale),2)
  simu = cbind(simu[,1]+1*sin(t_scale*length(t_scale)*pi),simu[,2]+1*cos(t_scale*length(t_scale)*pi))
  
  ucomp[i,]=simu[,1]
  vcomp[i,]=simu[,2]
}
#shift
for(i in (samplesize):(samplesize)){
  simu <- matrix(RFsimulate(modelsimple, t_scale),length(t_scale),2)
  simu = cbind(simu[,1]+rep(4,length(t_scale)),simu[,2]+rep(0,length(t_scale)))
  
  ucomp[i,]=simu[,1]
  vcomp[i,]=simu[,2]
}


library(scatterplot3d)
##plotting the simulated launches
s=scatterplot3d(ucomp[1,],vcomp[1,],t_scale,color="black", type = "l",
                xlim=c(-5,5),ylim=c(-5,5),xlab = expression(paste(X[1])), ylab = "",zlab = "t", angle=10)

dims <- par("usr")
x <- dims[1]+ 0.9*diff(dims[1:2])-2
y <- dims[3]+ 0.08*diff(dims[3:4]) - 0.35
text(x,y, expression(paste(X[2])), srt=45)


for(i in 2:(samplesize-2)){
  s$points3d(ucomp[i,],vcomp[i,],t_scale,col="black",type = "l")
}
for(i in (samplesize-1):(samplesize-1)){
  s$points3d(ucomp[i,],vcomp[i,],t_scale,col="red",type = "l", lwd=2)
}
for(i in (samplesize):(samplesize)){
  s$points3d(ucomp[i,],vcomp[i,],t_scale,col="blue",type = "l", lwd=2)
}




######################## Figure 5 #######################################
# Boxplots for functional outlier simulation

modelsimple <- RMbiwm(nudiag = c(1, 2), nured = 1, rhored = 0.5, cdiag = c(0.2, 0.2), s = c(0.1, 0.1, 0.1))
t=seq(70, 700, length.out = 25) ## 1D locations
t_scale=(t-min(t))/(max(t)-min(t))

samplesize=1000
ucomp=matrix(0,nrow = samplesize,ncol = length(t_scale))
vcomp=matrix(0,nrow = samplesize,ncol = length(t_scale))


for(i in 1:(samplesize)){
  simu <- RFsimulate(modelsimple, x=t_scale)   ##by default all gaussian random fields have zero mean
  
  ucomp[i,]=simu@data[,1]
  vcomp[i,]=simu@data[,2]
  
}
#shape
for(i in (samplesize-19):(samplesize-10)){
  simu <- matrix(RFsimulate(modelsimple, t_scale),length(t_scale),2)
  simu = cbind(simu[,1]+1*sin(t_scale*length(t_scale)*pi),simu[,2]+1*cos(t_scale*length(t_scale)*pi))
  
  ucomp[i,]=simu[,1]
  vcomp[i,]=simu[,2]
}
#shift
for(i in (samplesize-9):(samplesize)){
  simu <- matrix(RFsimulate(modelsimple, t_scale),length(t_scale),2)
  simu = cbind(simu[,1]+rep(4,length(t_scale)),simu[,2]+rep(0,length(t_scale)))
  
  ucomp[i,]=simu[,1]
  vcomp[i,]=simu[,2]
}


### Figure 5(a)
###################### l2 distances with median launch###################
##  to find magnitude outliers
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
euc_median_launch =vector()
for(j in 1:1000){
  launch= cbind(ucomp[j,],vcomp[j,])
  euc_median_launch[j] = euc.dist(median_launch_sim[1,], launch[1,]) +euc.dist(median_launch_sim[2,], launch[2,])+euc.dist(median_launch_sim[3,], launch[3,])+
    euc.dist(median_launch_sim[4,], launch[4,])+euc.dist(median_launch_sim[5,], launch[5,])+euc.dist(median_launch_sim[6,], launch[6,])+
    euc.dist(median_launch_sim[7,], launch[7,])+euc.dist(median_launch_sim[8,], launch[8,])
}
names(euc_median_launch) = seq(1,1000,1)
summary(euc_median_launch)

bshift = boxplot(euc_median_launch, main="Detecting magnitude outliers", range = 3, ylab=expression(D[M]))
outliers_shift=as.numeric(names(bshift$out))



### Figure 5(b)
## to find shape outliers

### euclidean distance differnce
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

## all the launches
diff_euc_launch = vector()
for(j in 1:1000){
  
  launch= cbind(ucomp[j,],vcomp[j,])
  diff_euc_launch[j] = euc.dist(launch[2,],launch[1,])/(t_scale[2]-t_scale[1]) +euc.dist(launch[3,],launch[2,])/(t_scale[3]-t_scale[2]) +
    euc.dist(launch[4,],launch[3,])/(t_scale[4]-t_scale[3]) +euc.dist(launch[5,],launch[4,])/(t_scale[5]-t_scale[4]) +
    euc.dist(launch[6,],launch[5,])/(t_scale[6]-t_scale[5]) +euc.dist(launch[7,],launch[6,])/(t_scale[7]-t_scale[6]) +
    euc.dist(launch[8,],launch[7,])/(t_scale[8]-t_scale[7])
}
names(diff_euc_launch) =seq(1,1000,1)
summary(diff_euc_launch)
bshape=boxplot(diff_euc_launch, main="Detecting shape outliers", range=3, ylab=expression(D[S]))
outliers_shape=as.numeric(names(bshape$out))






######################## Figure 6 #######################################
### Bivariate Gaussian Random Field Simulation study

library(RandomFields)
t=c(70,100,200,250,300,350,400,500,700) ## 1D locations
probs=0.05
dirs=100

## scaling the locations (pressure level) to 0 to 1
t_scale=(t-min(t))/(max(t)-min(t))



##increasing the range parameter (1/a) from 0 to 1 keeping corelation betweeen two variables same
## a11  =a12=a22
modelsim <- RMbiwm(nudiag = c(1, 2), nured = 1, rhored = 0.5, cdiag = c(0.1, 0.1), s = c(0.1, 0.1, 0.1))

modelsim2 <- RMbiwm(nudiag = c(1, 2), nured = 1, rhored = 0.5, cdiag = c(0.1, 0.1), s = c(0.3, 0.3, 0.3))

modelsim3 <- RMbiwm(nudiag = c(1, 2), nured = 1, rhored = 0.5, cdiag = c(0.1, 0.1), s = c(0.5, 0.5, 0.5))


#######################################################################################################
#######################################################################################################
#######################################################################################################

#### obtaining mean square prediction error by repeating simulations

#range parameter = 0.1 (weak dependence)
samplesize=1000
simulationsize=100

ucomp=matrix(0,nrow = samplesize,ncol = 8)
vcomp=matrix(0,nrow = samplesize,ncol = 8)

u350=rep(0,samplesize)
v350=rep(0,samplesize)
listucomp=list(); listvcomp=list()
listu350=list(); listv350=list()


#### with no outliers
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    simu <- RFsimulate(modelsim, x=t_scale)   ##by default all gaussian random fields have zero mean
    
    u350[i]=simu@data[6,1]
    v350[i]=simu@data[6,2]
    
    ucomp[i,]=simu@data[-6,1]
    vcomp[i,]=simu@data[-6,2]
    
  }
  
  
  listucomp[[k]]=ucomp; listvcomp[[k]]=vcomp
  listu350[[k]]=u350; listv350[[k]]=v350
}


#### with 1% outliers
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    simu <- RFsimulate(modelsim, x=t_scale)   ##by default all gaussian random fields have zero mean
    
    u350[i]=simu@data[6,1]
    v350[i]=simu@data[6,2]
    
    ucomp[i,]=simu@data[-6,1]
    vcomp[i,]=simu@data[-6,2]
    
  }
  
  for(i in (samplesize-9):(samplesize)){
    simu <- matrix(RFsimulate(modelsim, t_scale),9,2)
    simu = cbind(simu[,1]+rep(4,9),simu[,2]+rep(0,9))
    u350[i]=simu[6,1]
    v350[i]=simu[6,2]
    
    ucomp[i,]=simu[-6,1]
    vcomp[i,]=simu[-6,2]
  }
  
  listucomp[[k]]=ucomp; listvcomp[[k]]=vcomp
  listu350[[k]]=u350; listv350[[k]]=v350
}


## with 5% outliers
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    simu <- RFsimulate(modelsim, x=t_scale)   ##by default all gaussian random fields have zero mean
    
    u350[i]=simu@data[6,1]
    v350[i]=simu@data[6,2]
    
    ucomp[i,]=simu@data[-6,1]
    vcomp[i,]=simu@data[-6,2]
    
  }
  
  for(i in (samplesize-49):(samplesize)){
    simu <- matrix(RFsimulate(modelsim, t_scale),9,2)
    simu = cbind(simu[,1]+rep(4,9),simu[,2]+rep(0,9))
    u350[i]=simu[6,1]
    v350[i]=simu[6,2]
    
    ucomp[i,]=simu[-6,1]
    vcomp[i,]=simu[-6,2]
  }
  
  listucomp[[k]]=ucomp; listvcomp[[k]]=vcomp
  listu350[[k]]=u350; listv350[[k]]=v350
}


## with 10% outliers
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    simu <- RFsimulate(modelsim, x=t_scale)   ##by default all gaussian random fields have zero mean
    
    u350[i]=simu@data[6,1]
    v350[i]=simu@data[6,2]
    
    ucomp[i,]=simu@data[-6,1]
    vcomp[i,]=simu@data[-6,2]
    
  }
  
  for(i in (samplesize-99):(samplesize)){
    simu <- matrix(RFsimulate(modelsim, t_scale),9,2)
    simu = cbind(simu[,1]+rep(4,9),simu[,2]+rep(0,9))
    u350[i]=simu[6,1]
    v350[i]=simu[6,2]
    
    ucomp[i,]=simu[-6,1]
    vcomp[i,]=simu[-6,2]
  }
  
  listucomp[[k]]=ucomp; listvcomp[[k]]=vcomp
  listu350[[k]]=u350; listv350[[k]]=v350
}


### the following code can be run for any one case (no outliers/1% outliers/5% outliers/10% outliers) at a time

library(scatterplot3d)
##plotting the simulated launches
s=scatterplot3d(listucomp[[1]][1,],listvcomp[[1]][1,],t[-6],color="black", type = "l",
                xlim=c(-10,10),ylim=c(-10,10),zlim = c(70,800),xlab = "u", ylab = "v",zlab = "Pressure level",angle=10)

for(i in 2:samplesize){
  s$points3d(listucomp[[1]][i,],listvcomp[[1]][i,],t[-6],col="black",type = "l")
}

t_wo350=c(70,100,200,250,300,400,500,700)
t_wo350_scale=t_scale[-6]
############################## predicting conditional means ###################################
sigma22= RFcovmatrix(modelsim,t_wo350_scale) #16X16
sigma= RFcovmatrix(modelsim,c(350,t_wo350_scale))

sigma12=rbind(c(sigma[1,2:9],sigma[1,11:18]), c(sigma[9,2:9],sigma[9,11:18]))

mu=matrix(0,nrow = samplesize, ncol=2)
cmu=matrix(0, nrow = simulationsize, ncol=2)
obsmean=matrix(0, nrow = simulationsize, ncol=2)
truemean = matrix(0, nrow = simulationsize, ncol=2)
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    a=as.matrix(c(listucomp[[k]][i,],listvcomp[[k]][i,]))
    mu[i,]=sigma12%*%solve(sigma22)%*%a
  }
  cmu[k,]=apply(mu,2,mean)
  obsmean[k,]=apply(cbind(listu350[[k]],listv350[[k]]),2,mean)
}

diff_cmu=(cmu-truemean)^2
mse_cm=sum(sqrt(diff_cmu[,1]+diff_cmu[,2]))/simulationsize



### directional quantile envelope at t*=350##########
#####################################################

inln <- function(a,b) ## compute intersection of two lines
{
  inln=-c(b[2]*a[3]-a[2]*b[3],-b[1]*a[3]+a[1]*b[3])/(a[1]*b[2]-a[2]*b[1])
  inln
}

actln <- function(a,b,c) ## tell if a line is active
{
  aa=inln(a,b)
  cc=inln(a,c)
  if (sum((aa-cc)*c[1:2])>0) actln=1 else actln=0
  actln
}

dertln <- function(x,y)  ## determine lines
{
  k=length(x)
  
  if (y==1) bln=x[k] else bln=x[y-1]
  
  if (y==k) aln=x[1] else aln=x[y+1]
  
  cln=x[y]
  
  c(bln,cln,aln)
}

mqli <- function(x, prob, dirs=100) ## generate directions
{
  ang = seq(-pi,pi,length=dirs+1)
  dir = cbind(cos(ang[1:dirs]),sin(ang[1:dirs]))
  cbind(dir,-apply(dir,1,function(u) quantile(x %*% u, prob, type=1)))
}

abcline <- function(a,b,c,...) ## plot lines in terms of coefficients
{
  if (b == 0) abline(v=-c/a,...) else  abline(-c/b,-a/b,...)
}

median350s=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize)
{
  probs=0.495
  stop=FALSE
  while(!stop){
    qln = mqli(cbind(listu350[[k]],listv350[[k]]),prob=probs,dirs=dirs)
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    
    
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]
      b=qln[abc[2],]
      c=qln[abc[3],]
      
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        break}
      
      tt=actln(a,b,c)
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
    
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts350s=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts350s[i,]=inln(actlns[i,],actlns[j,])
  }
  if(sum(abs(intpts350s) > 1e+10 )>=2){
    intpts350s=intpts350s[-which(abs(intpts350s[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median350s[k,]=c(mean(intpts350s[,1]),mean(intpts350s[,2]))
  
}



######################################## kriging #########################################################
median350hk=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize){
  
  probs=0.495
  stop=FALSE
  while(!stop){
    dirs=100
    ###hyperplane kriging
    qln70s = mqli(cbind(listucomp[[k]][,1],listvcomp[[k]][,1]),prob=probs,dirs=dirs); qln100s = mqli(cbind(listucomp[[k]][,2],listvcomp[[k]][,2]),prob=probs,dirs=dirs);qln200s = mqli(cbind(listucomp[[k]][,3],listvcomp[[k]][,3]),prob=probs,dirs=dirs);qln250s = mqli(cbind(listucomp[[k]][,4],listvcomp[[k]][,4]),prob=probs,dirs=dirs);qln300s = mqli(cbind(listucomp[[k]][,5],listvcomp[[k]][,5]),prob=probs,dirs=dirs);qln400s = mqli(cbind(listucomp[[k]][,6],listvcomp[[k]][,6]),prob=probs,dirs=dirs);qln500s = mqli(cbind(listucomp[[k]][,7],listvcomp[[k]][,7]),prob=probs,dirs=dirs);qln700s = mqli(cbind(listucomp[[k]][,8],listvcomp[[k]][,8]),prob=probs,dirs=dirs)
    dirq70s= -qln70s[,3]; dirq100s= -qln100s[,3]; dirq200s= -qln200s[,3]; dirq250s= -qln250s[,3]; dirq300s= -qln300s[,3]; dirq400s= -qln400s[,3]; dirq500s= -qln500s[,3]; dirq700s= -qln700s[,3]
    
    qs=matrix(0,100,length(t_scale[-6]))
    for(s in 1:dirs)
    {
      qs[s,]=c(dirq70s[s],dirq100s[s],dirq200s[s],dirq250s[s],dirq300s[s],dirq400s[s],dirq500s[s],dirq700s[s])
      
    }
    ##univariate kriging with matern covaraince function
    qpred350s=rep(0,dirs)
    #qpred120s=rep(0,dirs)
    for(s in 1:dirs){
      pars.model=RMmatern(nu=NA,scale=NA, var=NA)
      df=data.frame(qs[s,])
      pars <- RFfit(pars.model, x=t_scale[-6], dim = 1, data = df, sub.methods="self")
      
      #print(pars)
      data.df=data.frame(x=t_scale[-6],v1=qs[s,])
      krig=RFinterpolate(pars,x=data.frame(x=0.44444444),data =data.df )
      qpred350s[s]=krig@data$v1
    }
    
    
    qln=cbind(qln70s[,1],qln70s[,2],-qpred350s)
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
      b=qln[abc[2],]
      c=qln[abc[3],]
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        break}
      tt=actln(a,b,c)           ## active line:1 else 0
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts350hk=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts350hk[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts350hk) > 1e+10 )>=2){
    intpts350hk=intpts350hk[-which(abs(intpts350hk[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median350hk[k,]=c(mean(intpts350hk[,1]), mean(intpts350hk[,2]))
}




################################## quantile regression #########################################
library(quantreg)
library(splines)
qdir <- function(dirs, x=NULL)
{
  ang = seq(0,2*pi,length=dirs+1)
  ang = ang[1:dirs]-ang[1]/2-ang[dirs]/2     ##angles ranging from -pi to pi
  cbind(cos(ang),sin(ang))
}


median350qr=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize){
  probs=0.495
  stop=FALSE
  while(!stop){
    ## using quantile regression
    horizontalspeeds=c(listucomp[[k]][,1],listucomp[[k]][,2],listucomp[[k]][,3],listucomp[[k]][,4],listucomp[[k]][,5],listucomp[[k]][,6],listucomp[[k]][,7],listucomp[[k]][,8])
    verticalspeeds=c(listvcomp[[k]][,1],listvcomp[[k]][,2],listvcomp[[k]][,3],listvcomp[[k]][,4],listvcomp[[k]][,5],listvcomp[[k]][,6],listvcomp[[k]][,7],listvcomp[[k]][,8])                     
    
    pressure=c(rep(70,samplesize), rep(100,samplesize),
               rep(200,samplesize), rep(250,samplesize),
               rep(300,samplesize), rep(400,samplesize),
               rep(500,samplesize), rep(700,samplesize))
    
    pressure_scale=c(rep(t_scale[1],samplesize), rep(t_scale[2],samplesize), #skipping 350
                     rep(t_scale[3],samplesize), rep(t_scale[4],samplesize),
                     rep(t_scale[5],samplesize), rep(t_scale[7],samplesize),
                     rep(t_scale[8],samplesize), rep(t_scale[8],samplesize))
    
    xx=horizontalspeeds      
    yy=verticalspeeds
    ttt =pressure_scale #case_2013[,"predicted.value_number_spikes_per_plant"] 
    
    days=c(0.44444444) #days=c(70,100,200,250,300,400,500,700)
    clrs=c('black','black','black')
    pchs=c(3,4,8)
    
    qd = qdir(dirs)
    df=6
    predmat = matrix(0,length(days),dirs)   
    xy = cbind(xx,yy)
    
    
    for (kk in 1:dirs)
      predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))
    ##predicting the tau^th directional quantile
    
    
    for (kk in 1:length(days))
    {
      qln=cbind(qd,-predmat[kk,])
      clnid=1
      cln=1
      stop=FALSE
      filn=c(1:dirs)
      ss=0
      sln=0
      while (!stop)
      {
        abc=dertln(filn,clnid)
        a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
        b=qln[abc[2],]
        c=qln[abc[3],]
        if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
          probs=probs-0.005
          break}
        tt=actln(a,b,c)           ## active line:1 else 0
        if (ss==0&tt==1) sln=filn[clnid]
        if (tt*ss==1&cln==sln) stop=TRUE
        #print('ss,tt,sln,cln')
        #print(c(ss,tt,sln,cln))
        ss=tt
        m=length(filn)
        if (tt==1)
        {
          if (clnid==m) {cln=filn[1] 
          clnid=1} else 
          {cln=filn[clnid+1] 
          clnid=clnid+1}
        } else { 
          clnidt=clnid
          if (clnid==1) {cln=filn[m] 
          clnid=m-1} else 
          {cln=filn[clnid-1] 
          clnid=clnid-1}
          filn=filn[-clnidt]
        }
      } 
    }
    
  }
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts350qr=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts350qr[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts350qr) > 1e+10 )>=2){
    intpts350qr=intpts350qr[-which(abs(intpts350qr[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  
  median350qr[k,]=c(mean(intpts350qr[,1]), mean(intpts350qr[,2]))
}




## comparing predicted medians with estimated medians
qr_m=(truemean-median350qr)^2
hk_m=(truemean-median350hk)^2

sum(sqrt(qr_m[,1]+qr_m[,2]))/(simulationsize)
sum(sqrt(hk_m[,1]+hk_m[,2]))/(simulationsize)






########################################################################################################## 
#### range parameter= 0.3 ##
ucomp=matrix(0,nrow = samplesize,ncol = 8)
vcomp=matrix(0,nrow = samplesize,ncol = 8)

u350=rep(0,samplesize)
v350=rep(0,samplesize)
listucomp=list(); listvcomp=list()
listu350=list(); listv350=list()

## with no outliers
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    simu <- RFsimulate(modelsim2, x=t_scale)   ##by default all gaussian random fields have zero mean
    
    u350[i]=simu@data[6,1]
    v350[i]=simu@data[6,2]
    
    ucomp[i,]=simu@data[-6,1]
    vcomp[i,]=simu@data[-6,2]
    
  }
  
  listucomp[[k]]=ucomp; listvcomp[[k]]=vcomp
  listu350[[k]]=u350; listv350[[k]]=v350
}

## with 1% outliers 
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    simu <- RFsimulate(modelsim2, x=t_scale)   ##by default all gaussian random fields have zero mean
    
    u350[i]=simu@data[6,1]
    v350[i]=simu@data[6,2]
    
    ucomp[i,]=simu@data[-6,1]
    vcomp[i,]=simu@data[-6,2]
    
  }
  
  for(i in (samplesize-9):(samplesize)){
    simu <- matrix(RFsimulate(modelsim2, t_scale),9,2)
    simu = cbind(simu[,1]+rep(4,9),simu[,2]+rep(0,9))
    u350[i]=simu[6,1]
    v350[i]=simu[6,2]
    
    ucomp[i,]=simu[-6,1]
    vcomp[i,]=simu[-6,2]
  }
  
  
  
  listucomp[[k]]=ucomp; listvcomp[[k]]=vcomp
  listu350[[k]]=u350; listv350[[k]]=v350
}

# with 5% outliers
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    simu <- RFsimulate(modelsim2, x=t_scale)   ##by default all gaussian random fields have zero mean
    
    u350[i]=simu@data[6,1]
    v350[i]=simu@data[6,2]
    
    ucomp[i,]=simu@data[-6,1]
    vcomp[i,]=simu@data[-6,2]
    
  }
  
  for(i in (samplesize-49):(samplesize)){
    simu <- matrix(RFsimulate(modelsim2, t_scale),9,2)
    simu = cbind(simu[,1]+rep(4,9),simu[,2]+rep(0,9))
    u350[i]=simu[6,1]
    v350[i]=simu[6,2]
    
    ucomp[i,]=simu[-6,1]
    vcomp[i,]=simu[-6,2]
  }
  
  
  
  listucomp[[k]]=ucomp; listvcomp[[k]]=vcomp
  listu350[[k]]=u350; listv350[[k]]=v350
}

## with 10% outliers
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    simu <- RFsimulate(modelsim2, x=t_scale)   ##by default all gaussian random fields have zero mean
    
    u350[i]=simu@data[6,1]
    v350[i]=simu@data[6,2]
    
    ucomp[i,]=simu@data[-6,1]
    vcomp[i,]=simu@data[-6,2]
    
  }
  
  for(i in (samplesize-49):(samplesize)){
    simu <- matrix(RFsimulate(modelsim2, t_scale),9,2)
    simu = cbind(simu[,1]+rep(4,9),simu[,2]+rep(0,9))
    u350[i]=simu[6,1]
    v350[i]=simu[6,2]
    
    ucomp[i,]=simu[-6,1]
    vcomp[i,]=simu[-6,2]
  }
  
  
  
  listucomp[[k]]=ucomp; listvcomp[[k]]=vcomp
  listu350[[k]]=u350; listv350[[k]]=v350
}



############################## predicting conditional means ###################################
sigma22= RFcovmatrix(modelsim2,t_wo350_scale) #16X16
sigma= RFcovmatrix(modelsim2,c(350,t_wo350_scale))
sigma12=rbind(c(sigma[1,2:9],sigma[1,11:18]), c(sigma[9,2:9],sigma[9,11:18]))

mu=matrix(0,nrow = samplesize, ncol=2)
cmu=matrix(0, nrow = simulationsize, ncol=2)
obsmean=matrix(0, nrow = simulationsize, ncol=2)
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    a=as.matrix(c(listucomp[[k]][i,],listvcomp[[k]][i,]))
    mu[i,]=sigma12%*%solve(sigma22)%*%a
  }
  cmu[k,]=apply(mu,2,mean)
  obsmean[k,]=apply(cbind(listu350[[k]],listv350[[k]]),2,mean)
}

diff_cmu=(cmu-truemean)^2
mse_cm_0.3=sum(sqrt(diff_cmu[,1]+diff_cmu[,2]))/simulationsize


### directional quantile envelope at t*=350##########
#####################################################

inln <- function(a,b) ## compute intersection of two lines
{
  inln=-c(b[2]*a[3]-a[2]*b[3],-b[1]*a[3]+a[1]*b[3])/(a[1]*b[2]-a[2]*b[1])
  inln
}

actln <- function(a,b,c) ## tell if a line is active
{
  aa=inln(a,b)
  cc=inln(a,c)
  if (sum((aa-cc)*c[1:2])>0) actln=1 else actln=0
  actln
}

dertln <- function(x,y)  ## determine lines
{
  k=length(x)
  
  if (y==1) bln=x[k] else bln=x[y-1]
  
  if (y==k) aln=x[1] else aln=x[y+1]
  
  cln=x[y]
  
  c(bln,cln,aln)
}

mqli <- function(x, prob, dirs=100) ## generate directions
{
  ang = seq(-pi,pi,length=dirs+1)
  dir = cbind(cos(ang[1:dirs]),sin(ang[1:dirs]))
  cbind(dir,-apply(dir,1,function(u) quantile(x %*% u, prob, type=1)))
}

abcline <- function(a,b,c,...) ## plot lines in terms of coefficients
{
  if (b == 0) abline(v=-c/a,...) else  abline(-c/b,-a/b,...)
}

median350s0.3=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize)
{
  probs=0.495
  stop=FALSE
  while(!stop){
    qln = mqli(cbind(listu350[[k]],listv350[[k]]),prob=probs,dirs=dirs)
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    
    
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]
      b=qln[abc[2],]
      c=qln[abc[3],]
      
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        break}
      
      tt=actln(a,b,c)
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
    
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts350s=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts350s[i,]=inln(actlns[i,],actlns[j,])
  }
  if(sum(abs(intpts350s) > 1e+10 )>=2){
    intpts350s=intpts350s[-which(abs(intpts350s[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median350s0.3[k,]=c(mean(intpts350s[,1]),mean(intpts350s[,2]))
  
}



######################################## kriging #########################################################
median350hk0.3=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize){
  
  probs=0.495
  stop=FALSE
  while(!stop){
    dirs=100
    ###hyperplane kriging
    qln70s = mqli(cbind(listucomp[[k]][,1],listvcomp[[k]][,1]),prob=probs,dirs=dirs); qln100s = mqli(cbind(listucomp[[k]][,2],listvcomp[[k]][,2]),prob=probs,dirs=dirs);qln200s = mqli(cbind(listucomp[[k]][,3],listvcomp[[k]][,3]),prob=probs,dirs=dirs);qln250s = mqli(cbind(listucomp[[k]][,4],listvcomp[[k]][,4]),prob=probs,dirs=dirs);qln300s = mqli(cbind(listucomp[[k]][,5],listvcomp[[k]][,5]),prob=probs,dirs=dirs);qln400s = mqli(cbind(listucomp[[k]][,6],listvcomp[[k]][,6]),prob=probs,dirs=dirs);qln500s = mqli(cbind(listucomp[[k]][,7],listvcomp[[k]][,7]),prob=probs,dirs=dirs);qln700s = mqli(cbind(listucomp[[k]][,8],listvcomp[[k]][,8]),prob=probs,dirs=dirs)
    dirq70s= -qln70s[,3]; dirq100s= -qln100s[,3]; dirq200s= -qln200s[,3]; dirq250s= -qln250s[,3]; dirq300s= -qln300s[,3]; dirq400s= -qln400s[,3]; dirq500s= -qln500s[,3]; dirq700s= -qln700s[,3]
    
    qs=matrix(0,100,length(t_scale[-6]))
    for(s in 1:dirs)
    {
      qs[s,]=c(dirq70s[s],dirq100s[s],dirq200s[s],dirq250s[s],dirq300s[s],dirq400s[s],dirq500s[s],dirq700s[s])
      
    }
    ##univariate kriging with matern covaraince function
    qpred350s=rep(0,dirs)
    #qpred120s=rep(0,dirs)
    for(s in 1:dirs){
      pars.model=RMmatern(nu=NA,scale=NA, var=NA)
      df=data.frame(qs[s,])
      pars <- RFfit(pars.model, x=t_scale[-6], dim = 1, data = df, sub.methods = "self")
      
      #print(pars)
      data.df=data.frame(x=t_scale[-6],v1=qs[s,])
      krig=RFinterpolate(pars,x=data.frame(x=0.44444444),data =data.df )
      qpred350s[s]=krig@data$v1
    }
    
    
    qln=cbind(qln70s[,1],qln70s[,2],-qpred350s)
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
      b=qln[abc[2],]
      c=qln[abc[3],]
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        break}
      tt=actln(a,b,c)           ## active line:1 else 0
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts350hk=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts350hk[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts350hk) > 1e+10 )>=2){
    intpts350hk=intpts350hk[-which(abs(intpts350hk[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median350hk0.3[k,]=c(mean(intpts350hk[,1]), mean(intpts350hk[,2]))
}



################################## quantile regression #########################################
library(quantreg)
library(splines)
qdir <- function(dirs, x=NULL)
{
  ang = seq(0,2*pi,length=dirs+1)
  ang = ang[1:dirs]-ang[1]/2-ang[dirs]/2     ##angles ranging from -pi to pi
  cbind(cos(ang),sin(ang))
}

median350qr0.3=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize){
  probs=0.495
  stop=FALSE
  while(!stop){
    ## using quantile regression
    horizontalspeeds=c(listucomp[[k]][,1],listucomp[[k]][,2],listucomp[[k]][,3],listucomp[[k]][,4],listucomp[[k]][,5],listucomp[[k]][,6],listucomp[[k]][,7],listucomp[[k]][,8])
    verticalspeeds=c(listvcomp[[k]][,1],listvcomp[[k]][,2],listvcomp[[k]][,3],listvcomp[[k]][,4],listvcomp[[k]][,5],listvcomp[[k]][,6],listvcomp[[k]][,7],listvcomp[[k]][,8])                     
    
    pressure=c(rep(70,samplesize), rep(100,samplesize),
               rep(200,samplesize), rep(250,samplesize),
               rep(300,samplesize), rep(400,samplesize),
               rep(500,samplesize), rep(700,samplesize))
    pressure_scale=c(rep(t_scale[1],samplesize), rep(t_scale[2],samplesize), #skipping 350
                     rep(t_scale[3],samplesize), rep(t_scale[4],samplesize),
                     rep(t_scale[5],samplesize), rep(t_scale[7],samplesize),
                     rep(t_scale[8],samplesize), rep(t_scale[8],samplesize))
    
    
    xx=horizontalspeeds      
    yy=verticalspeeds
    ttt =pressure_scale #case_2013[,"predicted.value_number_spikes_per_plant"] 
    
    days=c(0.44444444) #days=c(70,100,200,250,300,400,500,700)
    clrs=c('black','black','black')
    pchs=c(3,4,8)
    
    qd = qdir(dirs)
    df=6
    predmat = matrix(0,length(days),dirs)   
    xy = cbind(xx,yy)
    
    
    for (kk in 1:dirs)
      predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))
    ##predicting the tau^th directional quantile
    
    for (kk in 1:length(days))
    {
      qln=cbind(qd,-predmat[kk,])
      clnid=1
      cln=1
      stop=FALSE
      filn=c(1:dirs)
      ss=0
      sln=0
      while (!stop)
      {
        abc=dertln(filn,clnid)
        a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
        b=qln[abc[2],]
        c=qln[abc[3],]
        if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
          probs=probs-0.005
          break}
        tt=actln(a,b,c)           ## active line:1 else 0
        if (ss==0&tt==1) sln=filn[clnid]
        if (tt*ss==1&cln==sln) stop=TRUE
        #print('ss,tt,sln,cln')
        #print(c(ss,tt,sln,cln))
        ss=tt
        m=length(filn)
        if (tt==1)
        {
          if (clnid==m) {cln=filn[1] 
          clnid=1} else 
          {cln=filn[clnid+1] 
          clnid=clnid+1}
        } else { 
          clnidt=clnid
          if (clnid==1) {cln=filn[m] 
          clnid=m-1} else 
          {cln=filn[clnid-1] 
          clnid=clnid-1}
          filn=filn[-clnidt]
        }
      } 
    }
    
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts350qr=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts350qr[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts350qr) > 1e+10 )>=2){
    intpts350qr=intpts350qr[-which(abs(intpts350qr[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median350qr0.3[k,]=c(mean(intpts350qr[,1]), mean(intpts350qr[,2]))
}


## comparing predicted medians with estimated medians
qr_m=(truemean-median350qr0.3)^2
hk_m=(truemean-median350hk0.3)^2

sum(sqrt(qr_m[,1]+qr_m[,2]))/(simulationsize)
sum(sqrt(hk_m[,1]+hk_m[,2]))/(simulationsize)




########################################################################################################## 
#### range parameter= 0.5 ##
ucomp=matrix(0,nrow = samplesize,ncol = 8)
vcomp=matrix(0,nrow = samplesize,ncol = 8)

u350=rep(0,samplesize)
v350=rep(0,samplesize)
listucomp=list(); listvcomp=list()
listu350=list(); listv350=list()

## with no outliers
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    simu <- RFsimulate(modelsim3, x=t_scale)   ##by default all gaussian random fields have zero mean
    
    u350[i]=simu@data[6,1]
    v350[i]=simu@data[6,2]
    
    ucomp[i,]=simu@data[-6,1]
    vcomp[i,]=simu@data[-6,2]
    
  }
  
  listucomp[[k]]=ucomp; listvcomp[[k]]=vcomp
  listu350[[k]]=u350; listv350[[k]]=v350
}

### with 1% outliers
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    simu <- RFsimulate(modelsim3, x=t_scale)   ##by default all gaussian random fields have zero mean
    
    u350[i]=simu@data[6,1]
    v350[i]=simu@data[6,2]
    
    ucomp[i,]=simu@data[-6,1]
    vcomp[i,]=simu@data[-6,2]
    
  }
  for(i in (samplesize-9):(samplesize)){
    simu <- matrix(RFsimulate(modelsim3, t_scale),9,2)
    simu = cbind(simu[,1]+rep(4,9),simu[,2]+rep(0,9))
    u350[i]=simu[6,1]
    v350[i]=simu[6,2]
    
    ucomp[i,]=simu[-6,1]
    vcomp[i,]=simu[-6,2]
  }
  
  
  
  listucomp[[k]]=ucomp; listvcomp[[k]]=vcomp
  listu350[[k]]=u350; listv350[[k]]=v350
}

## with 5% outliers
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    simu <- RFsimulate(modelsim3, x=t_scale)   ##by default all gaussian random fields have zero mean
    
    u350[i]=simu@data[6,1]
    v350[i]=simu@data[6,2]
    
    ucomp[i,]=simu@data[-6,1]
    vcomp[i,]=simu@data[-6,2]
    
  }
  for(i in (samplesize-49):(samplesize)){
    simu <- matrix(RFsimulate(modelsim3, t_scale),9,2)
    simu = cbind(simu[,1]+rep(4,9),simu[,2]+rep(0,9))
    u350[i]=simu[6,1]
    v350[i]=simu[6,2]
    
    ucomp[i,]=simu[-6,1]
    vcomp[i,]=simu[-6,2]
  }
  
  
  
  listucomp[[k]]=ucomp; listvcomp[[k]]=vcomp
  listu350[[k]]=u350; listv350[[k]]=v350
}

# with 10% outliers
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    simu <- RFsimulate(modelsim3, x=t_scale)   ##by default all gaussian random fields have zero mean
    
    u350[i]=simu@data[6,1]
    v350[i]=simu@data[6,2]
    
    ucomp[i,]=simu@data[-6,1]
    vcomp[i,]=simu@data[-6,2]
    
  }
  for(i in (samplesize-49):(samplesize)){
    simu <- matrix(RFsimulate(modelsim3, t_scale),9,2)
    simu = cbind(simu[,1]+rep(4,9),simu[,2]+rep(0,9))
    u350[i]=simu[6,1]
    v350[i]=simu[6,2]
    
    ucomp[i,]=simu[-6,1]
    vcomp[i,]=simu[-6,2]
  }
  
  
  
  listucomp[[k]]=ucomp; listvcomp[[k]]=vcomp
  listu350[[k]]=u350; listv350[[k]]=v350
}


############################## predicting conditional means ###################################
sigma22= RFcovmatrix(modelsim3,t_wo350_scale) #16X16
sigma= RFcovmatrix(modelsim3,c(350,t_wo350_scale))
sigma12=rbind(c(sigma[1,2:9],sigma[1,11:18]), c(sigma[9,2:9],sigma[9,11:18]))

mu=matrix(0,nrow = samplesize, ncol=2)
cmu=matrix(0, nrow = simulationsize, ncol=2)
obsmean=matrix(0, nrow = simulationsize, ncol=2)
for(k in 1:simulationsize){
  for(i in 1:samplesize){
    a=as.matrix(c(listucomp[[k]][i,],listvcomp[[k]][i,]))
    mu[i,]=sigma12%*%solve(sigma22)%*%a
  }
  cmu[k,]=apply(mu,2,mean)
  obsmean[k,]=apply(cbind(listu350[[k]],listv350[[k]]),2,mean)
}

diff_cmu=(cmu-truemean)^2
mse_cm_0.5=sum(sqrt(diff_cmu[,1]+diff_cmu[,2]))/simulationsize


### directional quantile envelope at t*=350##########
#####################################################

inln <- function(a,b) ## compute intersection of two lines
{
  inln=-c(b[2]*a[3]-a[2]*b[3],-b[1]*a[3]+a[1]*b[3])/(a[1]*b[2]-a[2]*b[1])
  inln
}

actln <- function(a,b,c) ## tell if a line is active
{
  aa=inln(a,b)
  cc=inln(a,c)
  if (sum((aa-cc)*c[1:2])>0) actln=1 else actln=0
  actln
}

dertln <- function(x,y)  ## determine lines
{
  k=length(x)
  
  if (y==1) bln=x[k] else bln=x[y-1]
  
  if (y==k) aln=x[1] else aln=x[y+1]
  
  cln=x[y]
  
  c(bln,cln,aln)
}

mqli <- function(x, prob, dirs=100) ## generate directions
{
  ang = seq(-pi,pi,length=dirs+1)
  dir = cbind(cos(ang[1:dirs]),sin(ang[1:dirs]))
  cbind(dir,-apply(dir,1,function(u) quantile(x %*% u, prob, type=1)))
}

abcline <- function(a,b,c,...) ## plot lines in terms of coefficients
{
  if (b == 0) abline(v=-c/a,...) else  abline(-c/b,-a/b,...)
}

median350s0.5=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize)
{
  probs=0.495
  stop=FALSE
  while(!stop){
    qln = mqli(cbind(listu350[[k]],listv350[[k]]),prob=probs,dirs=dirs)
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    
    
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]
      b=qln[abc[2],]
      c=qln[abc[3],]
      
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        break}
      
      tt=actln(a,b,c)
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
    
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts350s=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts350s[i,]=inln(actlns[i,],actlns[j,])
  }
  if(sum(abs(intpts350s) > 1e+10 )>=2){
    intpts350s=intpts350s[-which(abs(intpts350s[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median350s0.5[k,]=c(mean(intpts350s[,1]),mean(intpts350s[,2]))
  
}



######################################## kriging #########################################################
median350hk0.5=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize){
  
  probs=0.495
  stop=FALSE
  while(!stop){
    ###hyperplane kriging
    qln70s = mqli(cbind(listucomp[[k]][,1],listvcomp[[k]][,1]),prob=probs,dirs=dirs); qln100s = mqli(cbind(listucomp[[k]][,2],listvcomp[[k]][,2]),prob=probs,dirs=dirs);qln200s = mqli(cbind(listucomp[[k]][,3],listvcomp[[k]][,3]),prob=probs,dirs=dirs);qln250s = mqli(cbind(listucomp[[k]][,4],listvcomp[[k]][,4]),prob=probs,dirs=dirs);qln300s = mqli(cbind(listucomp[[k]][,5],listvcomp[[k]][,5]),prob=probs,dirs=dirs);qln400s = mqli(cbind(listucomp[[k]][,6],listvcomp[[k]][,6]),prob=probs,dirs=dirs);qln500s = mqli(cbind(listucomp[[k]][,7],listvcomp[[k]][,7]),prob=probs,dirs=dirs);qln700s = mqli(cbind(listucomp[[k]][,8],listvcomp[[k]][,8]),prob=probs,dirs=dirs)
    dirq70s= -qln70s[,3]; dirq100s= -qln100s[,3]; dirq200s= -qln200s[,3]; dirq250s= -qln250s[,3]; dirq300s= -qln300s[,3]; dirq400s= -qln400s[,3]; dirq500s= -qln500s[,3]; dirq700s= -qln700s[,3]
    
    qs=matrix(0,100,length(t_scale[-6]))
    for(s in 1:dirs)
    {
      qs[s,]=c(dirq70s[s],dirq100s[s],dirq200s[s],dirq250s[s],dirq300s[s],dirq400s[s],dirq500s[s],dirq700s[s])
      
    }
    ##univariate kriging with matern covaraince function
    qpred350s=rep(0,dirs)
    #qpred120s=rep(0,dirs)
    for(s in 1:dirs){
      pars.model=RMmatern(nu=NA,scale=NA, var=NA)
      df=data.frame(qs[s,])
      pars <- RFfit(pars.model, x=t_scale[-6], dim = 1, data = df, sub.methods = "self")
      
      #print(pars)
      data.df=data.frame(x=t_scale[-6],v1=qs[s,])
      krig=RFinterpolate(pars,x=data.frame(x=0.44444444),data =data.df )
      qpred350s[s]=krig@data$v1
    }
    
    
    qln=cbind(qln70s[,1],qln70s[,2],-qpred350s)
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
      b=qln[abc[2],]
      c=qln[abc[3],]
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        break}
      tt=actln(a,b,c)           ## active line:1 else 0
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts350hk=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts350hk[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts350hk) > 1e+10 )>=2){
    intpts350hk=intpts350hk[-which(abs(intpts350hk[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median350hk0.5[k,]=c(mean(intpts350hk[,1]), mean(intpts350hk[,2]))
}



################################## quantile regression #########################################
library(quantreg)
library(splines)
qdir <- function(dirs, x=NULL)
{
  ang = seq(0,2*pi,length=dirs+1)
  ang = ang[1:dirs]-ang[1]/2-ang[dirs]/2     ##angles ranging from -pi to pi
  cbind(cos(ang),sin(ang))
}

median350qr0.5=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize){
  probs=0.495
  stop=FALSE
  while(!stop){
    ## using quantile regression
    horizontalspeeds=c(listucomp[[k]][,1],listucomp[[k]][,2],listucomp[[k]][,3],listucomp[[k]][,4],listucomp[[k]][,5],listucomp[[k]][,6],listucomp[[k]][,7],listucomp[[k]][,8])
    verticalspeeds=c(listvcomp[[k]][,1],listvcomp[[k]][,2],listvcomp[[k]][,3],listvcomp[[k]][,4],listvcomp[[k]][,5],listvcomp[[k]][,6],listvcomp[[k]][,7],listvcomp[[k]][,8])                     
    
    pressure=c(rep(70,samplesize), rep(100,samplesize),
               rep(200,samplesize), rep(250,samplesize),
               rep(300,samplesize), rep(400,samplesize),
               rep(500,samplesize), rep(700,samplesize))
    pressure_scale=c(rep(t_scale[1],samplesize), rep(t_scale[2],samplesize), #skipping 350
                     rep(t_scale[3],samplesize), rep(t_scale[4],samplesize),
                     rep(t_scale[5],samplesize), rep(t_scale[7],samplesize),
                     rep(t_scale[8],samplesize), rep(t_scale[8],samplesize))
    
    
    xx=horizontalspeeds      
    yy=verticalspeeds
    ttt =pressure_scale #case_2013[,"predicted.value_number_spikes_per_plant"] 
    
    days=c(0.44444444) #days=c(70,100,200,250,300,400,500,700)
    clrs=c('black','black','black')
    pchs=c(3,4,8)
    
    qd = qdir(dirs)
    df=6
    predmat = matrix(0,length(days),dirs)   
    xy = cbind(xx,yy)
    
    
    for (kk in 1:dirs)
      predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))
    ##predicting the tau^th directional quantile
    
    for (kk in 1:length(days))
    {
      qln=cbind(qd,-predmat[kk,])
      clnid=1
      cln=1
      stop=FALSE
      filn=c(1:dirs)
      ss=0
      sln=0
      while (!stop)
      {
        abc=dertln(filn,clnid)
        a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
        b=qln[abc[2],]
        c=qln[abc[3],]
        if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
          probs=probs-0.005
          break}
        tt=actln(a,b,c)           ## active line:1 else 0
        if (ss==0&tt==1) sln=filn[clnid]
        if (tt*ss==1&cln==sln) stop=TRUE
        #print('ss,tt,sln,cln')
        #print(c(ss,tt,sln,cln))
        ss=tt
        m=length(filn)
        if (tt==1)
        {
          if (clnid==m) {cln=filn[1] 
          clnid=1} else 
          {cln=filn[clnid+1] 
          clnid=clnid+1}
        } else { 
          clnidt=clnid
          if (clnid==1) {cln=filn[m] 
          clnid=m-1} else 
          {cln=filn[clnid-1] 
          clnid=clnid-1}
          filn=filn[-clnidt]
        }
      } 
    }
    
  }
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts350qr=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts350qr[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts350qr) > 1e+10 )>=2){
    intpts350qr=intpts350qr[-which(abs(intpts350qr[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  
  median350qr0.5[k,]=c(mean(intpts350qr[,1]), mean(intpts350qr[,2]))
}


## comparing predicted medians with estimated medians
qr_m=(truemean-median350qr0.5)^2
hk_m=(truemean-median350hk0.5)^2

sum(sqrt(qr_m[,1]+qr_m[,2]))/(simulationsize)
sum(sqrt(hk_m[,1]+hk_m[,2]))/(simulationsize)



 ####### gaussian simulation results summary for plot####
cm = c(0.007031, 0.01415, 0.01397)
kriging = c(0.013309, 0.013710, 0.01276)
qr = c(0.02582, 0.01858, 0.01402)

cm1 = c(0.039623, 0.040836, 0.040115)
kriging1 = c(0.012358, 0.014351, 0.014631)
qr1 = c(0.022505, 0.017095, 0.018804)

cm5 = c(0.200174, 0.199150, 0.2004259)
kriging5 = c(0.027588, 0.0297647, 0.028629)
qr5 = c(0.0305115, 0.0292915, 0.029193)

cm10 = c(0.3997267, 0.3990274, 0.4012094)
kriging10 = c(0.0598997, 0.0614743, 0.0610913)
qr10 = c(0.0564853, 0.0538005, 0.0531118)

par(oma = c(4, 1, 1, 1))

par(mfrow=c(2,2))
data = matrix(data= c(cm, kriging, qr), nrow=3, ncol=3, byrow = T)
rownames(data) = c("cm", "kr", "qr")
colnames(data) = c("0.1", "0.3", "0.5")
barplot(data, main="No outliers",
        xlab=expression(paste("Range parameter"), beta), ylab = "MSPE", col=c("darkblue","green","red"),
        beside=TRUE, ylim = c(0,0.4))



data1 = matrix(data= c(cm1, kriging1, qr1), nrow=3, ncol=3, byrow = T)
rownames(data1) = c("cm", "kr", "qr")
colnames(data1) = c("0.1", "0.3", "0.5")
barplot(data1, main="1% outliers",
        xlab=expression(paste("Range parameter"), beta),ylab = "MSPE", col=c("darkblue","green","red"),
        beside=TRUE, ylim=c(0,0.4))


data5 = matrix(data= c(cm5, kriging5, qr5), nrow=3, ncol=3, byrow = T)
rownames(data5) = c("cm", "kr", "qr")
colnames(data5) = c("0.1", "0.3", "0.5")
barplot(data5, main="5% outliers",
        xlab=expression(paste("Range parameter"), beta),ylab = "MSPE", col=c("darkblue","green","red"),
        beside=TRUE, ylim=c(0,0.4))

expression(paste("Range parameter"), beta)

data10 = matrix(data= c(cm10, kriging10, qr10), nrow=3, ncol=3, byrow = T)
rownames(data10) = c("cm", "kr", "qr")
colnames(data10) = c("0.1", "0.3", "0.5")
barplot(data10, main="10% outliers",
        xlab=expression(paste("Range parameter"), beta),ylab = "MSPE", col=c("darkblue","green","red"),
        beside=TRUE, ylim=c(0,0.4))

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", 
       legend = c("Conditional mean", "Kriging", "Quantile regression"), bty="n",
       fill = c("darkblue","green","red"), cex = 0.9, xpd = TRUE, horiz = TRUE, inset = c(0,0))





######################## Figure 7 #######################################
#Application- Comparing observed and predicted directional quntile envelopes

par(mfrow=c(1,2))


### estimating and predicting envelope at p=0.05
#################################################################################
################## connected envelope for 0.05 ##################################
## for pressure level 200

dirs=100

int_dir_t1=matrix(0,nrow=dirs,ncol=2)

qln = mqli(cbind(hst1c[,3],vst1c[,3]),prob=0.05,dirs=dirs)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


################ hyperplane kriging #########
probs = 0.05
int_dir_hk_t1=matrix(0,nrow=dirs,ncol=2)

qln70 = mqli(as.matrix(comb70),prob=probs,dirs=dirs); qln100 = mqli(as.matrix(comb100),prob=probs,dirs=dirs); qln200 = mqli(as.matrix(comb200),prob=probs,dirs=dirs); qln250 = mqli(as.matrix(comb250),prob=probs,dirs=dirs); qln300 = mqli(as.matrix(comb300),prob=probs,dirs=dirs);  qln400 = mqli(as.matrix(comb400),prob=probs,dirs=dirs); qln500 = mqli(as.matrix(comb500),prob=probs,dirs=dirs); qln700 = mqli(as.matrix(comb700),prob=probs,dirs=dirs)
dirq70= -qln70[,3]; dirq100= -qln100[,3]; dirq200= -qln200[,3]; dirq250= -qln250[,3]; dirq300= -qln300[,3]; dirq400= -qln400[,3]; dirq500= -qln500[,3]; dirq700= -qln700[,3]

qs=matrix(0,100,length(t[-3]))
for(s in 1:dirs)
{
  qs[s,]=c(dirq70[s],dirq100[s],dirq250[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
  
}
##univariate kriging with matern covaraince function
qpred_k=rep(0,dirs)

for(s in 1:dirs){
  pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean = NA)
  df=data.frame(qs[s,])
  pars <- RFfit(pars.model, x=t[-3], dim = 1, data = df)
  
  data.df=data.frame(x=t[-3],v1=qs[s,])
  krig=RFinterpolate(pars,x=data.frame(x=t[3]),data =data.df )
  qpred_k[s]=krig@data$v1
}


qln=cbind(qln70[,1],qln70[,2],-qpred_k)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_hk_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_hk_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


#################### quantile regression #############################
int_dir_qr_t1=matrix(0,nrow=dirs,ncol=2)
probs=0.05

horizontalspeeds=c(comb70[,1],comb100[,1],comb250[,1],comb300[,1],comb400[,1],comb500[,1],comb700[,1])
verticalspeeds=c(comb70[,2],comb100[,2],comb250[,2],comb300[,2], comb400[,2],comb500[,2],comb700[,2])                       
pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(250,sum(time1)),rep(300,sum(time1)),rep(400,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))

xx=horizontalspeeds      
yy=verticalspeeds
ttt =pressure 

days=c(t[3]) #days=c(70,100,200,250,300,400,500,700)
clrs=c('black','black','black')
pchs=c(3,4,8)

qd = qdir(dirs)
df=6
predmat = matrix(0,length(days),dirs)   
xy = cbind(xx,yy)

for (kk in 1:dirs){
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))}
##predicting the tau^th directional quantile

qln=cbind(qd,-predmat[1,])

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_qr_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_qr_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}

plot(comb200[,1], comb200[,2],col="grey60", pch=20, cex=0.3, main="1962-1986 (Pressure level 200)",
     xlab="Horizontal speed", ylab="Vertical speed")
polygon(int_dir_t1, border = "blue", lwd=2)
polygon(int_dir_hk_t1, border = "green", lwd=2)
polygon(int_dir_qr_t1, border = "red", lwd=2)
legend("topright", legend=c("Estimated", "Kriging","Quantile regression"),
       col=c("blue", "green", "red"), lwd=2, cex=0.6)




### estimating and predicting envelope at p=0.05
#################################################################################
################## connected envelope for 0.05 ##################################
## for pressure level 300

dirs=100

int_dir_t1=matrix(0,nrow=dirs,ncol=2)

qln = mqli(cbind(hst1c[,5],vst1c[,5]),prob=0.05,dirs=dirs)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


################ hyperplane kriging #########
probs = 0.05
int_dir_hk_t1=matrix(0,nrow=dirs,ncol=2)

qln70 = mqli(as.matrix(comb70),prob=probs,dirs=dirs); qln100 = mqli(as.matrix(comb100),prob=probs,dirs=dirs); qln200 = mqli(as.matrix(comb200),prob=probs,dirs=dirs); qln250 = mqli(as.matrix(comb250),prob=probs,dirs=dirs); qln300 = mqli(as.matrix(comb300),prob=probs,dirs=dirs);  qln400 = mqli(as.matrix(comb400),prob=probs,dirs=dirs); qln500 = mqli(as.matrix(comb500),prob=probs,dirs=dirs); qln700 = mqli(as.matrix(comb700),prob=probs,dirs=dirs)
dirq70= -qln70[,3]; dirq100= -qln100[,3]; dirq200= -qln200[,3]; dirq250= -qln250[,3]; dirq300= -qln300[,3]; dirq400= -qln400[,3]; dirq500= -qln500[,3]; dirq700= -qln700[,3]

qs=matrix(0,100,length(t[-4]))
for(s in 1:dirs)
{
  qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq250[s],dirq400[s],dirq500[s],dirq700[s])
  
}
##univariate kriging with matern covaraince function
qpred_k=rep(0,dirs)

for(s in 1:dirs){
  pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean = NA)
  df=data.frame(qs[s,])
  pars <- RFfit(pars.model, x=t[-5], dim = 1, data = df)
  
  data.df=data.frame(x=t[-5],v1=qs[s,])
  krig=RFinterpolate(pars,x=data.frame(x=t[5]),data =data.df )
  qpred_k[s]=krig@data$v1
}


qln=cbind(qln70[,1],qln70[,2],-qpred_k)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_hk_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_hk_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


#################### quantile regression #############################
int_dir_qr_t1=matrix(0,nrow=dirs,ncol=2)
probs=0.05

horizontalspeeds=c(comb70[,1],comb100[,1],comb200[,1],comb250[,1],comb400[,1],comb500[,1],comb700[,1])
verticalspeeds=c(comb70[,2],comb100[,2],comb200[,2],comb250[,2], comb400[,2],comb500[,2],comb700[,2])                       
pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(200,sum(time1)),rep(250,sum(time1)),rep(400,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))

xx=horizontalspeeds      
yy=verticalspeeds
ttt =pressure 

days=c(t[5]) #days=c(70,100,200,250,300,400,500,700)
clrs=c('black','black','black')
pchs=c(3,4,8)

qd = qdir(dirs)
df=6
predmat = matrix(0,length(days),dirs)   
xy = cbind(xx,yy)


for (kk in 1:dirs){
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))}
##predicting the tau^th directional quantile

qln=cbind(qd,-predmat[1,])

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_qr_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_qr_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}

plot(comb300[,1], comb300[,2],col="grey60", pch=20, cex=0.3, main="1962-1986 (Pressure level 300)",
     xlab="Horizontal speed", ylab="Vertical speed")
polygon(int_dir_t1, border = "blue", lwd=2)
polygon(int_dir_hk_t1, border = "green", lwd=2)
polygon(int_dir_qr_t1, border = "red", lwd=2)
legend("topright", legend=c("Estimated", "Kriging","Quantile regression"),
       col=c("blue", "green", "red"), lwd=2, cex=0.6)







############################# Table 1 #############################
# Simulation motivated by real data
t =c(70,100,200,250,300,400,500,700)
t_wo300 = c(70,100,200,250,400,500,700)
hs = cbind(u70[complete],u100[complete],u200[complete],u250[complete],u300[complete],u400[complete],u500[complete],u700[complete])
vs = cbind(v70[complete],v100[complete],v200[complete],v250[complete], v300[complete],v400[complete],v500[complete],v700[complete])                       


library(quantreg)
library(splines)

indicator = c(rep(0,length(horizontalspeeds)), rep(1,length(verticalspeeds)))
y_speed = c(horizontalspeeds, verticalspeeds)
p = c(pressure, pressure)

model_true= lm(y_speed ~ indicator + bs(p, df=6))

predicted = predict(model_true, list(indicator=c(rep(0,8), rep(1,8)), p=c(t,t)) )
predmatrix = matrix(predicted, nrow=8, ncol=2)


library(RandomFields)
####################### alpha=0 #######################################

simulationsize=10
sample_size=1000
sim_true_mat_ucomp = matrix(0, nrow = sample_size, ncol = 8)
sim_true_mat_vcomp = matrix(0, nrow = sample_size, ncol = 8)

listucomptrue=list(); listvcomptrue=list()
for(k in 1:simulationsize){
  
  for(i in 1:sample_size){
    sim_true_mat_ucomp[i,] = predmatrix[,1] + matrix(rsn(8,alpha=0), nrow = 1, ncol=8)
    sim_true_mat_vcomp[i,] = predmatrix[,2] + matrix(rsn(8,alpha=0), nrow = 1, ncol=8)
  }
  
  
  listucomptrue[[k]]=sim_true_mat_ucomp; listvcomptrue[[k]]=sim_true_mat_vcomp
}

t_wo300=c(70,100,200,250,400,500,700)


############################## predicting conditional means ###################################
cmean_true=matrix(0, nrow = simulationsize, ncol=2)
observedmean_true=matrix(0, nrow = simulationsize, ncol=2)
cm_true= matrix(0, nrow = sample_size, ncol=2)

for(k in 1: simulationsize){
  
  for(i in 1:sample_size){
    
    
    cm_true[i,] = c(krig@data$v1, krig@data$v2)
    
  }
  observedmean_true[k,] = apply(cbind(listucomptrue[[k]][,5],listvcomptrue[[k]][,5]), 2 , mean)
  cmean_true[k,] = apply(cm_true,2,mean)
  
}

diff_truemean = (cmean_true - matrix(rep(predmatrix[5,],10), nrow = simulationsize, ncol=2, byrow = T))^2
mse_truemean = sum(sqrt(diff_truemean[,1]+diff_truemean[,2]))/simulationsize

### directional quantile envelope at t*=300##########
#####################################################

inln <- function(a,b) ## compute intersection of two lines
{
  inln=-c(b[2]*a[3]-a[2]*b[3],-b[1]*a[3]+a[1]*b[3])/(a[1]*b[2]-a[2]*b[1])
  inln
}


actln <- function(a,b,c) ## tell if a line is active
{
  aa=inln(a,b)
  cc=inln(a,c)
  if (sum((aa-cc)*c[1:2])>0) actln=1 else actln=0
  actln
}

dertln <- function(x,y)  ## determine lines
{
  k=length(x)
  
  if (y==1) bln=x[k] else bln=x[y-1]
  
  if (y==k) aln=x[1] else aln=x[y+1]
  
  cln=x[y]
  
  c(bln,cln,aln)
}

mqli <- function(x, prob, dirs=100) ## generate directions
{
  ang = seq(-pi,pi,length=dirs+1)
  dir = cbind(cos(ang[1:dirs]),sin(ang[1:dirs]))
  cbind(dir,-apply(dir,1,function(u) quantile(x %*% u, prob, type=1)))
}

abcline <- function(a,b,c,...) ## plot lines in terms of coefficients
{
  if (b == 0) abline(v=-c/a,...) else  abline(-c/b,-a/b,...)
}
dirs=100

median300strue_inf=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize)
{
  probs=0.495
  stop=FALSE
  while(!stop){
    qln = mqli(cbind(listucomptrue[[k]][,5],listvcomptrue[[k]][,5]),prob=probs,dirs=dirs) ##5th one represents 300
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    
    
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]
      b=qln[abc[2],]
      c=qln[abc[3],]
      
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        break}
      
      tt=actln(a,b,c)
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
    
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts300s=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts300s[i,]=inln(actlns[i,],actlns[j,])
  }
  if(sum(abs(intpts300s) > 1e+10 )>=2){
    intpts300s=intpts300s[-which(abs(intpts300s[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median300strue_inf[k,]=c(mean(intpts300s[,1]),mean(intpts300s[,2]))
  
}


######################################## kriging #########################################################
median300hktrue_inf=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize){
  
  probs=0.495
  stop=FALSE
  while(!stop){
    dirs=100
    ###hyperplane kriging
    qln70s = mqli(cbind(listucomptrue[[k]][,1],listvcomptrue[[k]][,1]),prob=probs,dirs=dirs); qln100s = mqli(cbind(listucomptrue[[k]][,2],listvcomptrue[[k]][,2]),prob=probs,dirs=dirs);qln200s = mqli(cbind(listucomptrue[[k]][,3],listvcomptrue[[k]][,3]),prob=probs,dirs=dirs); qln250s = mqli(cbind(listucomptrue[[k]][,4],listvcomptrue[[k]][,4]),prob=probs,dirs=dirs); qln400s = mqli(cbind(listucomptrue[[k]][,6],listvcomptrue[[k]][,6]),prob=probs,dirs=dirs);qln500s = mqli(cbind(listucomptrue[[k]][,7],listvcomptrue[[k]][,7]),prob=probs,dirs=dirs);qln700s = mqli(cbind(listucomptrue[[k]][,8],listvcomptrue[[k]][,8]),prob=probs,dirs=dirs)
    dirq70s= -qln70s[,3]; dirq100s= -qln100s[,3]; dirq200s= -qln200s[,3]; dirq250s= -qln250s[,3]; dirq400s= -qln400s[,3]; dirq500s= -qln500s[,3]; dirq700s= -qln700s[,3]
    
    qs=matrix(0,100,length(t_wo300))
    for(s in 1:dirs)
    {
      qs[s,]=c(dirq70s[s],dirq100s[s],dirq200s[s],dirq250s[s],dirq400s[s],dirq500s[s],dirq700s[s])
      
    }
    
    ##univariate kriging with matern covaraince function
    qpred300s=rep(0,dirs)
    for(s in 1:dirs){
      pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean=NA)
      df=data.frame(qs[s,])
      pars <- RFfit(pars.model, x=t_wo300, dim = 1, data = df)
      #print(pars)
      data.df=data.frame(x=t_wo300,v1=qs[s,])
      krig=RFinterpolate(pars,x=data.frame(x=300),data =data.df )
      qpred300s[s]=krig@data$v1
    }
    
    
    qln=cbind(qln70s[,1],qln70s[,2],-qpred300s)
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
      b=qln[abc[2],]
      c=qln[abc[3],]
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        break}
      tt=actln(a,b,c)           ## active line:1 else 0
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts300hk=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts300hk[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts300hk) > 1e+10 )>=2){
    intpts300hk=intpts300hk[-which(abs(intpts300hk[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median300hktrue_inf[k,]=c(mean(intpts300hk[,1]), mean(intpts300hk[,2]))
}


################################## quantile regression #########################################
library(quantreg)
library(splines)
qdir <- function(dirs, x=NULL)
{
  ang = seq(0,2*pi,length=dirs+1)
  ang = ang[1:dirs]-ang[1]/2-ang[dirs]/2     ##angles ranging from -pi to pi
  cbind(cos(ang),sin(ang))
}

median300qrtrue_inf=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize){
  probs=0.495
  stop=FALSE
  while(!stop){
    ## using quantile regression
    horizontalspeeds=c(listucomptrue[[k]][,1],listucomptrue[[k]][,2],listucomptrue[[k]][,3],listucomptrue[[k]][,4],listucomptrue[[k]][,6],listucomptrue[[k]][,7],listucomptrue[[k]][,8])
    verticalspeeds=c(listvcomptrue[[k]][,1],listvcomptrue[[k]][,2],listvcomptrue[[k]][,3],listvcomptrue[[k]][,4],listvcomptrue[[k]][,6],listvcomptrue[[k]][,7],listvcomptrue[[k]][,8])                     
    
    pressure=c(rep(70,sample_size), rep(100,sample_size),
               rep(200,sample_size), rep(250,sample_size),
               rep(400,sample_size),
               rep(500,sample_size), rep(700,sample_size))
    
    
    xx=horizontalspeeds      
    yy=verticalspeeds
    ttt =pressure #case_2013[,"predicted.value_number_spikes_per_plant"] 
    
    days=c(300) 
    clrs=c('black','black','black')
    pchs=c(3,4,8)
    
    qd = qdir(dirs)
    df=6
    predmat = matrix(0,length(days),dirs)   
    xy = cbind(xx,yy)
    
    
    for (kk in 1:dirs)
      predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))
    ##predicting the tau^th directional quantile
    
    for (kk in 1:length(days))
    {
      qln=cbind(qd,-predmat[kk,])
      clnid=1
      cln=1
      stop=FALSE
      filn=c(1:dirs)
      ss=0
      sln=0
      while (!stop)
      {
        abc=dertln(filn,clnid)
        a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
        b=qln[abc[2],]
        c=qln[abc[3],]
        if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
          probs=probs-0.005
          break}
        tt=actln(a,b,c)           ## active line:1 else 0
        if (ss==0&tt==1) sln=filn[clnid]
        if (tt*ss==1&cln==sln) stop=TRUE
        #print('ss,tt,sln,cln')
        #print(c(ss,tt,sln,cln))
        ss=tt
        m=length(filn)
        if (tt==1)
        {
          if (clnid==m) {cln=filn[1] 
          clnid=1} else 
          {cln=filn[clnid+1] 
          clnid=clnid+1}
        } else { 
          clnidt=clnid
          if (clnid==1) {cln=filn[m] 
          clnid=m-1} else 
          {cln=filn[clnid-1] 
          clnid=clnid-1}
          filn=filn[-clnidt]
        }
      } 
    }
    
  }
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts300qr=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts300qr[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts300qr) > 1e+10 )>=2){
    intpts300qr=intpts300qr[-which(abs(intpts300qr[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  
  median300qrtrue_inf[k,]=c(mean(intpts300qr[,1]), mean(intpts300qr[,2]))
}



## comparing predicted medians with estimated medians
qr_mtrue_normal=(median300strue_inf-median300qrtrue_inf)^2
hk_mtrue_normal=(median300strue_inf-median300hktrue_inf)^2

sum(sqrt(qr_mtrue_normal[,1]+qr_mtrue_normal[,2]))/simulationsize
sum(sqrt(hk_mtrue_normal[,1]+hk_mtrue_normal[,2]))/simulationsize





####################### alpha =1 #######################################

sample_size=1000
sim_true_mat_ucomp = matrix(0, nrow = sample_size, ncol = 8)
sim_true_mat_vcomp = matrix(0, nrow = sample_size, ncol = 8)

listucomptrue=list(); listvcomptrue=list()
for(k in 1:simulationsize){
  
  for(i in 1:sample_size){
    sim_true_mat_ucomp[i,] = predmatrix[,1] + matrix(rsn(8,alpha = 1), nrow = 1, ncol=8)#matrix(rt(8, 5), nrow = 1, ncol=8)
    sim_true_mat_vcomp[i,] = predmatrix[,2] + matrix(rsn(8,alpha = 1), nrow = 1, ncol=8)#matrix(rt(8, 5), nrow = 1, ncol=8)
  }
  
  
  listucomptrue[[k]]=sim_true_mat_ucomp; listvcomptrue[[k]]=sim_true_mat_vcomp
}

t_wo300=c(70,100,200,250,400,500,700)


############################## predicting conditional means ###################################
cmean_true5=matrix(0, nrow = simulationsize, ncol=2)
observedmean_true5=matrix(0, nrow = simulationsize, ncol=2)
cm_true5= matrix(0, nrow = sample_size, ncol=2)

for(k in 1: simulationsize){
  
  for(i in 1:sample_size){
    
    pars.model <-  RMbiwm(nudiag = c(NA, NA), scale = NA,cdiag = c(NA, NA), rhored = NA)
    df=cbind(listucomptrue[[k]][i,-5],listvcomptrue[[k]][i,-5])
    pars <- RFfit(pars.model, x=t_wo300, dim = 1, data = df)
    data.df = data.frame(x=t_wo300, v1=listucomptrue[[k]][i,-5], v2=listvcomptrue[[k]][i,-5])
    
    krig=RFinterpolate(pars,x=data.frame(x=300),data =data.df )
    
    cm_true5[i,] = c(krig@data$v1, krig@data$v2)
    
  }
  observedmean_true5[k,] = apply(cbind(listucomptrue[[k]][,5],listvcomptrue[[k]][,5]), 2 , mean)
  cmean_true5[k,] = apply(cm_true5,2,mean)
  
}

diff_truemean5 = (cmean_true5 - matrix(rep(predmatrix[5,],10), nrow = simulationsize, ncol=2, byrow = T))^2
mse_truemean5 = sum(sqrt(diff_truemean5[,1]+diff_truemean5[,2]))/simulationsize

### directional quantile envelope at t*=300##########
#####################################################

inln <- function(a,b) ## compute intersection of two lines
{
  inln=-c(b[2]*a[3]-a[2]*b[3],-b[1]*a[3]+a[1]*b[3])/(a[1]*b[2]-a[2]*b[1])
  inln
}


actln <- function(a,b,c) ## tell if a line is active
{
  aa=inln(a,b)
  cc=inln(a,c)
  if (sum((aa-cc)*c[1:2])>0) actln=1 else actln=0
  actln
}

dertln <- function(x,y)  ## determine lines
{
  k=length(x)
  
  if (y==1) bln=x[k] else bln=x[y-1]
  
  if (y==k) aln=x[1] else aln=x[y+1]
  
  cln=x[y]
  
  c(bln,cln,aln)
}

mqli <- function(x, prob, dirs=100) ## generate directions
{
  ang = seq(-pi,pi,length=dirs+1)
  dir = cbind(cos(ang[1:dirs]),sin(ang[1:dirs]))
  cbind(dir,-apply(dir,1,function(u) quantile(x %*% u, prob, type=1)))
}

abcline <- function(a,b,c,...) ## plot lines in terms of coefficients
{
  if (b == 0) abline(v=-c/a,...) else  abline(-c/b,-a/b,...)
}

median300strue_5=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize)
{
  probs=0.495
  stop=FALSE
  while(!stop){
    qln = mqli(cbind(listucomptrue[[k]][,5],listvcomptrue[[k]][,5]),prob=probs,dirs=dirs) ##5th one represents 300
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    
    
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]
      b=qln[abc[2],]
      c=qln[abc[3],]
      
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        break}
      
      tt=actln(a,b,c)
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
    
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts300s=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts300s[i,]=inln(actlns[i,],actlns[j,])
  }
  if(sum(abs(intpts300s) > 1e+10 )>=2){
    intpts300s=intpts300s[-which(abs(intpts300s[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median300strue_5[k,]=c(mean(intpts300s[,1]),mean(intpts300s[,2]))
  
}


######################################## kriging #########################################################
median300hktrue_5=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize){
  
  probs=0.495
  stop=FALSE
  while(!stop){
    dirs=100
    ###hyperplane kriging
    qln70s = mqli(cbind(listucomptrue[[k]][,1],listvcomptrue[[k]][,1]),prob=probs,dirs=dirs); qln100s = mqli(cbind(listucomptrue[[k]][,2],listvcomptrue[[k]][,2]),prob=probs,dirs=dirs);qln200s = mqli(cbind(listucomptrue[[k]][,3],listvcomptrue[[k]][,3]),prob=probs,dirs=dirs); qln250s = mqli(cbind(listucomptrue[[k]][,4],listvcomptrue[[k]][,4]),prob=probs,dirs=dirs); qln400s = mqli(cbind(listucomptrue[[k]][,6],listvcomptrue[[k]][,6]),prob=probs,dirs=dirs);qln500s = mqli(cbind(listucomptrue[[k]][,7],listvcomptrue[[k]][,7]),prob=probs,dirs=dirs);qln700s = mqli(cbind(listucomptrue[[k]][,8],listvcomptrue[[k]][,8]),prob=probs,dirs=dirs)
    dirq70s= -qln70s[,3]; dirq100s= -qln100s[,3]; dirq200s= -qln200s[,3]; dirq250s= -qln250s[,3]; dirq400s= -qln400s[,3]; dirq500s= -qln500s[,3]; dirq700s= -qln700s[,3]
    
    qs=matrix(0,100,length(t_wo300))
    for(s in 1:dirs)
    {
      qs[s,]=c(dirq70s[s],dirq100s[s],dirq200s[s],dirq250s[s],dirq400s[s],dirq500s[s],dirq700s[s])
      
    }
    
    ##univariate kriging with matern covaraince function
    qpred300s=rep(0,dirs)
    for(s in 1:dirs){
      pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean=NA)
      df=data.frame(qs[s,])
      pars <- RFfit(pars.model, x=t_wo300, dim = 1, data = df)
      
      #print(pars)
      data.df=data.frame(x=t_wo300,v1=qs[s,])
      krig=RFinterpolate(pars,x=data.frame(x=300),data =data.df )
      qpred300s[s]=krig@data$v1
    }
    
    
    qln=cbind(qln70s[,1],qln70s[,2],-qpred300s)
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
      b=qln[abc[2],]
      c=qln[abc[3],]
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        break}
      tt=actln(a,b,c)           ## active line:1 else 0
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts300hk=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts300hk[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts300hk) > 1e+10 )>=2){
    intpts300hk=intpts300hk[-which(abs(intpts300hk[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median300hktrue_5[k,]=c(mean(intpts300hk[,1]), mean(intpts300hk[,2]))
}


################################## quantile regression #########################################
library(quantreg)
library(splines)
qdir <- function(dirs, x=NULL)
{
  ang = seq(0,2*pi,length=dirs+1)
  ang = ang[1:dirs]-ang[1]/2-ang[dirs]/2     ##angles ranging from -pi to pi
  cbind(cos(ang),sin(ang))
}

median300qrtrue_5=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize){
  probs=0.495
  stop=FALSE
  while(!stop){
    ## using quantile regression
    horizontalspeeds=c(listucomptrue[[k]][,1],listucomptrue[[k]][,2],listucomptrue[[k]][,3],listucomptrue[[k]][,4],listucomptrue[[k]][,6],listucomptrue[[k]][,7],listucomptrue[[k]][,8])
    verticalspeeds=c(listvcomptrue[[k]][,1],listvcomptrue[[k]][,2],listvcomptrue[[k]][,3],listvcomptrue[[k]][,4],listvcomptrue[[k]][,6],listvcomptrue[[k]][,7],listvcomptrue[[k]][,8])                     
    
    pressure=c(rep(70,sample_size), rep(100,sample_size),
               rep(200,sample_size), rep(250,sample_size),
               rep(400,sample_size),
               rep(500,sample_size), rep(700,sample_size))
    
    
    xx=horizontalspeeds      
    yy=verticalspeeds
    ttt =pressure #case_2013[,"predicted.value_number_spikes_per_plant"] 
    
    days=c(300) 
    clrs=c('black','black','black')
    pchs=c(3,4,8)
    
    qd = qdir(dirs)
    df=6
    predmat = matrix(0,length(days),dirs)   
    xy = cbind(xx,yy)
    
    
    for (kk in 1:dirs)
      predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))
    ##predicting the tau^th directional quantile
    
    for (kk in 1:length(days))
    {
      qln=cbind(qd,-predmat[kk,])
      clnid=1
      cln=1
      stop=FALSE
      filn=c(1:dirs)
      ss=0
      sln=0
      while (!stop)
      {
        abc=dertln(filn,clnid)
        a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
        b=qln[abc[2],]
        c=qln[abc[3],]
        if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
          probs=probs-0.005
          break}
        tt=actln(a,b,c)           ## active line:1 else 0
        if (ss==0&tt==1) sln=filn[clnid]
        if (tt*ss==1&cln==sln) stop=TRUE
        #print('ss,tt,sln,cln')
        #print(c(ss,tt,sln,cln))
        ss=tt
        m=length(filn)
        if (tt==1)
        {
          if (clnid==m) {cln=filn[1] 
          clnid=1} else 
          {cln=filn[clnid+1] 
          clnid=clnid+1}
        } else { 
          clnidt=clnid
          if (clnid==1) {cln=filn[m] 
          clnid=m-1} else 
          {cln=filn[clnid-1] 
          clnid=clnid-1}
          filn=filn[-clnidt]
        }
      } 
    }
    
  }
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts300qr=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts300qr[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts300qr) > 1e+10 )>=2){
    intpts300qr=intpts300qr[-which(abs(intpts300qr[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  
  median300qrtrue_5[k,]=c(mean(intpts300qr[,1]), mean(intpts300qr[,2]))
}



## comparing predicted medians with estimated medians
qr_mtrue5=(median300strue_5-median300qrtrue_5)^2
hk_mtrue5=(median300strue_5-median300hktrue_5)^2

sum(sqrt(qr_mtrue5[,1]+qr_mtrue5[,2]))/simulationsize
sum(sqrt(hk_mtrue5[,1]+hk_mtrue5[,2]))/simulationsize





####################### alpha =2 #######################################


sim_true_mat_ucomp = matrix(0, nrow = sample_size, ncol = 8)
sim_true_mat_vcomp = matrix(0, nrow = sample_size, ncol = 8)

listucomptrue=list(); listvcomptrue=list()
for(k in 1:simulationsize){
  
  for(i in 1:sample_size){
    sim_true_mat_ucomp[i,] = predmatrix[,1] + matrix(rsn(8,alpha = 2), nrow = 1, ncol=8)#matrix(rt(8, 4), nrow = 1, ncol=8)
    sim_true_mat_vcomp[i,] = predmatrix[,2] + matrix(rsn(8,alpha = 2), nrow = 1, ncol=8)#matrix(rt(8, 4), nrow = 1, ncol=8)
  }
  
  
  listucomptrue[[k]]=sim_true_mat_ucomp; listvcomptrue[[k]]=sim_true_mat_vcomp
}

t_wo300=c(70,100,200,250,400,500,700)

############################## predicting conditional means ###################################
cmean_true4=matrix(0, nrow = simulationsize, ncol=2)
observedmean_true4=matrix(0, nrow = simulationsize, ncol=2)
cm_true4= matrix(0, nrow = sample_size, ncol=2)

for(k in 1: simulationsize){
  
  for(i in 1:sample_size){
    
    pars.model <-  RMbiwm(nudiag = c(NA, NA), scale = NA,cdiag = c(NA, NA), rhored = NA)
    df=cbind(listucomptrue[[k]][i,-5],listvcomptrue[[k]][i,-5])
    pars <- RFfit(pars.model, x=t_wo300, dim = 1, data = df)
    data.df = data.frame(x=t_wo300, v1=listucomptrue[[k]][i,-5], v2=listvcomptrue[[k]][i,-5])
    
    krig=RFinterpolate(pars,x=data.frame(x=300),data =data.df )
    
    cm_true4[i,] = c(krig@data$v1, krig@data$v2)
    
  }
  observedmean_true4[k,] = apply(cbind(listucomptrue[[k]][,5],listvcomptrue[[k]][,5]), 2 , mean)
  cmean_true4[k,] = apply(cm_true4,2,mean)
  
}



diff_truemean4 = (cmean_true4 - matrix(rep(predmatrix[5,],10), nrow = simulationsize, ncol=2, byrow = T))^2
mse_truemean4 = sum(sqrt(diff_truemean4[,1]+diff_truemean4[,2]))/simulationsize


### directional quantile envelope at t*=300##########
#####################################################

inln <- function(a,b) ## compute intersection of two lines
{
  inln=-c(b[2]*a[3]-a[2]*b[3],-b[1]*a[3]+a[1]*b[3])/(a[1]*b[2]-a[2]*b[1])
  inln
}


actln <- function(a,b,c) ## tell if a line is active
{
  aa=inln(a,b)
  cc=inln(a,c)
  if (sum((aa-cc)*c[1:2])>0) actln=1 else actln=0
  actln
}

dertln <- function(x,y)  ## determine lines
{
  k=length(x)
  
  if (y==1) bln=x[k] else bln=x[y-1]
  
  if (y==k) aln=x[1] else aln=x[y+1]
  
  cln=x[y]
  
  c(bln,cln,aln)
}

mqli <- function(x, prob, dirs=100) ## generate directions
{
  ang = seq(-pi,pi,length=dirs+1)
  dir = cbind(cos(ang[1:dirs]),sin(ang[1:dirs]))
  cbind(dir,-apply(dir,1,function(u) quantile(x %*% u, prob, type=1)))
}

abcline <- function(a,b,c,...) ## plot lines in terms of coefficients
{
  if (b == 0) abline(v=-c/a,...) else  abline(-c/b,-a/b,...)
}

median300strue_4=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize)
{
  probs=0.495
  stop=FALSE
  while(!stop){
    qln = mqli(cbind(listucomptrue[[k]][,5],listvcomptrue[[k]][,5]),prob=probs,dirs=dirs) ##5th one represents 300
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    
    
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]
      b=qln[abc[2],]
      c=qln[abc[3],]
      
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        break}
      
      tt=actln(a,b,c)
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
    
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts300s=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts300s[i,]=inln(actlns[i,],actlns[j,])
  }
  
  if(sum(abs(intpts300s) > 1e+10 )>=2){
    intpts300s=intpts300s[-which(abs(intpts300s[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median300strue_4[k,]=c(mean(intpts300s[,1]),mean(intpts300s[,2]))
  
}


######################################## kriging #########################################################
median300hktrue_4=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize){
  
  probs=0.495
  stop=FALSE
  while(!stop){
    dirs=100
    ###hyperplane kriging
    qln70s = mqli(cbind(listucomptrue[[k]][,1],listvcomptrue[[k]][,1]),prob=probs,dirs=dirs); qln100s = mqli(cbind(listucomptrue[[k]][,2],listvcomptrue[[k]][,2]),prob=probs,dirs=dirs);qln200s = mqli(cbind(listucomptrue[[k]][,3],listvcomptrue[[k]][,3]),prob=probs,dirs=dirs); qln250s = mqli(cbind(listucomptrue[[k]][,4],listvcomptrue[[k]][,4]),prob=probs,dirs=dirs); qln400s = mqli(cbind(listucomptrue[[k]][,6],listvcomptrue[[k]][,6]),prob=probs,dirs=dirs);qln500s = mqli(cbind(listucomptrue[[k]][,7],listvcomptrue[[k]][,7]),prob=probs,dirs=dirs);qln700s = mqli(cbind(listucomptrue[[k]][,8],listvcomptrue[[k]][,8]),prob=probs,dirs=dirs)
    dirq70s= -qln70s[,3]; dirq100s= -qln100s[,3]; dirq200s= -qln200s[,3]; dirq250s= -qln250s[,3]; dirq400s= -qln400s[,3]; dirq500s= -qln500s[,3]; dirq700s= -qln700s[,3]
    
    qs=matrix(0,100,length(t_wo300))
    for(s in 1:dirs)
    {
      qs[s,]=c(dirq70s[s],dirq100s[s],dirq200s[s],dirq250s[s],dirq400s[s],dirq500s[s],dirq700s[s])
      
    }
    
    ##univariate kriging with matern covaraince function
    qpred300s=rep(0,dirs)
    for(s in 1:dirs){
      pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean=NA)
      df=data.frame(qs[s,])
      pars <- RFfit(pars.model, x=t_wo300, dim = 1, data = df)
      
      #print(pars)
      data.df=data.frame(x=t_wo300,v1=qs[s,])
      krig=RFinterpolate(pars,x=data.frame(x=300),data =data.df )
      qpred300s[s]=krig@data$v1
    }
    
    
    qln=cbind(qln70s[,1],qln70s[,2],-qpred300s)
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
      b=qln[abc[2],]
      c=qln[abc[3],]
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        break}
      tt=actln(a,b,c)           ## active line:1 else 0
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts300hk=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts300hk[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts300hk) > 1e+10 )>=2){
    intpts300hk=intpts300hk[-which(abs(intpts300hk[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median300hktrue_4[k,]=c(mean(intpts300hk[,1]), mean(intpts300hk[,2]))
}


################################## quantile regression #########################################
library(quantreg)
library(splines)
qdir <- function(dirs, x=NULL)
{
  ang = seq(0,2*pi,length=dirs+1)
  ang = ang[1:dirs]-ang[1]/2-ang[dirs]/2     ##angles ranging from -pi to pi
  cbind(cos(ang),sin(ang))
}

median300qrtrue_4=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize){
  probs=0.495
  stop=FALSE
  while(!stop){
    ## using quantile regression
    horizontalspeeds=c(listucomptrue[[k]][,1],listucomptrue[[k]][,2],listucomptrue[[k]][,3],listucomptrue[[k]][,4],listucomptrue[[k]][,6],listucomptrue[[k]][,7],listucomptrue[[k]][,8])
    verticalspeeds=c(listvcomptrue[[k]][,1],listvcomptrue[[k]][,2],listvcomptrue[[k]][,3],listvcomptrue[[k]][,4],listvcomptrue[[k]][,6],listvcomptrue[[k]][,7],listvcomptrue[[k]][,8])                     
    
    pressure=c(rep(70,sample_size), rep(100,sample_size),
               rep(200,sample_size), rep(250,sample_size),
               rep(400,sample_size),
               rep(500,sample_size), rep(700,sample_size))
    
    
    xx=horizontalspeeds      
    yy=verticalspeeds
    ttt =pressure #case_2013[,"predicted.value_number_spikes_per_plant"] 
    
    days=c(300) 
    clrs=c('black','black','black')
    pchs=c(3,4,8)
    
    qd = qdir(dirs)
    df=6
    predmat = matrix(0,length(days),dirs)   
    xy = cbind(xx,yy)
    
    
    for (kk in 1:dirs)
      predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))
    ##predicting the tau^th directional quantile
    
    for (kk in 1:length(days))
    {
      qln=cbind(qd,-predmat[kk,])
      clnid=1
      cln=1
      stop=FALSE
      filn=c(1:dirs)
      ss=0
      sln=0
      while (!stop)
      {
        abc=dertln(filn,clnid)
        a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
        b=qln[abc[2],]
        c=qln[abc[3],]
        if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
          probs=probs-0.005
          break}
        tt=actln(a,b,c)           ## active line:1 else 0
        if (ss==0&tt==1) sln=filn[clnid]
        if (tt*ss==1&cln==sln) stop=TRUE
        #print('ss,tt,sln,cln')
        #print(c(ss,tt,sln,cln))
        ss=tt
        m=length(filn)
        if (tt==1)
        {
          if (clnid==m) {cln=filn[1] 
          clnid=1} else 
          {cln=filn[clnid+1] 
          clnid=clnid+1}
        } else { 
          clnidt=clnid
          if (clnid==1) {cln=filn[m] 
          clnid=m-1} else 
          {cln=filn[clnid-1] 
          clnid=clnid-1}
          filn=filn[-clnidt]
        }
      } 
    }
    
  }
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts300qr=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts300qr[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts300qr) > 1e+10 )>=2){
    intpts300qr=intpts300qr[-which(abs(intpts300qr[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  
  median300qrtrue_4[k,]=c(mean(intpts300qr[,1]), mean(intpts300qr[,2]))
}



## comparing predicted medians with estimated medians
qr_mtrue4=(median300strue_4-median300qrtrue_4)^2
hk_mtrue4=(median300strue_4-median300hktrue_4)^2

sum(sqrt(qr_mtrue4[,1]+qr_mtrue4[,2]))/simulationsize
sum(sqrt(hk_mtrue4[,1]+hk_mtrue4[,2]))/simulationsize




####################### alpha = 3 #######################################

sim_true_mat_ucomp = matrix(0, nrow = sample_size, ncol = 8)
sim_true_mat_vcomp = matrix(0, nrow = sample_size, ncol = 8)

listucomptrue=list(); listvcomptrue=list()
for(k in 1:simulationsize){
  
  for(i in 1:sample_size){
    sim_true_mat_ucomp[i,] = predmatrix[,1] + matrix(rsn(8,alpha = 3), nrow = 1, ncol=8)#matrix(rt(8, 3), nrow = 1, ncol=8)
    sim_true_mat_vcomp[i,] = predmatrix[,2] + matrix(rsn(8,alpha = 3), nrow = 1, ncol=8)#matrix(rt(8, 3), nrow = 1, ncol=8)
  }
  
  
  listucomptrue[[k]]=sim_true_mat_ucomp; listvcomptrue[[k]]=sim_true_mat_vcomp
}

t_wo300=c(70,100,200,250,400,500,700)


############################## predicting conditional means ###################################
cmean_true3=matrix(0, nrow = simulationsize, ncol=2)
observedmean_true3=matrix(0, nrow = simulationsize, ncol=2)
cm_true3 = matrix(0, nrow = sample_size, ncol=2)

for(k in 1: simulationsize){
  
  for(i in 1:sample_size){
    
    pars.model <-  RMbiwm(nudiag = c(NA, NA), scale = NA,cdiag = c(NA, NA), rhored = NA)
    df=cbind(listucomptrue[[k]][i,-5],listvcomptrue[[k]][i,-5])
    pars <- RFfit(pars.model, x=t_wo300, dim = 1, data = df)
    data.df = data.frame(x=t_wo300, v1=listucomptrue[[k]][i,-5], v2=listvcomptrue[[k]][i,-5])
    
    krig=RFinterpolate(pars,x=data.frame(x=300),data =data.df )
    
    cm_true3[i,] = c(krig@data$v1, krig@data$v2)
    
  }
  observedmean_true3[k,] = apply(cbind(listucomptrue[[k]][,5],listvcomptrue[[k]][,5]), 2 , mean)
  cmean_true3[k,] = apply(cm_true3,2,mean)
  
}

diff_truemean3 = (cmean_true3 - matrix(rep(predmatrix[5,],10), nrow = simulationsize, ncol=2, byrow = T))^2
mse_truemean3 = sum(sqrt(diff_truemean3[,1]+diff_truemean3[,2]))/simulationsize

### directional quantile envelope at t*=300##########
#####################################################

inln <- function(a,b) ## compute intersection of two lines
{
  inln=-c(b[2]*a[3]-a[2]*b[3],-b[1]*a[3]+a[1]*b[3])/(a[1]*b[2]-a[2]*b[1])
  inln
}


actln <- function(a,b,c) ## tell if a line is active
{
  aa=inln(a,b)
  cc=inln(a,c)
  if (sum((aa-cc)*c[1:2])>0) actln=1 else actln=0
  actln
}

dertln <- function(x,y)  ## determine lines
{
  k=length(x)
  
  if (y==1) bln=x[k] else bln=x[y-1]
  
  if (y==k) aln=x[1] else aln=x[y+1]
  
  cln=x[y]
  
  c(bln,cln,aln)
}

mqli <- function(x, prob, dirs=100) ## generate directions
{
  ang = seq(-pi,pi,length=dirs+1)
  dir = cbind(cos(ang[1:dirs]),sin(ang[1:dirs]))
  cbind(dir,-apply(dir,1,function(u) quantile(x %*% u, prob, type=1)))
}

abcline <- function(a,b,c,...) ## plot lines in terms of coefficients
{
  if (b == 0) abline(v=-c/a,...) else  abline(-c/b,-a/b,...)
}

median300strue_3=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize)
{
  probs=0.495
  stop=FALSE
  while(!stop){
    qln = mqli(cbind(listucomptrue[[k]][,5],listvcomptrue[[k]][,5]),prob=probs,dirs=dirs) ##5th one represents 300
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    
    
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]
      b=qln[abc[2],]
      c=qln[abc[3],]
      
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        break}
      
      tt=actln(a,b,c)
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
    
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts300s=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts300s[i,]=inln(actlns[i,],actlns[j,])
  }
  if(sum(abs(intpts300s) > 1e+10 )>=2){
    intpts300s=intpts300s[-which(abs(intpts300s[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median300strue_3[k,]=c(mean(intpts300s[,1]),mean(intpts300s[,2]))
  
}


######################################## kriging #########################################################
median300hktrue_3=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize){
  
  probs=0.495
  stop=FALSE
  while(!stop){
    dirs=100
    ###hyperplane kriging
    qln70s = mqli(cbind(listucomptrue[[k]][,1],listvcomptrue[[k]][,1]),prob=probs,dirs=dirs); qln100s = mqli(cbind(listucomptrue[[k]][,2],listvcomptrue[[k]][,2]),prob=probs,dirs=dirs);qln200s = mqli(cbind(listucomptrue[[k]][,3],listvcomptrue[[k]][,3]),prob=probs,dirs=dirs); qln250s = mqli(cbind(listucomptrue[[k]][,4],listvcomptrue[[k]][,4]),prob=probs,dirs=dirs); qln400s = mqli(cbind(listucomptrue[[k]][,6],listvcomptrue[[k]][,6]),prob=probs,dirs=dirs);qln500s = mqli(cbind(listucomptrue[[k]][,7],listvcomptrue[[k]][,7]),prob=probs,dirs=dirs);qln700s = mqli(cbind(listucomptrue[[k]][,8],listvcomptrue[[k]][,8]),prob=probs,dirs=dirs)
    dirq70s= -qln70s[,3]; dirq100s= -qln100s[,3]; dirq200s= -qln200s[,3]; dirq250s= -qln250s[,3]; dirq400s= -qln400s[,3]; dirq500s= -qln500s[,3]; dirq700s= -qln700s[,3]
    
    qs=matrix(0,100,length(t_wo300))
    for(s in 1:dirs)
    {
      qs[s,]=c(dirq70s[s],dirq100s[s],dirq200s[s],dirq250s[s],dirq400s[s],dirq500s[s],dirq700s[s])
      
    }
    
    ##univariate kriging with matern covaraince function
    qpred300s=rep(0,dirs)
    for(s in 1:dirs){
      pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean=NA)
      df=data.frame(qs[s,])
      pars <- RFfit(pars.model, x=t_wo300, dim = 1, data = df)
      #print(pars)
      data.df=data.frame(x=t_wo300,v1=qs[s,])
      krig=RFinterpolate(pars,x=data.frame(x=300),data =data.df )
      qpred300s[s]=krig@data$v1
    }
    
    
    qln=cbind(qln70s[,1],qln70s[,2],-qpred300s)
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
      b=qln[abc[2],]
      c=qln[abc[3],]
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        break}
      tt=actln(a,b,c)           ## active line:1 else 0
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts300hk=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts300hk[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts300hk) > 1e+10 )>=2){
    intpts300hk=intpts300hk[-which(abs(intpts300hk[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median300hktrue_3[k,]=c(mean(intpts300hk[,1]), mean(intpts300hk[,2]))
}


################################## quantile regression #########################################
library(quantreg)
library(splines)
qdir <- function(dirs, x=NULL)
{
  ang = seq(0,2*pi,length=dirs+1)
  ang = ang[1:dirs]-ang[1]/2-ang[dirs]/2     ##angles ranging from -pi to pi
  cbind(cos(ang),sin(ang))
}

median300qrtrue_3=matrix(0,nrow = simulationsize,ncol=2)
for(k in 1:simulationsize){
  probs=0.495
  stop=FALSE
  while(!stop){
    ## using quantile regression
    horizontalspeeds=c(listucomptrue[[k]][,1],listucomptrue[[k]][,2],listucomptrue[[k]][,3],listucomptrue[[k]][,4],listucomptrue[[k]][,6],listucomptrue[[k]][,7],listucomptrue[[k]][,8])
    verticalspeeds=c(listvcomptrue[[k]][,1],listvcomptrue[[k]][,2],listvcomptrue[[k]][,3],listvcomptrue[[k]][,4],listvcomptrue[[k]][,6],listvcomptrue[[k]][,7],listvcomptrue[[k]][,8])                     
    
    pressure=c(rep(70,sample_size), rep(100,sample_size),
               rep(200,sample_size), rep(250,sample_size),
               rep(400,sample_size),
               rep(500,sample_size), rep(700,sample_size))
    
    
    xx=horizontalspeeds      
    yy=verticalspeeds
    ttt =pressure #case_2013[,"predicted.value_number_spikes_per_plant"] 
    
    days=c(300) 
    clrs=c('black','black','black')
    pchs=c(3,4,8)
    
    qd = qdir(dirs)
    df=6
    predmat = matrix(0,length(days),dirs)   
    xy = cbind(xx,yy)
    
    
    for (kk in 1:dirs)
      predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))
    ##predicting the tau^th directional quantile
    
    for (kk in 1:length(days))
    {
      qln=cbind(qd,-predmat[kk,])
      clnid=1
      cln=1
      stop=FALSE
      filn=c(1:dirs)
      ss=0
      sln=0
      while (!stop)
      {
        abc=dertln(filn,clnid)
        a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
        b=qln[abc[2],]
        c=qln[abc[3],]
        if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
          probs=probs-0.005
          break}
        tt=actln(a,b,c)           ## active line:1 else 0
        if (ss==0&tt==1) sln=filn[clnid]
        if (tt*ss==1&cln==sln) stop=TRUE
        #print('ss,tt,sln,cln')
        #print(c(ss,tt,sln,cln))
        ss=tt
        m=length(filn)
        if (tt==1)
        {
          if (clnid==m) {cln=filn[1] 
          clnid=1} else 
          {cln=filn[clnid+1] 
          clnid=clnid+1}
        } else { 
          clnidt=clnid
          if (clnid==1) {cln=filn[m] 
          clnid=m-1} else 
          {cln=filn[clnid-1] 
          clnid=clnid-1}
          filn=filn[-clnidt]
        }
      } 
    }
    
  }
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts300qr=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts300qr[i,]=inln(actlns[i,],actlns[j,])
    
  }
  
  if(sum(abs(intpts300qr) > 1e+10 )>=2){
    intpts300qr=intpts300qr[-which(abs(intpts300qr[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median300qrtrue_3[k,]=c(mean(intpts300qr[,1]), mean(intpts300qr[,2]))
}



## comparing predicted medians with estimated medians
qr_mtrue3=(median300strue_3-median300qrtrue_3)^2
hk_mtrue3=(median300strue_3-median300hktrue_3)^2

sum(sqrt(qr_mtrue3[,1]+qr_mtrue3[,2]))/simulationsize
sum(sqrt(hk_mtrue3[,1]+hk_mtrue3[,2]))/simulationsize







############################# Table 2 #############################

complete=complete.cases(u70,v70,u100,v100,u200,v200,u250,v250,u300,v300,u400,v400,u500,v500,u700,v700)

horz70=u70; horz100=u100; horz200=u200; horz250=u250; horz300=u300;horz400=u400; horz500=u500; horz700=u700; horz850=u850; horz1000=u1000
vert70=v70; vert100=v100; vert200=v200; vert250=v250; vert300=v300;vert400=v400; vert500=v500; vert700=v700; vert850=v850; vert1000=v1000 

xx70=horz70[complete]; xx100=horz100[complete]; xx200=horz200[complete]; xx250=horz250[complete]; xx300=horz300[complete]; xx400=horz400[complete]; xx500=horz500[complete];xx700=horz700[complete]
yy70=vert70[complete]; yy100=vert100[complete]; yy200=vert200[complete]; yy250=vert250[complete]; yy300=vert300[complete]; yy400=vert400[complete]; yy500=vert500[complete];yy700=vert700[complete]
t= c(70,100,200,250,300,400,500,700)

time1=year[complete]<=1986  ## 5259 launches
time2=year[complete]>=1987  ## 6193 launches

##### for time period 1

####### centering the data using tukey median 
library(DepthProc)
comb70 = data.frame(xx70[time1],yy70[time1]); comb100 = data.frame(xx100[time1],yy100[time1]); comb200 = data.frame(xx200[time1],yy200[time1]); comb250 = data.frame(xx250[time1],yy250[time1]); comb300 = data.frame(xx300[time1],yy300[time1]); comb400 = data.frame(xx400[time1],yy400[time1]); comb500 = data.frame(xx500[time1],yy500[time1]); comb700 = data.frame(xx700[time1],yy700[time1])
data_depth70=depth(u=comb70, method ="Tukey" ); data_depth100=depth(u=comb100, method ="Tukey" ); data_depth200=depth(u=comb200, method ="Tukey" ); data_depth250=depth(u=comb250, method ="Tukey" ); data_depth300=depth(u=comb300, method ="Tukey" ); data_depth400=depth(u=comb400, method ="Tukey" ); data_depth500=depth(u=comb500, method ="Tukey" ); data_depth700=depth(u=comb700, method ="Tukey" )
median70=comb70[which(data_depth70 == max(data_depth70))[1],]; median100=comb100[which(data_depth100 == max(data_depth100))[1],]; median200=comb200[which(data_depth200 == max(data_depth200))[1],]; median250=comb250[which(data_depth250 == max(data_depth250))[1],]; median300=comb300[which(data_depth300 == max(data_depth300))[1],]; median400=comb400[which(data_depth400 == max(data_depth400))[1],]; median500=comb500[which(data_depth500 == max(data_depth500))[1],]; median700=comb700[which(data_depth700 == max(data_depth700))[1],]
comb70c= matrix(0, sum(time1), 2);comb100c= matrix(0, sum(time1), 2);comb200c= matrix(0, sum(time1), 2);comb250c= matrix(0, sum(time1), 2);comb300c= matrix(0, sum(time1), 2);comb400c= matrix(0, sum(time1), 2);comb500c= matrix(0, sum(time1), 2);comb700c= matrix(0, sum(time1), 2)
for(i in 1:sum(time1)){
  comb70[i,]= comb70[i,]- median70
  comb100[i,]= comb100[i,]- median100
  comb200[i,]= comb200[i,]- median200
  comb250[i,]= comb250[i,]- median250
  comb300[i,]= comb300[i,]- median300
  comb400[i,]= comb400[i,]- median400
  comb500[i,]= comb500[i,]- median500
  comb700[i,]= comb700[i,]- median700
}

hst1c = cbind(comb70[,1],comb100[,1],comb200[,1],comb250[,1],comb300[,1],comb400[,1],comb500[,1],comb700[,1])
vst1c = cbind(comb70[,2],comb100[,2],comb200[,2],comb250[,2], comb300[,2],comb400[,2],comb500[,2],comb700[,2])

############################## predicting conditional means ###################################
cmean=matrix(0, nrow = length(t)-3, ncol=2)
observedmean=matrix(0, nrow = length(t)-3, ncol=2)
cm = matrix(0, nrow = sum(time1), ncol=2)


### leaving k^th pressure level out
for(k in 2:(length(t)-2)){
  
  for(i in 1:sum(time1)){
    
    pars.model <-  RMbiwm(nudiag = c(NA, NA), scale = NA,cdiag = c(NA, NA), rhored = NA)
    df=cbind(hst1c[i,-k],vst1c[i,-k])
    pars <- RFfit(pars.model, x=t[-k], dim = 1, data = df, sub.methods = "self")
    data.df = data.frame(x=t[-k], v1=hst1c[i,-k], v2=vst1c[i,-k])
    
    krig=RFinterpolate(pars,x=data.frame(x=t[k]),data =data.df )
    
    cm[i,] = c(krig@data$v1, krig@data$v2)
    
  }
  
  cmean[k-1,] = apply(cm,2,mean)
  observedmean[k-1,] = apply(cbind(hst1c[,k],vst1c[,k]),2,mean)
  
}

for(k in 2:(length(t)-2)){
  observedmean[k-1,] = apply(cbind(hst1c[,k],vst1c[,k]),2,mean)
  
}


# comparing with center of the distribution (estimated median)
diff_cmu=(cmean-median[-6,])^2
mse=sum(sqrt(diff_cmu[,1]+diff_cmu[,2]))/(length(t)-3)
mse


### leaving out k^th pressure level #####
median= matrix(0, nrow= length(t)-2, ncol=2)
median_hk= matrix(0, nrow= length(t)-2, ncol=2)
median_qr= matrix(0, nrow= length(t)-2, ncol=2)


### directional quantile envelope at pressure level k ##########
#####################################################

inln <- function(a,b) ## compute intersection of two lines
{
  inln=-c(b[2]*a[3]-a[2]*b[3],-b[1]*a[3]+a[1]*b[3])/(a[1]*b[2]-a[2]*b[1])
  inln
}

actln <- function(a,b,c) ## tell if a line is active
{
  aa=inln(a,b)
  cc=inln(a,c)
  if (sum((aa-cc)*c[1:2])>0) actln=1 else actln=0
  actln
}

dertln <- function(x,y)  ## determine lines
{
  k=length(x)
  
  if (y==1) bln=x[k] else bln=x[y-1]
  
  if (y==k) aln=x[1] else aln=x[y+1]
  
  cln=x[y]
  
  c(bln,cln,aln)
}

mqli <- function(x, prob, dirs=100) ## generate directions
{
  ang = seq(-pi,pi,length=dirs+1)
  dir = cbind(cos(ang[1:dirs]),sin(ang[1:dirs]))
  cbind(dir,-apply(dir,1,function(u) quantile(x %*% u, prob, type=1)))
}

abcline <- function(a,b,c,...) ## plot lines in terms of coefficients
{
  if (b == 0) abline(v=-c/a,...) else  abline(-c/b,-a/b,...)
}



dirs=100
for(k in 2:(length(t)-1)){
  
  probs= 0.5
  
  stop=FALSE
  while(!stop){
    print(probs)
    qln = mqli(cbind(hst1c[,k],vst1c[,k]),prob=probs,dirs=dirs)
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    
    
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]
      b=qln[abc[2],]
      c=qln[abc[3],]
      
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        
        probs=probs-0.005
        break}
      
      tt=actln(a,b,c)
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
    
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts_k=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts_k[i,]=inln(actlns[i,],actlns[j,])
  }
  
  if(sum(abs(intpts_k) > 1e+10 )>=2){
    intpts_k=intpts_k[-which(abs(intpts_k[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median[k-1,]=c(mean(intpts_k[,1]),mean(intpts_k[,2]))
  
}


######################################## kriging #########################################################

for(k in 2:(length(t)-1)){
  
  probs= 0.48
  
  stop=FALSE
  while(!stop){
    probs
    dirs=100
    ###hyperplane kriging
    qln70 = mqli(as.matrix(comb70),prob=probs,dirs=dirs); qln100 = mqli(as.matrix(comb100),prob=probs,dirs=dirs); qln200 = mqli(as.matrix(comb200),prob=probs,dirs=dirs); qln250 = mqli(as.matrix(comb250),prob=probs,dirs=dirs); qln300 = mqli(as.matrix(comb300),prob=probs,dirs=dirs);  qln400 = mqli(as.matrix(comb400),prob=probs,dirs=dirs); qln500 = mqli(as.matrix(comb500),prob=probs,dirs=dirs); qln700 = mqli(as.matrix(comb700),prob=probs,dirs=dirs)
    dirq70= -qln70[,3]; dirq100= -qln100[,3]; dirq200= -qln200[,3]; dirq250= -qln250[,3]; dirq300= -qln300[,3]; dirq400= -qln400[,3]; dirq500= -qln500[,3]; dirq700= -qln700[,3]
    
    qs=matrix(0,100,length(t[-k]))
    for(s in 1:dirs)
    {
      if(k==1){
        qs[s,]=c(dirq100[s],dirq200[s],dirq250[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
      }else if(k==2){
        qs[s,]=c(dirq70[s],dirq200[s],dirq250[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
        
      }else if(k==3){
        qs[s,]=c(dirq70[s],dirq100[s],dirq250[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
        
      }else if(k==4){
        qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
        
      }else if(k==5){
        qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq250[s],dirq400[s],dirq500[s],dirq700[s])
        
      }else if(k==6){
        qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq250[s],dirq300[s],dirq500[s],dirq700[s])
        
      }else if(k==7){
        qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq250[s],dirq300[s],dirq400[s],dirq700[s])
        
      }else{
        qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq250[s],dirq300[s],dirq400[s],dirq500[s])
        
      }
      
    }
    ##univariate kriging with matern covaraince function
    qpred_k=rep(0,dirs)
    
    for(s in 1:dirs){
      pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean=NA)
      df=data.frame(qs[s,])
      pars <- RFfit(pars.model, x=t[-k], dim = 1, data = df, sub.methods=c("self"))
      
      data.df=data.frame(x=t[-k],v1=qs[s,])
      krig=RFinterpolate(pars,x=data.frame(x=t[k]),data =data.df )
      qpred_k[s]=krig@data$v1
    }
    
    
    qln=cbind(qln70[,1],qln70[,2],-qpred_k)
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
      b=qln[abc[2],]
      c=qln[abc[3],]
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        
        break}
      tt=actln(a,b,c)           ## active line:1 else 0
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts_k_hk=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts_k_hk[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts_k_hk) > 1e+10 )>=2){
    intpts_k_hk=intpts_k_hk[-which(abs(intpts_k_hk[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median_hk[k-1,]=c(mean(intpts_k_hk[,1]), mean(intpts_k_hk[,2]))
  
}


################################## quantile regression #########################################
library(quantreg)
library(splines)
qdir <- function(dirs, x=NULL)
{
  ang = seq(0,2*pi,length=dirs+1)
  ang = ang[1:dirs]-ang[1]/2-ang[dirs]/2     ##angles ranging from -pi to pi
  cbind(cos(ang),sin(ang))
}


for(k in 2: (length(t)-1)){
  
  if(k==1){
    horizontalspeeds=c(comb100[,1],comb200[,1],comb250[,1],comb300[,1],comb400[,1],comb500[,1],comb700[,1])
    verticalspeeds=c(comb100[,2],comb200[,2],comb250[,2],comb300[,2], comb400[,2],comb500[,2],comb700[,2])                       
    pressure=c( rep(100,sum(time1)), rep(200,sum(time1)), rep(250,sum(time1)),rep(300,sum(time1)),rep(400,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))
  }else if(k==2){
    horizontalspeeds=c(comb70[,1],comb200[,1],comb250[,1],comb300[,1],comb400[,1],comb500[,1],comb700[,1])
    verticalspeeds=c(comb70[,2],comb200[,2],comb250[,2],comb300[,2], comb400[,2],comb500[,2],comb700[,2])                       
    pressure=c(rep(70,sum(time1)), rep(200,sum(time1)), rep(250,sum(time1)),rep(300,sum(time1)),rep(400,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))
    
  }else if(k==3){
    horizontalspeeds=c(comb70[,1],comb100[,1],comb250[,1],comb300[,1],comb400[,1],comb500[,1],comb700[,1])
    verticalspeeds=c(comb70[,2],comb100[,2],comb250[,2],comb300[,2], comb400[,2],comb500[,2],comb700[,2])                       
    pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(250,sum(time1)),rep(300,sum(time1)),rep(400,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))
    
  }else if(k==4){
    horizontalspeeds=c(comb70[,1],comb100[,1],comb200[,1],comb300[,1],comb400[,1],comb500[,1],comb700[,1])
    verticalspeeds=c(comb70[,2],comb100[,2],comb200[,2],comb300[,2], comb400[,2],comb500[,2],comb700[,2])                       
    pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(200,sum(time1)),rep(300,sum(time1)),rep(400,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))
    
  }else if(k==5){
    horizontalspeeds=c(comb70[,1],comb100[,1],comb200[,1],comb250[,1],comb400[,1],comb500[,1],comb700[,1])
    verticalspeeds=c(comb70[,2],comb100[,2],comb200[,2],comb250[,2], comb400[,2],comb500[,2],comb700[,2])                       
    pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(200,sum(time1)), rep(250,sum(time1)),rep(400,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))
    
  }else if(k==6){
    horizontalspeeds=c(comb70[,1],comb100[,1],comb200[,1],comb250[,1],comb300[,1],comb500[,1],comb700[,1])
    verticalspeeds=c(comb70[,2],comb100[,2],comb200[,2],comb250[,2],comb300[,2],comb500[,2],comb700[,2])                       
    pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(200,sum(time1)), rep(250,sum(time1)),rep(300,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))
    
  }else if(k==7){
    horizontalspeeds=c(comb70[,1],comb100[,1],comb200[,1],comb250[,1],comb300[,1],comb400[,1],comb700[,1])
    verticalspeeds=c(comb70[,2],comb100[,2],comb200[,2],comb250[,2],comb300[,2], comb400[,2],comb700[,2])                       
    pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(200,sum(time1)), rep(250,sum(time1)),rep(300,sum(time1)),rep(400,sum(time1)), rep(700,sum(time1)))
    
  }else{
    horizontalspeeds=c(comb70[,1],comb100[,1],comb200[,1],comb250[,1],comb300[,1],comb400[,1],comb500[,1])
    verticalspeeds=c(comb70[,2],comb100[,2],comb200[,2],comb250[,2],comb300[,2], comb400[,2],comb500[,2])                       
    pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(200,sum(time1)), rep(250,sum(time1)),rep(300,sum(time1)),rep(400,sum(time1)),rep(500,sum(time1)))
    
  }
  
  probs= 0.48
  
  stop=FALSE
  while(!stop){
    ## using quantile regression
    xx=horizontalspeeds      
    yy=verticalspeeds
    ttt =pressure 
    
    days=c(t[k]) #days=c(70,100,200,250,300,400,500,700)
    clrs=c('black','black','black')
    pchs=c(3,4,8)
    
    qd = qdir(dirs)
    df=6
    predmat = matrix(0,length(days),dirs)   
    xy = cbind(xx,yy)
    
    
    for (kk in 1:dirs){
      predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))}
    ##predicting the tau^th directional quantile
    
    qln=cbind(qd,-predmat[1,])
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
      b=qln[abc[2],]
      c=qln[abc[3],]
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        
        break}
      tt=actln(a,b,c)           ## active line:1 else 0
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
    
    
  }
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts_k_qr=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts_k_qr[i,]=inln(actlns[i,],actlns[j,])
    
  }
  
  if(sum(abs(intpts_k_qr) > 1e+10 )>=2){
    intpts_k_qr=intpts_k_qr[-which(abs(intpts_k_qr[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median_qr[k-1,]=c(mean(intpts_k_qr[,1]), mean(intpts_k_qr[,2]))
  
}



## mean square prediction error by leaving one pressure level out at a time
diff_hkc=(median[-6,]-median_hk[-6,])^2
mspe_hkc=sum(sqrt(diff_hkc[,1]+diff_hkc[,2]))/(length(t)-3)

diff_qrc=(median[-6,]-median_qr[-6,])^2
mspe_qrc=sum(sqrt(diff_qrc[,1]+diff_qrc[,2]))/(length(t)-3)

mspe_hkc
mspe_qrc



### for time period 2
####### centering the data using tukey median 
library(DepthProc)
comb70_t2 = data.frame(xx70[time2],yy70[time2]); comb100_t2 = data.frame(xx100[time2],yy100[time2]); comb200_t2 = data.frame(xx200[time2],yy200[time2]); comb250_t2 = data.frame(xx250[time2],yy250[time2]); comb300_t2 = data.frame(xx300[time2],yy300[time2]); comb400_t2 = data.frame(xx400[time2],yy400[time2]); comb500_t2 = data.frame(xx500[time2],yy500[time2]); comb700_t2 = data.frame(xx700[time2],yy700[time2])
data_depth70=depth(u=comb70_t2, method ="Tukey" ); data_depth100=depth(u=comb100_t2, method ="Tukey" ); data_depth200=depth(u=comb200_t2, method ="Tukey" ); data_depth250=depth(u=comb250_t2, method ="Tukey" ); data_depth300=depth(u=comb300_t2, method ="Tukey" ); data_depth400=depth(u=comb400_t2, method ="Tukey" ); data_depth500=depth(u=comb500_t2, method ="Tukey" ); data_depth700=depth(u=comb700_t2, method ="Tukey" )
median70=comb70_t2[which(data_depth70 == max(data_depth70))[1],]; median100=comb100_t2[which(data_depth100 == max(data_depth100))[1],]; median200=comb200_t2[which(data_depth200 == max(data_depth200))[1],]; median250=comb250_t2[which(data_depth250 == max(data_depth250))[1],]; median300=comb300_t2[which(data_depth300 == max(data_depth300))[1],]; median400=comb400_t2[which(data_depth400 == max(data_depth400))[1],]; median500=comb500_t2[which(data_depth500 == max(data_depth500))[1],]; median700=comb700_t2[which(data_depth700 == max(data_depth700))[1],]
for(i in 1:sum(time2)){
  comb70_t2[i,]= comb70_t2[i,]- median70
  comb100_t2[i,]= comb100_t2[i,]- median100
  comb200_t2[i,]= comb200_t2[i,]- median200
  comb250_t2[i,]= comb250_t2[i,]- median250
  comb300_t2[i,]= comb300_t2[i,]- median300
  comb400_t2[i,]= comb400_t2[i,]- median400
  comb500_t2[i,]= comb500_t2[i,]- median500
  comb700_t2[i,]= comb700_t2[i,]- median700
}


hst2c = cbind(comb70_t2[,1],comb100_t2[,1],comb200_t2[,1],comb250_t2[,1],comb300_t2[,1],comb400_t2[,1],comb500_t2[,1],comb700_t2[,1])
vst2c = cbind(comb70_t2[,2],comb100_t2[,2],comb200_t2[,2],comb250_t2[,2], comb300_t2[,2],comb400_t2[,2],comb500_t2[,2],comb700_t2[,2])

plot(comb100_t2[,1], comb100_t2[,2])

#### estimation of envelope and prediction
library(RandomFields)
library(quantreg)
library(splines)

############################## predicting conditional means ###################################
cmean_t2 = matrix(0, nrow = length(t)-3, ncol=2)
observedmean_t2 = matrix(0, nrow = length(t)-3, ncol=2)
cm_t2 = matrix(0, nrow = sum(time2), ncol=2)


### leaving k^th pressure level out
for(k in 2:(length(t)-2)){
  
  for(i in 1:sum(time2)){
    
    pars.model <-  RMbiwm(nudiag = c(NA, NA), scale = NA,cdiag = c(NA, NA), rhored = NA)
    df=cbind(hst2c[i,-k],vst2c[i,-k])
    pars <- RFfit(pars.model, x=t[-k], dim = 1, data = df, sub.methods = "self")
    data.df = data.frame(x=t[-k], v1=hst2c[i,-k], v2=vst2c[i,-k])
    
    krig=RFinterpolate(pars,x=data.frame(x=t[k]),data =data.df )
    
    cm_t2[i,] = c(krig@data$v1, krig@data$v2)
    
  }
  
  cmean_t2[k-1,] = apply(cm_t2,2,mean)
  observedmean_t2[k-1,] = apply(cbind(hst2c[,k],vst2c[,k]),2,mean)
  
}


# comparing with center of the distribution (estimated median)
diff_cmu_t2 = (cmean_t2-median_t2[-6,])^2
mse_t2 = sum(sqrt(diff_cmu_t2[,1]+diff_cmu_t2[,2]))/(length(t)-3)
mse_t2


### leaving out k^th pressure level #####
median_t2= matrix(0, nrow= length(t)-2, ncol=2)
median_hk_t2= matrix(0, nrow= length(t)-2, ncol=2)
median_qr_t2= matrix(0, nrow= length(t)-2, ncol=2)


### directional quantile envelope at pressure level k ##########
#####################################################

inln <- function(a,b) ## compute intersection of two lines
{
  inln=-c(b[2]*a[3]-a[2]*b[3],-b[1]*a[3]+a[1]*b[3])/(a[1]*b[2]-a[2]*b[1])
  inln
}

actln <- function(a,b,c) ## tell if a line is active
{
  aa=inln(a,b)
  cc=inln(a,c)
  if (sum((aa-cc)*c[1:2])>0) actln=1 else actln=0
  actln
}

dertln <- function(x,y)  ## determine lines
{
  k=length(x)
  
  if (y==1) bln=x[k] else bln=x[y-1]
  
  if (y==k) aln=x[1] else aln=x[y+1]
  
  cln=x[y]
  
  c(bln,cln,aln)
}

mqli <- function(x, prob, dirs=100) ## generate directions
{
  ang = seq(-pi,pi,length=dirs+1)
  dir = cbind(cos(ang[1:dirs]),sin(ang[1:dirs]))
  cbind(dir,-apply(dir,1,function(u) quantile(x %*% u, prob, type=1)))
}

abcline <- function(a,b,c,...) ## plot lines in terms of coefficients
{
  if (b == 0) abline(v=-c/a,...) else  abline(-c/b,-a/b,...)
}



dirs=100
for(k in 2:(length(t)-1)){
  
  probs= 0.5
  
  stop=FALSE
  while(!stop){
    print(probs)
    qln = mqli(cbind(hst2c[,k],vst2c[,k]),prob=probs,dirs=dirs)
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    
    
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]
      b=qln[abc[2],]
      c=qln[abc[3],]
      
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        
        probs=probs-0.005
        break}
      
      tt=actln(a,b,c)
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
    
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts_k=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts_k[i,]=inln(actlns[i,],actlns[j,])
  }
  
  if(sum(abs(intpts_k) > 1e+10 )>=2){
    intpts_k=intpts_k[-which(abs(intpts_k[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median_t2[k-1,]=c(mean(intpts_k[,1]),mean(intpts_k[,2]))
  
}


######################################## kriging #########################################################

for(k in 2:(length(t)-1)){
  
  probs= 0.48
  
  stop=FALSE
  while(!stop){
    probs
    dirs=100
    ###hyperplane kriging
    qln70 = mqli(as.matrix(comb70_t2),prob=probs,dirs=dirs); qln100 = mqli(as.matrix(comb100_t2),prob=probs,dirs=dirs); qln200 = mqli(as.matrix(comb200_t2),prob=probs,dirs=dirs); qln250 = mqli(as.matrix(comb250_t2),prob=probs,dirs=dirs); qln300 = mqli(as.matrix(comb300_t2),prob=probs,dirs=dirs);  qln400 = mqli(as.matrix(comb400_t2),prob=probs,dirs=dirs); qln500 = mqli(as.matrix(comb500_t2),prob=probs,dirs=dirs); qln700 = mqli(as.matrix(comb700_t2),prob=probs,dirs=dirs)
    dirq70= -qln70[,3]; dirq100= -qln100[,3]; dirq200= -qln200[,3]; dirq250= -qln250[,3]; dirq300= -qln300[,3]; dirq400= -qln400[,3]; dirq500= -qln500[,3]; dirq700= -qln700[,3]
    
    qs=matrix(0,100,length(t[-k]))
    for(s in 1:dirs)
    {
      if(k==1){
        qs[s,]=c(dirq100[s],dirq200[s],dirq250[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
      }else if(k==2){
        qs[s,]=c(dirq70[s],dirq200[s],dirq250[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
        
      }else if(k==3){
        qs[s,]=c(dirq70[s],dirq100[s],dirq250[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
        
      }else if(k==4){
        qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
        
      }else if(k==5){
        qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq250[s],dirq400[s],dirq500[s],dirq700[s])
        
      }else if(k==6){
        qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq250[s],dirq300[s],dirq500[s],dirq700[s])
        
      }else if(k==7){
        qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq250[s],dirq300[s],dirq400[s],dirq700[s])
        
      }else{
        qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq250[s],dirq300[s],dirq400[s],dirq500[s])
        
      }
      
    }
    ##univariate kriging with matern covaraince function
    qpred_k=rep(0,dirs)
    
    for(s in 1:dirs){
      pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean=NA)
      df=data.frame(qs[s,])
      pars <- RFfit(pars.model, x=t[-k], dim = 1, data = df, sub.methods=c("self"))
      
      data.df=data.frame(x=t[-k],v1=qs[s,])
      krig=RFinterpolate(pars,x=data.frame(x=t[k]),data =data.df )
      qpred_k[s]=krig@data$v1
    }
    
    
    qln=cbind(qln70[,1],qln70[,2],-qpred_k)
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
      b=qln[abc[2],]
      c=qln[abc[3],]
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        
        break}
      tt=actln(a,b,c)           ## active line:1 else 0
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
  }
  
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts_k_hk=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts_k_hk[i,]=inln(actlns[i,],actlns[j,])
    
  }
  if(sum(abs(intpts_k_hk) > 1e+10 )>=2){
    intpts_k_hk=intpts_k_hk[-which(abs(intpts_k_hk[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median_hk_t2[k-1,]=c(mean(intpts_k_hk[,1]), mean(intpts_k_hk[,2]))
  
}


################################## quantile regression #########################################
library(quantreg)
library(splines)
qdir <- function(dirs, x=NULL)
{
  ang = seq(0,2*pi,length=dirs+1)
  ang = ang[1:dirs]-ang[1]/2-ang[dirs]/2     ##angles ranging from -pi to pi
  cbind(cos(ang),sin(ang))
}


for(k in 2: (length(t)-1)){
  
  if(k==1){
    horizontalspeeds=c(comb100_t2[,1],comb200_t2[,1],comb250_t2[,1],comb300_t2[,1],comb400_t2[,1],comb500_t2[,1],comb700_t2[,1])
    verticalspeeds=c(comb100_t2[,2],comb200_t2[,2],comb250_t2[,2],comb300_t2[,2], comb400_t2[,2],comb500_t2[,2],comb700_t2[,2])                       
    pressure=c( rep(100,sum(time2)), rep(200,sum(time2)), rep(250,sum(time2)),rep(300,sum(time2)),rep(400,sum(time2)),rep(500,sum(time2)), rep(700,sum(time2)))
  }else if(k==2){
    horizontalspeeds=c(comb70_t2[,1],comb200_t2[,1],comb250_t2[,1],comb300_t2[,1],comb400_t2[,1],comb500_t2[,1],comb700_t2[,1])
    verticalspeeds=c(comb70_t2[,2],comb200_t2[,2],comb250_t2[,2],comb300_t2[,2], comb400_t2[,2],comb500_t2[,2],comb700_t2[,2])                       
    pressure=c(rep(70,sum(time2)), rep(200,sum(time2)), rep(250,sum(time2)),rep(300,sum(time2)),rep(400,sum(time2)),rep(500,sum(time2)), rep(700,sum(time2)))
    
  }else if(k==3){
    horizontalspeeds=c(comb70_t2[,1],comb100_t2[,1],comb250_t2[,1],comb300_t2[,1],comb400_t2[,1],comb500_t2[,1],comb700_t2[,1])
    verticalspeeds=c(comb70_t2[,2],comb100_t2[,2],comb250_t2[,2],comb300_t2[,2], comb400_t2[,2],comb500_t2[,2],comb700_t2[,2])                       
    pressure=c(rep(70,sum(time2)), rep(100,sum(time2)), rep(250,sum(time2)),rep(300,sum(time2)),rep(400,sum(time2)),rep(500,sum(time2)), rep(700,sum(time2)))
    
  }else if(k==4){
    horizontalspeeds=c(comb70_t2[,1],comb100_t2[,1],comb200_t2[,1],comb300_t2[,1],comb400_t2[,1],comb500_t2[,1],comb700_t2[,1])
    verticalspeeds=c(comb70_t2[,2],comb100_t2[,2],comb200_t2[,2],comb300_t2[,2], comb400_t2[,2],comb500_t2[,2],comb700_t2[,2])                       
    pressure=c(rep(70,sum(time2)), rep(100,sum(time2)), rep(200,sum(time2)),rep(300,sum(time2)),rep(400,sum(time2)),rep(500,sum(time2)), rep(700,sum(time2)))
    
  }else if(k==5){
    horizontalspeeds=c(comb70_t2[,1],comb100_t2[,1],comb200_t2[,1],comb250_t2[,1],comb400_t2[,1],comb500_t2[,1],comb700_t2[,1])
    verticalspeeds=c(comb70_t2[,2],comb100_t2[,2],comb200_t2[,2],comb250_t2[,2], comb400_t2[,2],comb500_t2[,2],comb700_t2[,2])                       
    pressure=c(rep(70,sum(time2)), rep(100,sum(time2)), rep(200,sum(time2)), rep(250,sum(time2)),rep(400,sum(time2)),rep(500,sum(time2)), rep(700,sum(time2)))
    
  }else if(k==6){
    horizontalspeeds=c(comb70_t2[,1],comb100_t2[,1],comb200_t2[,1],comb250_t2[,1],comb300_t2[,1],comb500_t2[,1],comb700_t2[,1])
    verticalspeeds=c(comb70_t2[,2],comb100_t2[,2],comb200_t2[,2],comb250_t2[,2],comb300_t2[,2],comb500_t2[,2],comb700_t2[,2])                       
    pressure=c(rep(70,sum(time2)), rep(100,sum(time2)), rep(200,sum(time2)), rep(250,sum(time2)),rep(300,sum(time2)),rep(500,sum(time2)), rep(700,sum(time2)))
    
  }else if(k==7){
    horizontalspeeds=c(comb70_t2[,1],comb100_t2[,1],comb200_t2[,1],comb250_t2[,1],comb300_t2[,1],comb400_t2[,1],comb700_t2[,1])
    verticalspeeds=c(comb70_t2[,2],comb100_t2[,2],comb200_t2[,2],comb250_t2[,2],comb300_t2[,2], comb400_t2[,2],comb700_t2[,2])                       
    pressure=c(rep(70,sum(time2)), rep(100,sum(time2)), rep(200,sum(time2)), rep(250,sum(time2)),rep(300,sum(time2)),rep(400,sum(time2)), rep(700,sum(time2)))
    
  }else{
    horizontalspeeds=c(comb70_t2[,1],comb100_t2[,1],comb200_t2[,1],comb250_t2[,1],comb300_t2[,1],comb400_t2[,1],comb500_t2[,1])
    verticalspeeds=c(comb70_t2[,2],comb100_t2[,2],comb200_t2[,2],comb250_t2[,2],comb300_t2[,2], comb400_t2[,2],comb500_t2[,2])                       
    pressure=c(rep(70,sum(time2)), rep(100,sum(time2)), rep(200,sum(time2)), rep(250,sum(time2)),rep(300,sum(time2)),rep(400,sum(time2)),rep(500,sum(time2)))
    
  }
  
  probs= 0.48
  
  stop=FALSE
  while(!stop){
    ## using quantile regression
    xx=horizontalspeeds      
    yy=verticalspeeds
    ttt =pressure 
    
    days=c(t[k]) #days=c(70,100,200,250,300,400,500,700)
    clrs=c('black','black','black')
    pchs=c(3,4,8)
    
    qd = qdir(dirs)
    df=6
    predmat = matrix(0,length(days),dirs)   
    xy = cbind(xx,yy)
    
    
    for (kk in 1:dirs){
      predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))}
    ##predicting the tau^th directional quantile
    
    qln=cbind(qd,-predmat[1,])
    clnid=1
    cln=1
    stop=FALSE
    filn=c(1:dirs)
    ss=0
    sln=0
    while (!stop)
    {
      abc=dertln(filn,clnid)
      a=qln[abc[1],]            ##cos  sin  cov corresponding to determined line
      b=qln[abc[2],]
      c=qln[abc[3],]
      if( (abs(a[1])==abs(c[1]) & abs(a[2])== abs(c[2])) | ( abs(a[1])== abs(b[1] & abs(a[2])==abs(b[2])) )) {
        probs=probs-0.005
        
        break}
      tt=actln(a,b,c)           ## active line:1 else 0
      if (ss==0&tt==1) sln=filn[clnid]
      if (tt*ss==1&cln==sln) stop=TRUE
      #print('ss,tt,sln,cln')
      #print(c(ss,tt,sln,cln))
      ss=tt
      m=length(filn)
      if (tt==1)
      {
        if (clnid==m) {cln=filn[1] 
        clnid=1} else 
        {cln=filn[clnid+1] 
        clnid=clnid+1}
      } else { 
        clnidt=clnid
        if (clnid==1) {cln=filn[m] 
        clnid=m-1} else 
        {cln=filn[clnid-1] 
        clnid=clnid-1}
        filn=filn[-clnidt]
      }
    } 
    
    
  }
  actlns=qln[filn,]
  m=dim(actlns)[1]
  intpts_k_qr=matrix(rep(0,2*m),nrow=m)
  for (i in 1:m) 
  {
    if (i==m) j=1 else j=i+1
    intpts_k_qr[i,]=inln(actlns[i,],actlns[j,])
    
  }
  
  if(sum(abs(intpts_k_qr) > 1e+10 )>=2){
    intpts_k_qr=intpts_k_qr[-which(abs(intpts_k_qr[,1])>1e+10),]}  ## taking care of parallel lines, intersection point is infinity
  
  median_qr_t2[k-1,]=c(mean(intpts_k_qr[,1]), mean(intpts_k_qr[,2]))
  
}



## mean square prediction error by leaving one pressure level out at a time
diff_hkc_t2=(median_t2[-6,]-median_hk_t2[-6,])^2
mspe_hkc_t2=sum(sqrt(diff_hkc_t2[,1]+diff_hkc_t2[,2]))/(length(t)-3)

diff_qrc_t2=(median_t2[-6,]-median_qr_t2[-6,])^2
mspe_qrc_t2=sum(sqrt(diff_qrc_t2[,1]+diff_qrc_t2[,2]))/(length(t)-3)

mspe_hkc_t2
mspe_qrc_t2







############################# Table 3 #############################

## time period 1
################## connected envelope for 0.05 ##################################
## for pressure level 100

dirs=100

int_dir_t1=matrix(0,nrow=dirs,ncol=2)

qln = mqli(cbind(hst1c[,2],vst1c[,2]),prob=0.05,dirs=dirs)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


################ hyperplane kriging #########
probs = 0.05
int_dir_hk_t1=matrix(0,nrow=dirs,ncol=2)

qln70 = mqli(as.matrix(comb70),prob=probs,dirs=dirs); qln100 = mqli(as.matrix(comb100),prob=probs,dirs=dirs); qln200 = mqli(as.matrix(comb200),prob=probs,dirs=dirs); qln250 = mqli(as.matrix(comb250),prob=probs,dirs=dirs); qln300 = mqli(as.matrix(comb300),prob=probs,dirs=dirs);  qln400 = mqli(as.matrix(comb400),prob=probs,dirs=dirs); qln500 = mqli(as.matrix(comb500),prob=probs,dirs=dirs); qln700 = mqli(as.matrix(comb700),prob=probs,dirs=dirs)
dirq70= -qln70[,3]; dirq100= -qln100[,3]; dirq200= -qln200[,3]; dirq250= -qln250[,3]; dirq300= -qln300[,3]; dirq400= -qln400[,3]; dirq500= -qln500[,3]; dirq700= -qln700[,3]

qs=matrix(0,100,length(t[-2]))
for(s in 1:dirs)
{
  qs[s,]=c(dirq70[s],dirq200[s],dirq250[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
  
}
##univariate kriging with matern covaraince function
qpred_k=rep(0,dirs)

for(s in 1:dirs){
  pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean = NA)
  df=data.frame(qs[s,])
  pars <- RFfit(pars.model, x=t[-2], dim = 1, data = df, sub.methods = "self")
  
  data.df=data.frame(x=t[-2],v1=qs[s,])
  krig=RFinterpolate(pars,x=data.frame(x=t[2]),data =data.df )
  qpred_k[s]=krig@data$v1
}


qln=cbind(qln70[,1],qln70[,2],-qpred_k)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_hk_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_hk_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


#################### quantile regression #############################
int_dir_qr_t1=matrix(0,nrow=dirs,ncol=2)
probs=0.05

horizontalspeeds=c(comb70[,1],comb200[,1],comb250[,1],comb300[,1],comb400[,1],comb500[,1],comb700[,1])
verticalspeeds=c(comb70[,2],comb200[,2],comb250[,2],comb300[,2], comb400[,2],comb500[,2],comb700[,2])                       
pressure=c(rep(70,sum(time1)), rep(200,sum(time1)), rep(250,sum(time1)),rep(300,sum(time1)),rep(400,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))

xx=horizontalspeeds      
yy=verticalspeeds
ttt =pressure 

days=c(t[2]) #days=c(70,100,200,250,300,400,500,700)
clrs=c('black','black','black')
pchs=c(3,4,8)

qd = qdir(dirs)
df=6
predmat = matrix(0,length(days),dirs)   
xy = cbind(xx,yy)


for (kk in 1:dirs){
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))}
##predicting the tau^th directional quantile

qln=cbind(qd,-predmat[1,])

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_qr_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_qr_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}



diff_hk_100_t1c=(int_dir_t1-int_dir_hk_t1)^2
mse_hk_100_t1c=sum(sqrt(diff_hk_100_t1c[,1]+diff_hk_100_t1c[,2]))/dirs

diff_qr_100_t1c=(int_dir_t1-int_dir_qr_t1)^2
mse_qr_100_t1c=sum(sqrt(diff_qr_100_t1c[,1]+diff_qr_100_t1c[,2]))/dirs


## for pressure level 200

dirs=100

int_dir_t1=matrix(0,nrow=dirs,ncol=2)

qln = mqli(cbind(hst1c[,3],vst1c[,3]),prob=0.05,dirs=dirs)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


################ hyperplane kriging #########
probs = 0.05
int_dir_hk_t1=matrix(0,nrow=dirs,ncol=2)

qln70 = mqli(as.matrix(comb70),prob=probs,dirs=dirs); qln100 = mqli(as.matrix(comb100),prob=probs,dirs=dirs); qln200 = mqli(as.matrix(comb200),prob=probs,dirs=dirs); qln250 = mqli(as.matrix(comb250),prob=probs,dirs=dirs); qln300 = mqli(as.matrix(comb300),prob=probs,dirs=dirs);  qln400 = mqli(as.matrix(comb400),prob=probs,dirs=dirs); qln500 = mqli(as.matrix(comb500),prob=probs,dirs=dirs); qln700 = mqli(as.matrix(comb700),prob=probs,dirs=dirs)
dirq70= -qln70[,3]; dirq100= -qln100[,3]; dirq200= -qln200[,3]; dirq250= -qln250[,3]; dirq300= -qln300[,3]; dirq400= -qln400[,3]; dirq500= -qln500[,3]; dirq700= -qln700[,3]

qs=matrix(0,100,length(t[-3]))
for(s in 1:dirs)
{
  qs[s,]=c(dirq70[s],dirq100[s],dirq250[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
  
}
##univariate kriging with matern covaraince function
qpred_k=rep(0,dirs)

for(s in 1:dirs){
  pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean = NA)
  df=data.frame(qs[s,])
  pars <- RFfit(pars.model, x=t[-3], dim = 1, data = df)
  
  data.df=data.frame(x=t[-3],v1=qs[s,])
  krig=RFinterpolate(pars,x=data.frame(x=t[3]),data =data.df )
  qpred_k[s]=krig@data$v1
}


qln=cbind(qln70[,1],qln70[,2],-qpred_k)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_hk_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_hk_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


#################### quantile regression #############################
int_dir_qr_t1=matrix(0,nrow=dirs,ncol=2)
probs=0.05

horizontalspeeds=c(comb70[,1],comb100[,1],comb250[,1],comb300[,1],comb400[,1],comb500[,1],comb700[,1])
verticalspeeds=c(comb70[,2],comb100[,2],comb250[,2],comb300[,2], comb400[,2],comb500[,2],comb700[,2])                       
pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(250,sum(time1)),rep(300,sum(time1)),rep(400,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))

xx=horizontalspeeds      
yy=verticalspeeds
ttt =pressure 

days=c(t[3]) #days=c(70,100,200,250,300,400,500,700)
clrs=c('black','black','black')
pchs=c(3,4,8)

qd = qdir(dirs)
df=6
predmat = matrix(0,length(days),dirs)   
xy = cbind(xx,yy)

for (kk in 1:dirs){
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))}
##predicting the tau^th directional quantile

qln=cbind(qd,-predmat[1,])

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_qr_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_qr_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


diff_hk_200_t1c=(int_dir_t1-int_dir_hk_t1)^2
mse_hk_200_t1c=sum(sqrt(diff_hk_200_t1c[,1]+diff_hk_200_t1c[,2]))/dirs

diff_qr_200_t1c=(int_dir_t1-int_dir_qr_t1)^2
mse_qr_200_t1c=sum(sqrt(diff_qr_200_t1c[,1]+diff_qr_200_t1c[,2]))/dirs


## for pressure level 250

dirs=100

int_dir_t1=matrix(0,nrow=dirs,ncol=2)

qln = mqli(cbind(hst1c[,4],vst1c[,4]),prob=0.05,dirs=dirs)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


################ hyperplane kriging #########
probs = 0.05
int_dir_hk_t1=matrix(0,nrow=dirs,ncol=2)

qln70 = mqli(as.matrix(comb70),prob=probs,dirs=dirs); qln100 = mqli(as.matrix(comb100),prob=probs,dirs=dirs); qln200 = mqli(as.matrix(comb200),prob=probs,dirs=dirs); qln250 = mqli(as.matrix(comb250),prob=probs,dirs=dirs); qln300 = mqli(as.matrix(comb300),prob=probs,dirs=dirs);  qln400 = mqli(as.matrix(comb400),prob=probs,dirs=dirs); qln500 = mqli(as.matrix(comb500),prob=probs,dirs=dirs); qln700 = mqli(as.matrix(comb700),prob=probs,dirs=dirs)
dirq70= -qln70[,3]; dirq100= -qln100[,3]; dirq200= -qln200[,3]; dirq250= -qln250[,3]; dirq300= -qln300[,3]; dirq400= -qln400[,3]; dirq500= -qln500[,3]; dirq700= -qln700[,3]

qs=matrix(0,100,length(t[-4]))
for(s in 1:dirs)
{
  qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
  
}
##univariate kriging with matern covaraince function
qpred_k=rep(0,dirs)

for(s in 1:dirs){
  pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean = NA)
  df=data.frame(qs[s,])
  pars <- RFfit(pars.model, x=t[-4], dim = 1, data = df)
  
  data.df=data.frame(x=t[-4],v1=qs[s,])
  krig=RFinterpolate(pars,x=data.frame(x=t[4]),data =data.df )
  qpred_k[s]=krig@data$v1
}


qln=cbind(qln70[,1],qln70[,2],-qpred_k)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_hk_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_hk_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


#################### quantile regression #############################
int_dir_qr_t1=matrix(0,nrow=dirs,ncol=2)
probs=0.05

horizontalspeeds=c(comb70[,1],comb100[,1],comb200[,1],comb300[,1],comb400[,1],comb500[,1],comb700[,1])
verticalspeeds=c(comb70[,2],comb100[,2],comb200[,2],comb300[,2], comb400[,2],comb500[,2],comb700[,2])                       
pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(200,sum(time1)),rep(300,sum(time1)),rep(400,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))

xx=horizontalspeeds      
yy=verticalspeeds
ttt =pressure 

days=c(t[4]) #days=c(70,100,200,250,300,400,500,700)
clrs=c('black','black','black')
pchs=c(3,4,8)

qd = qdir(dirs)
df=6
predmat = matrix(0,length(days),dirs)   
xy = cbind(xx,yy)


for (kk in 1:dirs){
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))}
##predicting the tau^th directional quantile

qln=cbind(qd,-predmat[1,])

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_qr_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_qr_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


diff_hk_250_t1c=(int_dir_t1-int_dir_hk_t1)^2
mse_hk_250_t1c=sum(sqrt(diff_hk_250_t1c[,1]+diff_hk_250_t1c[,2]))/dirs

diff_qr_250_t1c=(int_dir_t1-int_dir_qr_t1)^2
mse_qr_250_t1c=sum(sqrt(diff_qr_250_t1c[,1]+diff_qr_250_t1c[,2]))/dirs

################## connected envelope for 0.05 ##################################
## for pressure level 300

dirs=100

int_dir_t1=matrix(0,nrow=dirs,ncol=2)

qln = mqli(cbind(hst1c[,5],vst1c[,5]),prob=0.05,dirs=dirs)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


################ hyperplane kriging #########
probs = 0.05
int_dir_hk_t1=matrix(0,nrow=dirs,ncol=2)

qln70 = mqli(as.matrix(comb70),prob=probs,dirs=dirs); qln100 = mqli(as.matrix(comb100),prob=probs,dirs=dirs); qln200 = mqli(as.matrix(comb200),prob=probs,dirs=dirs); qln250 = mqli(as.matrix(comb250),prob=probs,dirs=dirs); qln300 = mqli(as.matrix(comb300),prob=probs,dirs=dirs);  qln400 = mqli(as.matrix(comb400),prob=probs,dirs=dirs); qln500 = mqli(as.matrix(comb500),prob=probs,dirs=dirs); qln700 = mqli(as.matrix(comb700),prob=probs,dirs=dirs)
dirq70= -qln70[,3]; dirq100= -qln100[,3]; dirq200= -qln200[,3]; dirq250= -qln250[,3]; dirq300= -qln300[,3]; dirq400= -qln400[,3]; dirq500= -qln500[,3]; dirq700= -qln700[,3]

qs=matrix(0,100,length(t[-4]))
for(s in 1:dirs)
{
  qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq250[s],dirq400[s],dirq500[s],dirq700[s])
  
}
##univariate kriging with matern covaraince function
qpred_k=rep(0,dirs)

for(s in 1:dirs){
  pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean = NA)
  df=data.frame(qs[s,])
  pars <- RFfit(pars.model, x=t[-5], dim = 1, data = df)
  
  data.df=data.frame(x=t[-5],v1=qs[s,])
  krig=RFinterpolate(pars,x=data.frame(x=t[5]),data =data.df )
  qpred_k[s]=krig@data$v1
}


qln=cbind(qln70[,1],qln70[,2],-qpred_k)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_hk_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_hk_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


#################### quantile regression #############################
int_dir_qr_t1=matrix(0,nrow=dirs,ncol=2)
probs=0.05

horizontalspeeds=c(comb70[,1],comb100[,1],comb200[,1],comb250[,1],comb400[,1],comb500[,1],comb700[,1])
verticalspeeds=c(comb70[,2],comb100[,2],comb200[,2],comb250[,2], comb400[,2],comb500[,2],comb700[,2])                       
pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(200,sum(time1)),rep(250,sum(time1)),rep(400,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))

xx=horizontalspeeds      
yy=verticalspeeds
ttt =pressure 

days=c(t[5]) #days=c(70,100,200,250,300,400,500,700)
clrs=c('black','black','black')
pchs=c(3,4,8)

qd = qdir(dirs)
df=6
predmat = matrix(0,length(days),dirs)   
xy = cbind(xx,yy)


for (kk in 1:dirs){
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))}
##predicting the tau^th directional quantile

qln=cbind(qd,-predmat[1,])

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_qr_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_qr_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}

diff_hk_300_t1c=(int_dir_t1-int_dir_hk_t1)^2
mse_hk_300_t1c=sum(sqrt(diff_hk_300_t1c[,1]+diff_hk_300_t1c[,2]))/dirs

diff_qr_300_t1c=(int_dir_t1-int_dir_qr_t1)^2
mse_qr_300_t1c=sum(sqrt(diff_qr_300_t1c[,1]+diff_qr_300_t1c[,2]))/dirs

mse_hk_300_t1c
mse_qr_300_t1c

## for pressure level 400

dirs=100

int_dir_t1=matrix(0,nrow=dirs,ncol=2)

qln = mqli(cbind(hst1c[,6],vst1c[,6]),prob=0.05,dirs=dirs)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


################ hyperplane kriging #########
probs = 0.05
int_dir_hk_t1=matrix(0,nrow=dirs,ncol=2)

qln70 = mqli(as.matrix(comb70),prob=probs,dirs=dirs); qln100 = mqli(as.matrix(comb100),prob=probs,dirs=dirs); qln200 = mqli(as.matrix(comb200),prob=probs,dirs=dirs); qln250 = mqli(as.matrix(comb250),prob=probs,dirs=dirs); qln300 = mqli(as.matrix(comb300),prob=probs,dirs=dirs);  qln400 = mqli(as.matrix(comb400),prob=probs,dirs=dirs); qln500 = mqli(as.matrix(comb500),prob=probs,dirs=dirs); qln700 = mqli(as.matrix(comb700),prob=probs,dirs=dirs)
dirq70= -qln70[,3]; dirq100= -qln100[,3]; dirq200= -qln200[,3]; dirq250= -qln250[,3]; dirq300= -qln300[,3]; dirq400= -qln400[,3]; dirq500= -qln500[,3]; dirq700= -qln700[,3]

qs=matrix(0,100,length(t[-4]))
for(s in 1:dirs)
{
  qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq250[s],dirq300[s],dirq500[s],dirq700[s])
  
}
##univariate kriging with matern covaraince function
qpred_k=rep(0,dirs)

for(s in 1:dirs){
  pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean = NA)
  df=data.frame(qs[s,])
  pars <- RFfit(pars.model, x=t[-6], dim = 1, data = df, sub.methods = "self")
  
  data.df=data.frame(x=t[-6],v1=qs[s,])
  krig=RFinterpolate(pars,x=data.frame(x=t[6]),data =data.df )
  qpred_k[s]=krig@data$v1
}


qln=cbind(qln70[,1],qln70[,2],-qpred_k)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_hk_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_hk_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


#################### quantile regression #############################
int_dir_qr_t1=matrix(0,nrow=dirs,ncol=2)
probs=0.05

horizontalspeeds=c(comb70[,1],comb100[,1],comb200[,1],comb250[,1],comb300[,1],comb500[,1],comb700[,1])
verticalspeeds=c(comb70[,2],comb100[,2],comb200[,2],comb250[,2], comb300[,2],comb500[,2],comb700[,2])                       
pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(200,sum(time1)),rep(250,sum(time1)),rep(300,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))

xx=horizontalspeeds      
yy=verticalspeeds
ttt =pressure 

days=c(t[6]) #days=c(70,100,200,250,300,400,500,700)
clrs=c('black','black','black')
pchs=c(3,4,8)

qd = qdir(dirs)
df=6
predmat = matrix(0,length(days),dirs)   
xy = cbind(xx,yy)


for (kk in 1:dirs){
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))}
##predicting the tau^th directional quantile

qln=cbind(qd,-predmat[1,])

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_qr_t1[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_qr_t1[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


diff_hk_400_t1c=(int_dir_t1-int_dir_hk_t1)^2
mse_hk_400_t1c=sum(sqrt(diff_hk_400_t1c[,1]+diff_hk_400_t1c[,2]))/dirs

diff_qr_400_t1c=(int_dir_t1-int_dir_qr_t1)^2
mse_qr_400_t1c=sum(sqrt(diff_qr_400_t1c[,1]+diff_qr_400_t1c[,2]))/dirs

mse_hk_100_t1c;mse_hk_200_t1c;mse_hk_250_t1c;mse_hk_300_t1c;mse_hk_400_t1c
mse_qr_100_t1c;mse_qr_200_t1c;mse_qr_250_t1c;mse_qr_300_t1c;mse_qr_400_t1c


### time period 2
#################################################################################
################## connected envelope for 0.05 ##################################
## for pressure level 100

dirs=100

int_dir_t2=matrix(0,nrow=dirs,ncol=2)

qln = mqli(cbind(hst2c[,2],vst2c[,2]),prob=0.05,dirs=dirs)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


################ hyperplane kriging #########
probs = 0.05
int_dir_hk_t2=matrix(0,nrow=dirs,ncol=2)

qln70 = mqli(as.matrix(comb70_t2),prob=probs,dirs=dirs); qln100 = mqli(as.matrix(comb100_t2),prob=probs,dirs=dirs); qln200 = mqli(as.matrix(comb200_t2),prob=probs,dirs=dirs); qln250 = mqli(as.matrix(comb250_t2),prob=probs,dirs=dirs); qln300 = mqli(as.matrix(comb300_t2),prob=probs,dirs=dirs);  qln400 = mqli(as.matrix(comb400_t2),prob=probs,dirs=dirs); qln500 = mqli(as.matrix(comb500_t2),prob=probs,dirs=dirs); qln700 = mqli(as.matrix(comb700_t2),prob=probs,dirs=dirs)
dirq70= -qln70[,3]; dirq100= -qln100[,3]; dirq200= -qln200[,3]; dirq250= -qln250[,3]; dirq300= -qln300[,3]; dirq400= -qln400[,3]; dirq500= -qln500[,3]; dirq700= -qln700[,3]

qs=matrix(0,100,length(t[-3]))
for(s in 1:dirs)
{
  qs[s,]=c(dirq70[s],dirq200[s],dirq250[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
  
}
##univariate kriging with matern covaraince function
qpred_k=rep(0,dirs)

for(s in 1:dirs){
  pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean = NA)
  df=data.frame(qs[s,])
  pars <- RFfit(pars.model, x=t[-2], dim = 1, data = df, sub.methods ="self")
  
  data.df=data.frame(x=t[-2],v1=qs[s,])
  krig=RFinterpolate(pars,x=data.frame(x=t[2]),data =data.df )
  qpred_k[s]=krig@data$v1
}


qln=cbind(qln70[,1],qln70[,2],-qpred_k)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_hk_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_hk_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


#################### quantile regression #############################
int_dir_qr_t2=matrix(0,nrow=dirs,ncol=2)
probs=0.05

horizontalspeeds=c(comb70_t2[,1],comb200_t2[,1],comb250_t2[,1],comb300_t2[,1],comb400_t2[,1],comb500_t2[,1],comb700_t2[,1])
verticalspeeds=c(comb70_t2[,2],comb200_t2[,2],comb250_t2[,2],comb300_t2[,2], comb400_t2[,2],comb500_t2[,2],comb700_t2[,2])                       
pressure=c(rep(70,sum(time2)), rep(200,sum(time2)), rep(250,sum(time2)),rep(300,sum(time2)),rep(400,sum(time2)),rep(500,sum(time2)), rep(700,sum(time2)))

xx=horizontalspeeds      
yy=verticalspeeds
ttt =pressure 

days=c(t[2]) #days=c(70,100,200,250,300,400,500,700)
clrs=c('black','black','black')
pchs=c(3,4,8)

qd = qdir(dirs)
df=6
predmat = matrix(0,length(days),dirs)   
xy = cbind(xx,yy)


for (kk in 1:dirs){
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))}
##predicting the tau^th directional quantile

qln=cbind(qd,-predmat[1,])

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_qr_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_qr_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


diff_hk_100_t2c=(int_dir_t2-int_dir_hk_t2)^2
mse_hk_100_t2c=sum(sqrt(diff_hk_100_t2c[,1]+diff_hk_100_t2c[,2]))/dirs

diff_qr_100_t2c=(int_dir_t2-int_dir_qr_t2)^2
mse_qr_100_t2c=sum(sqrt(diff_qr_100_t2c[,1]+diff_qr_100_t2c[,2]))/dirs


## for pressure level 200

dirs=100

int_dir_t2=matrix(0,nrow=dirs,ncol=2)

qln = mqli(cbind(hst2c[,3],vst2c[,3]),prob=0.05,dirs=dirs)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


################ hyperplane kriging #########
probs = 0.05
int_dir_hk_t2=matrix(0,nrow=dirs,ncol=2)

qln70 = mqli(as.matrix(comb70_t2),prob=probs,dirs=dirs); qln100 = mqli(as.matrix(comb100_t2),prob=probs,dirs=dirs); qln200 = mqli(as.matrix(comb200_t2),prob=probs,dirs=dirs); qln250 = mqli(as.matrix(comb250_t2),prob=probs,dirs=dirs); qln300 = mqli(as.matrix(comb300_t2),prob=probs,dirs=dirs);  qln400 = mqli(as.matrix(comb400_t2),prob=probs,dirs=dirs); qln500 = mqli(as.matrix(comb500_t2),prob=probs,dirs=dirs); qln700 = mqli(as.matrix(comb700_t2),prob=probs,dirs=dirs)
dirq70= -qln70[,3]; dirq100= -qln100[,3]; dirq200= -qln200[,3]; dirq250= -qln250[,3]; dirq300= -qln300[,3]; dirq400= -qln400[,3]; dirq500= -qln500[,3]; dirq700= -qln700[,3]

qs=matrix(0,100,length(t[-3]))
for(s in 1:dirs)
{
  qs[s,]=c(dirq70[s],dirq100[s],dirq250[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
  
}
##univariate kriging with matern covaraince function
qpred_k=rep(0,dirs)

for(s in 1:dirs){
  pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean = NA)
  df=data.frame(qs[s,])
  pars <- RFfit(pars.model, x=t[-3], dim = 1, data = df)
  
  data.df=data.frame(x=t[-3],v1=qs[s,])
  krig=RFinterpolate(pars,x=data.frame(x=t[3]),data =data.df )
  qpred_k[s]=krig@data$v1
}


qln=cbind(qln70[,1],qln70[,2],-qpred_k)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_hk_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_hk_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


#################### quantile regression #############################
int_dir_qr_t2=matrix(0,nrow=dirs,ncol=2)
probs=0.05

horizontalspeeds=c(comb70_t2[,1],comb100_t2[,1],comb250_t2[,1],comb300_t2[,1],comb400_t2[,1],comb500_t2[,1],comb700_t2[,1])
verticalspeeds=c(comb70_t2[,2],comb100_t2[,2],comb250_t2[,2],comb300_t2[,2], comb400_t2[,2],comb500_t2[,2],comb700_t2[,2])                       
pressure=c(rep(70,sum(time2)), rep(100,sum(time2)), rep(250,sum(time2)),rep(300,sum(time2)),rep(400,sum(time2)),rep(500,sum(time2)), rep(700,sum(time2)))

xx=horizontalspeeds      
yy=verticalspeeds
ttt =pressure 

days=c(t[3]) #days=c(70,100,200,250,300,400,500,700)
clrs=c('black','black','black')
pchs=c(3,4,8)

qd = qdir(dirs)
df=6
predmat = matrix(0,length(days),dirs)   
xy = cbind(xx,yy)


for (kk in 1:dirs){
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))}
##predicting the tau^th directional quantile

qln=cbind(qd,-predmat[1,])

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_qr_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_qr_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


diff_hk_200_t2c=(int_dir_t2-int_dir_hk_t2)^2
mse_hk_200_t2c=sum(sqrt(diff_hk_200_t2c[,1]+diff_hk_200_t2c[,2]))/dirs

diff_qr_200_t2c=(int_dir_t2-int_dir_qr_t2)^2
mse_qr_200_t2c=sum(sqrt(diff_qr_200_t2c[,1]+diff_qr_200_t2c[,2]))/dirs


## for pressure level 250

dirs=100

int_dir_t2=matrix(0,nrow=dirs,ncol=2)

qln = mqli(cbind(hst2c[,4],vst2c[,4]),prob=0.05,dirs=dirs)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


################ hyperplane kriging #########
probs = 0.05
int_dir_hk_t2=matrix(0,nrow=dirs,ncol=2)

qln70 = mqli(as.matrix(comb70),prob=probs,dirs=dirs); qln100 = mqli(as.matrix(comb100),prob=probs,dirs=dirs); qln200 = mqli(as.matrix(comb200),prob=probs,dirs=dirs); qln250 = mqli(as.matrix(comb250),prob=probs,dirs=dirs); qln300 = mqli(as.matrix(comb300),prob=probs,dirs=dirs);  qln400 = mqli(as.matrix(comb400),prob=probs,dirs=dirs); qln500 = mqli(as.matrix(comb500),prob=probs,dirs=dirs); qln700 = mqli(as.matrix(comb700),prob=probs,dirs=dirs)
dirq70= -qln70[,3]; dirq100= -qln100[,3]; dirq200= -qln200[,3]; dirq250= -qln250[,3]; dirq300= -qln300[,3]; dirq400= -qln400[,3]; dirq500= -qln500[,3]; dirq700= -qln700[,3]

qs=matrix(0,100,length(t[-4]))
for(s in 1:dirs)
{
  qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq300[s],dirq400[s],dirq500[s],dirq700[s])
  
}
##univariate kriging with matern covaraince function
qpred_k=rep(0,dirs)

for(s in 1:dirs){
  pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean = NA)
  df=data.frame(qs[s,])
  pars <- RFfit(pars.model, x=t[-4], dim = 1, data = df)
  
  data.df=data.frame(x=t[-4],v1=qs[s,])
  krig=RFinterpolate(pars,x=data.frame(x=t[4]),data =data.df )
  qpred_k[s]=krig@data$v1
}


qln=cbind(qln70[,1],qln70[,2],-qpred_k)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_hk_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_hk_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


#################### quantile regression #############################
int_dir_qr_t2=matrix(0,nrow=dirs,ncol=2)
probs=0.05

horizontalspeeds=c(comb70[,1],comb100[,1],comb200[,1],comb300[,1],comb400[,1],comb500[,1],comb700[,1])
verticalspeeds=c(comb70[,2],comb100[,2],comb200[,2],comb300[,2], comb400[,2],comb500[,2],comb700[,2])                       
pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(200,sum(time1)),rep(300,sum(time1)),rep(400,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))

xx=horizontalspeeds      
yy=verticalspeeds
ttt =pressure 

days=c(t[4]) #days=c(70,100,200,250,300,400,500,700)
clrs=c('black','black','black')
pchs=c(3,4,8)

qd = qdir(dirs)
df=6
predmat = matrix(0,length(days),dirs)   
xy = cbind(xx,yy)


for (kk in 1:dirs){
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))}
##predicting the tau^th directional quantile

qln=cbind(qd,-predmat[1,])

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_qr_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_qr_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}



diff_hk_250_t2c=(int_dir_t2-int_dir_hk_t2)^2
mse_hk_250_t2c=sum(sqrt(diff_hk_250_t2c[,1]+diff_hk_250_t2c[,2]))/dirs

diff_qr_250_t2c=(int_dir_t2-int_dir_qr_t2)^2
mse_qr_250_t2c=sum(sqrt(diff_qr_250_t2c[,1]+diff_qr_250_t2c[,2]))/dirs


## for pressure level 300

dirs=100

int_dir_t2=matrix(0,nrow=dirs,ncol=2)

qln = mqli(cbind(hst2c[,5],vst2c[,5]),prob=0.05,dirs=dirs)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


################ hyperplane kriging #########
probs = 0.05
int_dir_hk_t2=matrix(0,nrow=dirs,ncol=2)

qln70 = mqli(as.matrix(comb70),prob=probs,dirs=dirs); qln100 = mqli(as.matrix(comb100),prob=probs,dirs=dirs); qln200 = mqli(as.matrix(comb200),prob=probs,dirs=dirs); qln250 = mqli(as.matrix(comb250),prob=probs,dirs=dirs); qln300 = mqli(as.matrix(comb300),prob=probs,dirs=dirs);  qln400 = mqli(as.matrix(comb400),prob=probs,dirs=dirs); qln500 = mqli(as.matrix(comb500),prob=probs,dirs=dirs); qln700 = mqli(as.matrix(comb700),prob=probs,dirs=dirs)
dirq70= -qln70[,3]; dirq100= -qln100[,3]; dirq200= -qln200[,3]; dirq250= -qln250[,3]; dirq300= -qln300[,3]; dirq400= -qln400[,3]; dirq500= -qln500[,3]; dirq700= -qln700[,3]

qs=matrix(0,100,length(t[-4]))
for(s in 1:dirs)
{
  qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq250[s],dirq400[s],dirq500[s],dirq700[s])
  
}
##univariate kriging with matern covaraince function
qpred_k=rep(0,dirs)

for(s in 1:dirs){
  pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean = NA)
  df=data.frame(qs[s,])
  pars <- RFfit(pars.model, x=t[-5], dim = 1, data = df)
  
  data.df=data.frame(x=t[-5],v1=qs[s,])
  krig=RFinterpolate(pars,x=data.frame(x=t[5]),data =data.df )
  qpred_k[s]=krig@data$v1
}


qln=cbind(qln70[,1],qln70[,2],-qpred_k)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_hk_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_hk_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


#################### quantile regression #############################
int_dir_qr_t2=matrix(0,nrow=dirs,ncol=2)
probs=0.05

horizontalspeeds=c(comb70[,1],comb100[,1],comb200[,1],comb250[,1],comb400[,1],comb500[,1],comb700[,1])
verticalspeeds=c(comb70[,2],comb100[,2],comb200[,2],comb250[,2], comb400[,2],comb500[,2],comb700[,2])                       
pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(200,sum(time1)),rep(250,sum(time1)),rep(400,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))

xx=horizontalspeeds      
yy=verticalspeeds
ttt =pressure 

days=c(t[5]) #days=c(70,100,200,250,300,400,500,700)
clrs=c('black','black','black')
pchs=c(3,4,8)

qd = qdir(dirs)
df=6
predmat = matrix(0,length(days),dirs)   
xy = cbind(xx,yy)


for (kk in 1:dirs){
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))}
##predicting the tau^th directional quantile

qln=cbind(qd,-predmat[1,])

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_qr_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_qr_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


diff_hk_300_t2c=(int_dir_t2-int_dir_hk_t2)^2
mse_hk_300_t2c=sum(sqrt(diff_hk_300_t2c[,1]+diff_hk_300_t2c[,2]))/dirs

diff_qr_300_t2c=(int_dir_t2-int_dir_qr_t2)^2
mse_qr_300_t2c=sum(sqrt(diff_qr_300_t2c[,1]+diff_qr_300_t2c[,2]))/dirs

## for pressure level 400

dirs=100

int_dir_t2=matrix(0,nrow=dirs,ncol=2)

qln = mqli(cbind(hst2c[,6],vst2c[,6]),prob=0.05,dirs=dirs)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


################ hyperplane kriging #########
probs = 0.05
int_dir_hk_t2=matrix(0,nrow=dirs,ncol=2)

qln70 = mqli(as.matrix(comb70),prob=probs,dirs=dirs); qln100 = mqli(as.matrix(comb100),prob=probs,dirs=dirs); qln200 = mqli(as.matrix(comb200),prob=probs,dirs=dirs); qln250 = mqli(as.matrix(comb250),prob=probs,dirs=dirs); qln300 = mqli(as.matrix(comb300),prob=probs,dirs=dirs);  qln400 = mqli(as.matrix(comb400),prob=probs,dirs=dirs); qln500 = mqli(as.matrix(comb500),prob=probs,dirs=dirs); qln700 = mqli(as.matrix(comb700),prob=probs,dirs=dirs)
dirq70= -qln70[,3]; dirq100= -qln100[,3]; dirq200= -qln200[,3]; dirq250= -qln250[,3]; dirq300= -qln300[,3]; dirq400= -qln400[,3]; dirq500= -qln500[,3]; dirq700= -qln700[,3]

qs=matrix(0,100,length(t[-4]))
for(s in 1:dirs)
{
  qs[s,]=c(dirq70[s],dirq100[s],dirq200[s],dirq250[s],dirq300[s],dirq500[s],dirq700[s])
  
}
##univariate kriging with matern covaraince function
qpred_k=rep(0,dirs)

for(s in 1:dirs){
  pars.model=RMmatern(nu=NA,scale=NA, var=NA)+RMtrend(mean = NA)
  df=data.frame(qs[s,])
  pars <- RFfit(pars.model, x=t[-6], dim = 1, data = df, sub.methods = "self")
  
  data.df=data.frame(x=t[-6],v1=qs[s,])
  krig=RFinterpolate(pars,x=data.frame(x=t[6]),data =data.df )
  qpred_k[s]=krig@data$v1
}


qln=cbind(qln70[,1],qln70[,2],-qpred_k)

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_hk_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_hk_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}


#################### quantile regression #############################
int_dir_qr_t2=matrix(0,nrow=dirs,ncol=2)
probs=0.05

horizontalspeeds=c(comb70[,1],comb100[,1],comb200[,1],comb250[,1],comb300[,1],comb500[,1],comb700[,1])
verticalspeeds=c(comb70[,2],comb100[,2],comb200[,2],comb250[,2], comb300[,2],comb500[,2],comb700[,2])                       
pressure=c(rep(70,sum(time1)), rep(100,sum(time1)), rep(200,sum(time1)),rep(250,sum(time1)),rep(300,sum(time1)),rep(500,sum(time1)), rep(700,sum(time1)))

xx=horizontalspeeds      
yy=verticalspeeds
ttt =pressure 

days=c(t[6]) #days=c(70,100,200,250,300,400,500,700)
clrs=c('black','black','black')
pchs=c(3,4,8)

qd = qdir(dirs)
df=6
predmat = matrix(0,length(days),dirs)   
xy = cbind(xx,yy)


for (kk in 1:dirs){
  predmat[,kk] = predict(rq(xy %*% qd[kk,] ~ bs(ttt,df=df), tau=probs),list(ttt=days))}
##predicting the tau^th directional quantile

qln=cbind(qd,-predmat[1,])

#### calculate active lines
clnid=1
cln=1
stop=FALSE
filn=c(1:dirs)
ss=0
sln=0
while (!stop)
{
  abc=dertln(filn,clnid)
  a=qln[abc[1],]
  b=qln[abc[2],]
  c=qln[abc[3],]
  tt=actln(a,b,c)
  if (ss==0&tt==1) sln=filn[clnid]
  if (tt*ss==1&cln==sln) stop=TRUE
  #print('ss,tt,sln,cln')
  #print(c(ss,tt,sln,cln))
  ss=tt
  m=length(filn)
  if (tt==1)
  {
    if (clnid==m) {cln=filn[1] 
    clnid=1} else 
    {cln=filn[clnid+1] 
    clnid=clnid+1}
  } else { 
    clnidt=clnid
    if (clnid==1) {cln=filn[m] 
    clnid=m-1} else 
    {cln=filn[clnid-1] 
    clnid=clnid-1}
    filn=filn[-clnidt]
  }
} 

actlns=qln[filn,]
m=dim(actlns)[1]
intpts=matrix(rep(0,2*m),nrow=m)
for (i in 1:m) 
{
  if (i==m) j=1 else j=i+1
  intpts[i,]=inln(actlns[i,],actlns[j,])
}

activelines=rep(0,dirs)
activelines[filn]=filn
k=1
### intersection points corresponding to active directions
for(i in 1:100){
  if(activelines[i] != 0) {int_dir_qr_t2[i,]=intpts[k,]; k=k+1 }
  
}
### aprroximation of points corresponding to inactive directions
zero_start=c()
count=numeric()
index=1

i=1 
while(i <= 100){
  count[index]=0
  if(activelines[i]==0){ 
    zero_start[index]=i
    count[index]=count[index]+1
    
    if(i+1<=100){             ### in case the last direction is inactive
      for(j in (i+1):100){
        if(activelines[j]==0){
          count[index]=count[index]+1
          i=i+1
        } else {index= index+1
        break}
      }
    }
  }
  i=i+1
}


for(i in 1:length(zero_start)){
  
  ##taking care of if 1st direction is inactive
  if(zero_start[i]==1){
    x1=intpts[dim(intpts)[1],1]
    y1=intpts[dim(intpts)[1],2]
  }else {
    x1=intpts[which(filn== zero_start[i]-1),1]
    y1=intpts[which(filn==zero_start[i]-1),2]
  }
  
  #taking care of end point
  if(i==length(zero_start) & count[length(count)]!=0){
    x2=intpts[1,1]
    y2=intpts[1,2]
    
  }
  else{
    x2=intpts[which(filn==zero_start[i]+count[i]),1]
    y2=intpts[which(filn==zero_start[i]+count[i]),2]
  }
  
  if(x1!=x2){
    x_inter=runif(count[i],min(x1,x2),max(x1,x2))
    y_inter= y1 + ((y2-y1)/(x2-x1))*(x_inter-x1)
  } else {
    x_inter= rep(x1, count[i])
    y_inter= runif(count[i], min(y1,y2), max(y1,y2))
  }
  int_dir_qr_t2[zero_start[i]:(zero_start[i]+count[i]-1),]=t(rbind(x_inter,y_inter))
}



diff_hk_400_t2c=(int_dir_t2-int_dir_hk_t2)^2
mse_hk_400_t2c=sum(sqrt(diff_hk_400_t2c[,1]+diff_hk_400_t2c[,2]))/dirs

diff_qr_400_t2c=(int_dir_t2-int_dir_qr_t2)^2
mse_qr_400_t2c=sum(sqrt(diff_qr_400_t2c[,1]+diff_qr_400_t2c[,2]))/dirs


mse_hk_100_t2c;mse_hk_200_t2c;mse_hk_250_t2c;mse_hk_300_t2c;mse_hk_400_t2c
mse_qr_100_t2c;mse_qr_200_t2c;mse_qr_250_t2c;mse_qr_300_t2c;mse_qr_400_t2c



