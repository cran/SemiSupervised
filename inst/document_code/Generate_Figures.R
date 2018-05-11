################################################################################
##      Programmer: Mark Vere Culp
##      Date:       1-13-2017
##      Purpose:    This file prints all figures in the paper.
################################################################################


#####
## Install libraries
####

libloader<-function(txt){
  g<-try(library(package=txt,character.only=TRUE,quietly=TRUE),silent=TRUE)
  if(class(g)=="try-error"){
    stop(paste("Must install package: ",txt," in order to proceed.",sep=""))
  }
}
g<-lapply(c("spa","kernlab","MASS"),libloader)
if(!file.exists("source/data_reader.R")){stop("Must be in the main code directory (refer to README.txt)")}
source("source/data_reader.R")
#####
## Set up preliminaries
#####
paper=FALSE## Set flag for paper

if(!paper)pdf("all_figures.pdf")
####
## Plot Figure Two Moons (Figure 1a in the paper)
###
set.seed(100)
dat<-spa.sim(n=1000)
L<-as.vector(which(!is.na(dat[,1])))
x<-x.scaleL(dat[,-1],L)
Dij<-euclidDist(x)
g1<-s4pm(y~dG(Dij),data=dat,control=s4pm.control(adjust=0.0),gams=0.0)
g<-ksvm(y~.,data=data.frame(y=as.factor(dat[,1]),x)[L,],prob.model=TRUE)

lgrid<-500
zhat=matrix(0,lgrid,lgrid)
zsup=matrix(0,lgrid,lgrid)
xgrid1<-seq(min(x[,1]),max(x[,1]),length=lgrid)
xgrid2<-seq(min(x[,2]),max(x[,2]),length=lgrid)
gr<-gray(c(.3,.7))[c(rep(1,500),rep(2,500))]
t1<-proc.time()
for(i in 1:lgrid){
  x1<-data.frame(X1=xgrid1[i],X2=xgrid2)
  ngraph<-kgraph.predict(euclidDist(x1,x))
  zhat[,i]<-predict(g1,gnew=ngraph,type="vector")
  zsup[,i]<-predict(g,newdata=x1,type="prob")[,2]
}
cat("time=",(proc.time()-t1)/60,"\n")

if(paper)pdf("fig1a.pdf")
plot(x,col=gr,pch=16,cex=0.5,xlab="",ylab="",main="")
points(x,col=gr,pch=16,cex=0.5)
points(x[L,],col=gr[L],pch=16,cex=2)
points(x[L,],col=1,pch=1,cex=2)
contour(xgrid1,xgrid2,t(zhat),levels=0.5,method="edge",add=TRUE,lwd=2)
contour(xgrid1,xgrid2,t(zsup),levels=0.5,method="edge",add=TRUE,lwd=2,lty=2)
if(!paper){
 legend("topright",c("Semi-Supervsied","Supervsied SVM"),lty=c(1,2),lwd=2,bg="white")
 mtext("x1",side=1,padj=3.5)
 mtext("x2",side=2,padj=-3.5)
}
if(paper)dev.off()

####
## Plot Figure Extrapolation 1-D (Figure 1b in the paper)
####

extrapSim1d<-function(l=31,u=69,mux=c(0,0),rhox=-0.9,varx1=1,varx2=1,muz=c(1,0),rhoz=0,varz1=0.5,varz2=0.1){
  n=l+u
  xvar=cbind(c(varx1,rhox),c(rhox,varx2))
  zvar=cbind(c(varz1,rhoz),c(rhoz,varz2))
  x=scale(mvrnorm(n,mux,xvar),scale=F)
  z=mvrnorm(n,muz,zvar)
  x=rbind(x[1:l,],z[(l+1):n,])
  
  mls<-apply(x[1:l,],2,mean)
  vls<-apply(x[1:l,],2,var)
  x=scale(x,center=mls,scale=sqrt(vls))
  xL=x[1:l,]/sqrt(l-1)
  xU=x[(l+1):n,]/sqrt(n-l)
  list(x=rbind(xL,xU),xL=xL,xU=xU,L=1:l,U=(l+1):n,mls=mls,vls=vls)
}
set.seed(100)
dat<-extrapSim1d(muz=c(0,3),varz1=0.05)
if(paper)pdf("fig1b.pdf")
plot(dat$x,col=c("white","gray")[c(rep(1,31),rep(2,69))],pch=16,cex=1,ylab="",xlab="",xlim=c(-0.5,0.5),ylim=c(-0.5,1.0))
points(dat$x,cex=1)
if(!paper){
  legend("topright",c("Labeled","Unlabeled"),fill=c("white","gray"))
  mtext("x1",side=1,padj=3.5)
  mtext("x2",side=2,padj=-3.5)
  b1<-eigen(t(dat$xU)%*%dat$xU%*%solve(t(dat$xL)%*%dat$xL))
  vs1=-b1$vec[,1]*0.9
  vs2=-b1$vec[,2]*0.5
  arrows(0,0,vs1[1],vs1[2],lwd=2)
  arrows(0,0,vs2[1],vs2[2],lwd=2)
  text(0.25,-0.15,expression(nu[1]),cex=1.5)
  text(-0.03,0.2,expression(nu[2]),cex=1.5)
}
if(paper)dev.off()

####
## Plot Timings (Figure 2a,b in the paper)
####
x<-read.csv("csvs/perf_results.csv",T)
times<-sapply(c(3,9:10),function(i)tapply(x[,i],x[,1],mean,na.rm=TRUE))
times<-times[order(times[,1]),]


for(i in 1:2){
  if(paper)pdf(paste("fig_time",i,".pdf",sep=""))
  f<-function(x)log(x)
  mn<-range(f(times[,-1]))
  plot(0,1,type="n",xlim=range(f(times[,1])),ylim=c(mn[1],mn[2]),xlab="",ylab="",axes=FALSE,frame=FALSE)
  ivec=rep(FALSE,10)
  ind=c(1,5:10)
  ivec[ind]=TRUE
  points(f(times[ivec,1]),f(times[ivec,i+1]),type="b",pch="",cex=1.5,lwd=2)
  text(f(times[ivec,1]),f(times[ivec,i+1]),paste(ind),cex=1.5)
  points(f(times[!ivec,1]),f(times[!ivec,i+1]),type="b",pch="", col="red",lty=2,cex=1.5,lwd=2)
  text(f(times[!ivec,1]),f(times[!ivec,i+1]),paste(which(!ivec)),cex=1.5)
  axis(1,lwd=2,cex.axis=1.25)
  axis(2,lwd=2,cex.axis=1.25)
  if(!paper){
    title(c("S4PM","AGRAPH")[i])
    mtext("log(n)",side=1,padj=3.5,cex=1.25)
    mtext("log(Time)",side=2,padj=-3.5,cex=1.25)
    legend("topleft",c("Regression","Classification"),lwd=2,lty=c(1,2),col=c(1,2))
  }
  box(lwd=2)
  if(paper)dev.off()
}

if(!paper)dev.off()

q(save="no")
