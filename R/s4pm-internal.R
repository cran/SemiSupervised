################################################################################
##
##  Programmer: Mark Vere Culp
##
##  Date:        August 26, 2016
##
##  Description: s4pm internal functionality. These are the internal
##               cross validation and fitting functions for s4pm.
##
################################################################################

s4pm.cv<-function(graph,y,x=NULL,weights,tunes=matrix(c(NA,0.1,0.1,0.0),1),type=0L,control=SemiSupervised.control()){
  n=graph[[4]]
  m=length(y)
  
  if(missing(weights))weights=rep(1.0,n)

  set.seed(control$cv.seed)
  cv.fold=control$cv.fold
  all.folds <- cv.folds(m,folds=cv.fold)
  cv.fold=length(all.folds)
  all.folds <- lapply(1:cv.fold,function(i)sort(all.folds[[i]]))
  if(type==1L){
    if(control$cv.cl){
      cl=c(-1.0,1.0)
      y=cl[as.numeric(as.factor(y))]
      force1<-NULL
      for(j in 1:length(cl)){
        force1<- c(force1,which(y==cl[j])[1])
      }
      force1=sort(force1)
      all.folds<-	lapply(1:cv.fold,function(i)unique(sort(c(force1,all.folds[[i]]))))
    }
  }
  omits<-lapply(1:cv.fold,function(i)setdiff(1:m,all.folds[[i]]))
  
  l.thresh=as.integer(floor(control$l.thresh))
  l.eps=as.numeric(control$l.eps)
  if(is.vector(tunes)){
    tunes=matrix(tunes,1)
  }
  
  cv.str<-.Call("cv_s4pm_fit",all.folds,omits,graph,x,as.numeric(y),weights,
  tunes,l.thresh,l.eps,type,PACKAGE="SemiSupervised")
  tunes=cbind(tunes,sqrt(cv.str[[1]]/m),cv.str[[2]],cv.str[[3]])
  tunes=as.data.frame(tunes)
  names(tunes)=c("hs","lam1","lam2","gam","cv.errs","cv.convs","cv.internal")
  tunes=as.matrix(tunes)
  return(list(tunes=tunes,folds=all.folds))
}

s4pm.fit<-function(graph,y,x=NULL,weights,tunes=c(NA,0.1,0.1,0.0),type=0L,est.only=FALSE, control=SemiSupervised.control()){
  xdat=!is.null(x)
  n=graph[[4]]
  m=length(y)
  
  p=0
  if(xdat){
    p=dim(x)[2]
  }
  
  if(missing(weights))weights=rep(1.0,n)
  
  l.thresh=as.integer(control$l.thresh)
  l.eps=as.numeric(control$l.eps)
  a1=as.integer(est.only)
  if(a1==1L){
    a1=2L
  }
  
  eta=rep(0.0,n)
  if(xdat){
      fit<-.Call("s4pm_fit",graph,x,y,as.numeric(weights),m,as.numeric(tunes),a1,type,
      l.thresh,l.eps,eta,PACKAGE="SemiSupervised")
  }else{
    fit<-.Call("s4pm_fit",graph,NULL,y,as.numeric(weights),m,as.numeric(tunes),a1,type,
    l.thresh,l.eps,eta,PACKAGE="SemiSupervised")
  }
 
  if(est.only){
    return(eta)
  }
  fitinfo<-list()
  fitinfo[[1]]=eta
  fitinfo[[2]]=matrix(fit[[1]],n,1+as.numeric(xdat))
  fitinfo[[3]]=fit[[3]]
  fitinfo[[4]]=fit[[9]]
  fitinfo[[5]]=fit[[2]]
  fitinfo[[6]]=fit[[6]]
  fitinfo[[7]]=fit[[7]]
  fitinfo[[8]]=fit[[8]]
  fitinfo[[9]]=fit[[5]]
  fitinfo[[10]]=fit[[4]]
  fitinfo[[11]]=mean(eta[1:m])
  fitinfo[[12]]=xdat
  names(fitinfo)<-c("fit","fhat","beta","x.vars","weights","edf","dims","resp.info",
                    "error.code","conv.info","intercept","xdat")
  
  
  xvars=fitinfo$x.vars
  a1<-which(as.numeric(xvars)>0)
  if(length(a1)>0){
    bet=rep(0.0,length(xvars))
    bet[a1]=fitinfo$beta
    fitinfo$beta=bet
    dimnames(fitinfo[[2]])[[2]]=c("f1","f2")
  }
  return(fitinfo)
}

