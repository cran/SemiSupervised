################################################################################
##
##  Programmer: Mark Vere Culp
##
##  Date:        August 26, 2016
##
##  Description: Anchor Graph internal functionality. These are the internal
##               cross validation and fitting functions for agraph.
##
################################################################################
agraph.cv<-function(z,x=NULL,W,y,weights,tunes=matrix(c(1.0,0.1,0.1),1),type=0L,control=SemiSupervised.control()){
  n=dim(z)[1]
  m=length(y)
  if(missing(weights))weights=rep(1.0,n)
  
  set.seed(control$cv.seed)
  all.folds <- cv.folds(m,folds=control$cv.fold)
  cv.fold=length(all.folds)
  all.folds <- lapply(1:cv.fold,function(i)sort(all.folds[[i]]))
  if(nlevels(as.factor(y))<3L){
    if(control$cv.cl){
      cl=c(-1.0,1.0)
      y=cl[as.numeric(as.factor(y))]
      force1<-NULL
      for(j in 1:length(control$cl)){
        force1<- c(force1,which(y==cl[j])[1])
      }
      force1=sort(force1)
      all.folds<-	lapply(1:cv.fold,function(i)unique(sort(c(force1,all.folds[[i]]))))
    }
  }
  omits<-lapply(1:cv.fold,function(i)setdiff(1:m,all.folds[[i]]))
  l.thresh=as.integer(control$l.thresh)
  l.eps=as.numeric(control$l.eps)
  type=as.integer(type)
  
  if(is.vector(tunes)){
    tunes=matrix(tunes,1)
  }

  
  if(is.null(x)){
    cv.str<-.Call("AREG_CV",all.folds,omits,z,NULL,W,as.numeric(y),weights,tunes,type,l.thresh,l.eps,as.numeric(control$stability),PACKAGE="SemiSupervised")
  }else{
    cv.str<-.Call("AREG_CV",all.folds,omits,z,x,W,as.numeric(y),weights,tunes,type,l.thresh,l.eps,as.numeric(control$stability),PACKAGE="SemiSupervised")
  }
  tunes=cbind(tunes,sqrt(cv.str[[1]]/m),cv.str[[2]],cv.str[[3]])
  tunes=as.data.frame(tunes)
  names(tunes)=c("lam1","lam2","gam1","cv.errs","cv.convs","cv.internal")
  tunes=as.matrix(tunes)
  return(list(tunes=tunes,folds=all.folds))
}

agraph.fit<-function(z,x=NULL,W,y,weights,tunes=c(1.0,0.1,0.1),type=0L,est.only=FALSE,control=SemiSupervised.control()){
  n=dim(z)[1]
  m=length(y)
  L=1:m
  if(missing(weights))weights=rep(1.0,n)
  
  l.thresh=as.integer(control$l.thresh)
  l.eps=as.numeric(control$l.eps)
  adjust=as.numeric(control$adjust)
  
  eta=rep(0.0,n)
  fit<-.Call("ARIDGE",z,x,W,as.numeric(y),weights,m,tunes,type,l.thresh,l.eps,
  as.integer(est.only),as.numeric(control$stability),as.numeric(eta),PACKAGE="SemiSupervised")
  
  if(est.only){
    return(eta)
  }
  fitinfo<-list()
  fitinfo[[1]]=eta
  fitinfo[[2]]=fit[[1]]
  fitinfo[[3]]<-fit[[2]]
  fitinfo[[4]]<-fit[[7]]
  fitinfo[[5]]<-fit[[4]]
  fitinfo[[6]]<-fit[[5]]
  fitinfo[[7]]<-fit[[6]]
  fitinfo[[8]]<-fit[[3]]
  fitinfo[[9]]<-mean(eta[L])
  fitinfo[[10]]<-FALSE
  
  names(fitinfo)=c("fit","coef","x.vars","weights","edf","dims","resp.info","error.code","intercept","xdat")
  if(fitinfo$edf<0.001){
    fitinfo$edf=0.001
  }
  if(!is.null(x)){
    fitinfo$xdat<-TRUE
  }
  return(fitinfo)
}
