################################################################################
##
##  Programmer: Mark Vere Culp
##
##  Date:        August 26, 2016
##
##  Description: jtharm internal functionality. These are the internal
##               cross validation and fitting functions for jtharm.
##
################################################################################
jtharm.cv<-function(graph,y,weights,tunes=c(NA,1.0,0.0),type=0L,control=SemiSupervised.control()){
  n=graph[[4]]
  m=length(y)
  if(missing(weights))weights=rep(1.0,n)
  
  set.seed(control$cv.seed)
  cv.fold=control$cv.fold
  all.folds <- cv.folds(m,folds=cv.fold)
  cv.fold=length(all.folds)
  all.folds <- lapply(1:cv.fold,function(i)sort(all.folds[[i]]))
  if(nlevels(as.factor(y))<3L){
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
  if(is.vector(tunes)){
    tunes=matrix(tunes,1)
  }

  cv.str<-.Call("cv_jtharm_fit",all.folds,omits,graph,as.numeric(y),weights,tunes,as.numeric(control$stability),type,PACKAGE="SemiSupervised")
  tunes=cbind(tunes,sqrt(cv.str[[1]]/m),cv.str[[2]])
  tunes=as.data.frame(tunes)
  names(tunes)=c("hs","lam","gam","cv.errs","cv.convs")
  tunes=as.matrix(tunes)
  return(list(tunes=tunes,folds=all.folds))

}

jtharm.fit<-function(graph,y,weights,tunes=c(NA,1.0,0.0),type=0L,est.only=FALSE,stab=0.0,control=SemiSupervised.control()){
  n=graph[[4]]
  m=length(y)
  tunes=as.numeric(tunes)
  if(missing(weights))weights=rep(1.0,n)
  eta=rep(0.0,n);
  if(is.null(control$stability)){
    control$stability=0.0;
  }
  
  
  fit<-.Call("jt_harm_fit", graph, y, as.numeric(weights), as.numeric(tunes[2]), as.numeric(tunes[3]),as.numeric(control$stability), as.integer(est.only), as.numeric(eta),PACKAGE="SemiSupervised")

  if(est.only){
    return(eta)
  }
  fitinfo<-list()
  fitinfo[[1]]<-eta
  fitinfo[[2]]<-fit[[1]]
  fitinfo[[3]]<-fit[[2]]
  fitinfo[[4]]<-fit[[3]]
  fitinfo[[5]]<-fit[[4]]
  fitinfo[[6]]<-fit[[5]]
  names(fitinfo)<-c("fit","fitted_response","edf","dims","error.code","resp.info")
  return(fitinfo)
}
