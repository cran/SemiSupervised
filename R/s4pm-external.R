################################################################################
##
##  Programmer: Mark Vere Culp
##
##  Date:        August 26, 2016
##
##  Description: s4pm functionality. These functions define the s4pm
##               object. All functions in this file are exported.
##               These are all documented in the R help files.
##
################################################################################


## Definition of s4pm methods. These are all documented.

setGeneric("s4pm", function(x, ...) standardGeneric("s4pm"))
setMethod("s4pm",signature(x="formula"),function(x,data,metric= c("cosine","euclidean"),...,
est.only=FALSE,control=SemiSupervised.control()){
  cl <- match.call()
  cl[[1]]<-as.name("s4pm")
  
  env.parent=parent.frame()
  mf <- match.call(expand.dots = FALSE)
  
  m <- match(c("x", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  exg<-extract.xg(x,env.parent)
  control$dissimilar=exg[[4]]
  if(exg[[1]]<3L){
    dat=NULL
    if(!is.null(exg[[2]])){
      mt<-terms(exg[[2]],data=data)
      attr(mt,"intercept")=0
      mf$formula<-mt
      mf$na.action=na.pass
      exg[[2]]=list(terms=mt,frame=mf)
      mf.x<-eval(mf,env.parent)
      dat=model.matrix(mt,mf.x)
    }
    mt<-terms(exg[[3]])
    attr(mt,"intercept")=0
    mf$formula<-mt
    mf$na.action=na.pass
    exg[[3]]=list(terms=mt,frame=mf)
    mf.g<-eval(mf,env.parent)
    metric=NA
    y=model.response(mf.g,"any")
    control$U.as.anchor=FALSE
  }else{
    mt<-terms(x,data=data)
    attr(mt,"intercept")=0
    mf$formula<-mt
    mf$na.action=na.pass
    exg[[2]]=list(terms=mt,frame=mf)
    mf1<-eval(mf,env.parent)
    y=model.response(mf1,"any")
    dat=model.matrix(mt,mf1)
    if(missing(metric))metric="cosine"
    control$dissimilar=TRUE
  }
  L=as.vector(which(!is.na(y)))
  U=as.vector(which(is.na(y)))
  
  n=length(y)
  u=n-length(L)
  if(u<control$U.as.anchor.thresh){
    control$U.as.anchor=FALSE
  }
  if(length(U)==0){
    control$U.as.anchor=FALSE
  }
  if(!is.null(dat)){
    dat=x.scaleL(dat,L)
    x.scaling=list(center=attr(dat,"center"),scale=attr(dat,"scale"))
  }
  if(is.null(control$k)){
    control$k=6L
  }
  
  
  if(!control$U.as.anchor){
    if(exg[[1]]<3L){
      graph=as.matrix(model.matrix(mt,mf.g))
    }else{
      graph=as.matrix(dG(dat,k=control$k,nok=control$nok,metric=metric))
    }
    obj=s4pm.default(dat,y,graph,...,est.only=est.only,control=control)
    if(est.only){
      return(obj)
    }
    if(exg[[1]]>2L){
      slot(obj,".fitinfo")$as.default=FALSE
    }
  }else{
    anchor=getAnchor(dat,list(anchor.seed=control$anchor.seed,k=control$U.as.anchor.thresh,
    iter.max=control$iter.max))
    ndat=rbind(dat[L,,drop=FALSE],anchor)
    tgraph=as.matrix(dG(ndat,k=control$k,nok=control$nok,metric=metric))
    obj=s4pm.default(ndat,y[L],tgraph,...,est.only=FALSE,control=control)
  }
  
  slot(obj,".call")<-cl
  exg[[5]]=env.parent
  slot(obj,".terminfo")<-exg
  slot(obj,".fitinfo")$metric=metric
  if(control$U.as.anchor){
    graph=knnGraph(dat,ndat,k=control$k,nok=control$nok,metric=metric)
    dat=as.data.frame(dat)
    names(dat)<-attr(exg[[2]]$terms,"term.labels")
    p=dim(dat)[2]
    attr(dat,"scaled:center")=rep(0.0,p)
    attr(dat,"scaled:scale")=rep(1.0,p)
    
    fhat<-predict(obj,xnew=dat,gnew=graph,type="terms")
    attr(dat,"scaled:center")=x.scaling$center
    attr(dat,"scaled:scale")=x.scaling$scale
    slot(obj,"fit")<-fhat[,1]
    slot(obj,".fitinfo")$fhat=fhat[,-1]
    slot(obj,".fitinfo")$dims[1]=n
    slot(obj,".respinfo")$L=L
    slot(obj,".respinfo")$U=U
    slot(obj,".fitinfo")$as.default=FALSE
  }
  SemiSupervised::gmatrix(obj)<-graph
  SemiSupervised::xmatrix(obj)<-dat
  return(obj)
})
setMethod("s4pm",signature(x="NULL"),function(x,y,graph,...){
  obj<-s4pm.default(NULL,y,graph,...)
  SemiSupervised::gmatrix(obj)=graph
  SemiSupervised::xmatrix(obj)=NULL
  obj
})

setMethod("s4pm",signature(x="matrix"),function(x,y,graph,...){
  m=length(y)
  n=dim(x)[1]
  if(m<n){
    L=1:m
  }else{
    L=as.vector(which(!is.na(y)))
  }
  x=x.scaleL(x,L)
  obj<-s4pm.default(x,y,graph,...)
  SemiSupervised::gmatrix(obj)=graph
  SemiSupervised::xmatrix(obj)=x
  obj
})
setMethod("s4pm",signature(x="vector"),function(x,y,graph,...){
  x=t(t(x))
  m=length(y)
  n=dim(x)[1]
  if(m<n){
    L=1:m
  }else{
    L=as.vector(which(!is.na(y)))
  }
  x=x.scaleL(x,L)
  obj<-s4pm.default(x,y,graph,...)
  SemiSupervised::xmatrix(obj)=x
  SemiSupervised::gmatrix(obj)=graph
  obj
})
setMethod("s4pm",signature(x="data.frame"),function(x,y,graph,...){
  nms=names(x)
  m=length(y)
  n=dim(x)[1]
  if(m<n){
    L=1:m
  }else{
    L=as.vector(which(!is.na(y)))
  }
  x=x.scaleL(x,L)
  obj=s4pm.default(x,y,graph,...)
  obj@.fitinfo$colnames=nms
  SemiSupervised::xmatrix(obj)=x
  SemiSupervised::gmatrix(obj)=graph
  obj
})

## s4pm default function. This provides the main functionality for the s4pm.
## This function is exported but is not expected to be called directly.


"s4pm.default"<-function(x,y,graph,weights,hs,lams,gams,type=c("r","c"),est.only=FALSE,control=SemiSupervised.control()){
  cl <- match.call()
  cl[[1]]<-as.name("s4pm")
  
  lev<-ord<-L<-U<-NULL
  
  
  sanity.init(method="s4pm")
  create.tunes(method="s4pm")
  
  obj<-run.model(method="s4pm")
  obj@.fitinfo$fit[ord]=obj@.fitinfo$fit
  slot(obj,"fit")<-obj@.fitinfo$fit
  if(est.only){
    return(obj@fit)
  }
  obj@.fitinfo$fhat[ord,]<-obj@.fitinfo$fhat
  obj@.fitinfo$weights[ord]<-obj@.fitinfo$weights
  muy=obj@.fitinfo$resp.info[1]
  ssy=obj@.fitinfo$resp.info[2]
  slot(obj,".call")<-cl
  slot(obj,".respinfo")<-list(L=as.vector(L),U=as.vector(U),lev=lev,muy=muy,ssy=ssy)
  obj@.fitinfo$resp.info<-NULL
  slot(obj,"hparm")<-obj@.fitinfo$h
  slot(obj,"lparm")<-obj@.fitinfo$lam
  slot(obj,"gparm")<-obj@.fitinfo$gam
  xmatrix(obj)=NULL
  obj
}

## Generic functions for class s4pm are defined next.

setMethod("show", signature(object = "s4pm"),function(object){
  if(class(object)[1]!="s4pm"){
    stop("Error: jtharm is not of type s4pm")
  }
  cl <- object@.call
  ind=grep("s4pm",cl)
  if(length(ind)==0){
    cat("Empty s4pm object\n")
    return()
  }
  type=object@type
  n=object@.fitinfo$dims[1]
  m=object@.fitinfo$dims[2]
  cv.err=object@.cv_str$opt.row[5]
  edf=measures(object)[3]
  df=m-edf
  dis=object@.control$dissimilar
  if(df<0)df=0.0
  
  cat("S4PM Fit with (n,|L|)=(",n,",",m,") or ",round(m/n*100,0),"% labeled\n")
  cat("\nPerformance Estimates:\nk-CV: ",round(cv.err,3)," GCV: ",round(measures(object)[2],3)," DF: ",round(df,3))
  cat("\n\nFit Estimates:\n")
  if(dis){
    cat("Graph Kernel h: ",round(hparm(object),3)," ")
  }
  cat("Lagrangians: ",round(lparm(object)[1],3)[1]," ",round(gparm(object),3)," ")
  if(object@.fitinfo$xdat){
    cat("Safe-Lagrangian",round(lparm(object)[2],3)," ")
  }
  cat("\n\n")
  invisible(object)
})

setMethod("predict", signature(object = "s4pm"), function(object,xnew,gnew,type=c("vector","response","prob","terms"),pow=1,...){
  if(class(object)!="s4pm"){
    stop("Error: object must be of type s4pm")
  }
  if(missing(type)){
    type="response"
  }
  default.mode=object@.fitinfo$as.default
  xdat=object@.fitinfo$xdat
  
  form.struct<-object@.terminfo
  if(!is.null(form.struct)){
    if(!is.null(form.struct[[2]])){
      tt=form.struct[[2]]$terms
      Terms <- delete.response(tt)
      m <- model.frame(Terms, xnew, na.action = na.pass)
      xnew <- model.matrix(Terms, m)
    }
  }

  if(!missing(xnew)){
    xnew=x.scaleL(xnew,sanity=TRUE,sanity.only=TRUE)
    if(!is.null(SemiSupervised::xmatrix(object))){
      xnew=scale(xnew,center=attr(SemiSupervised::xmatrix(object),"scaled:center"),scale=attr(SemiSupervised::xmatrix(object),"scaled:scale"))
    }
    p=dim(xnew)[2]
    
    if(p!=length(object@.fitinfo$x.vars)){
      stop("Error: supplied xnew does not have correct number of variables")
    }
  }else{
    if(xdat)stop("Error: Need xnew since model was fit in safe mode")
  }
  if(missing(gnew)){
    if(!default.mode){
      if(missing(xnew))stop("Error: Must input xnew or gnew")
      gnew=knnGraph(xnew,object@xmatrix,object@.control$k,object@.control$nok,object@.fitinfo$metric)
    }else{
      stop("Error: missing gnew")
    }
  }
  gnew=x.scaleL(gnew,sanity.only=TRUE)
  n=dim(gnew)[1]
  np=dim(gnew)[2]
  if(dim(object)[1]!=np){
    stop("Error: The dims gnew do not match: Expect: (",dim(object)[1],") Inputted: (",np,")")
  }
  if(missing(xnew)& xdat){
    stop("Error: object xnew must be supplied whenever s4pm was fit with x data")
  }
  if(!missing(xnew)){
    if(n!=dim(xnew)[1]){
      stop("Error: xnew and gnew must have the same number of rows")
    }
  }

  ctype=object@type
  yvec=object@fit
  if(object@.control$dissimilar) gnew=exp(-gnew/object@hparm)
  if(xdat){
    f2=rep(0,dim(object)[1])
    f2=object@.fitinfo$fhat[,"f2"]
    yvec=(yvec-f2)
    bet=object@.fitinfo$beta
  }else{
    xnew=NULL
    bet=NULL
  }
  term=FALSE
  if(type=="terms"){
    term=TRUE
  }

  fit=inter.predict(yvec,gnew,xnew,bet,pow,term)
  if(term){
    return(fit)
  }
  if(ctype=="r" & type=="response"){
    type="vector"
  }

  if(type=="prob" & ctype=="r"){
    type="vector"
  }
  if(type=="vector"|ctype=="r"){
    return(fit)
  }
  if(type=="prob"){
    return(exp(fit)/(1+exp(fit)))
  }
  ty<-as.factor(object@.respinfo$lev[as.numeric(fit>0.0)+1])
  return(ty)
})
