################################################################################
##
##  Programmer: Mark Vere Culp
##
##  Date:        August 26, 2016
##
##  Description: Anchor Graph functionality. These functions define the anchor
##               graph object. All functions in this file are exported.
##               These are all documented in the R help files.
##
################################################################################


## Definition of agraph methods. These are all documented.

setGeneric("agraph", function(x, ...) standardGeneric("agraph"))
setMethod("agraph",signature(x="formula"),function(x,data,metric= c("cosine","euclidean"),...,
est.only=FALSE,control=SemiSupervised.control()){
  cl <- match.call()
  cl[[1]]<-as.name("agraph")
  
  env.parent=parent.frame()
  mf <- match.call(expand.dots = FALSE)
  
  m <- match(c("x", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  exg<-extract.a(x,env.parent)
  if(exg[[1]]<2L){
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
    if(!is.integer(control$LAE.thresh)){
      control$LAE.thresh=as.integer(ceiling(control$LAE.thresh))
    }
    if(missing(metric)){
      metric="cosine"
    }
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
    x.scaling=list(center=attr(dat,"scaled:center"),scale=attr(dat,"scaled:scale"))
  }
  
  if(exg[[1]]<2L){
    graph=eval(exg[[5]],env.parent)
  }else{
    graph<-AnchorGraph(dat,metric=metric,control=control)
  }
  
  if(!control$U.as.anchor){
    obj=agraph.default(graph,dat,y,...,est.only=est.only,control=control)
    if(est.only){
      return(obj)
    }
    if(exg[[1]]>1L){
      slot(obj,".fitinfo")$as.default=FALSE
    }
  }else{
    tgraph=graph
    tgraph$Z=rbind(graph$Z[L,,drop=FALSE],fitLAE(graph$anchor,graph$anchor,metric=metric,control=control))
    tgraph$rL=graph$rL
    ndat=rbind(dat[L,,drop=FALSE],graph$anchor)
    obj=agraph.default(tgraph,ndat,y[L],...,est.only=FALSE,control=control)
    slot(obj,".fitinfo")$as.default=FALSE
  }
  slot(obj,".call")<-cl
  exg[[6]]=env.parent
  slot(obj,".terminfo")<-exg
  
  if(control$U.as.anchor){
    slot(obj,".fitinfo")$dims[1]=n
    slot(obj,".respinfo")$L=L
    slot(obj,".respinfo")$U=U
    dat=as.data.frame(dat)
    names(dat)<-attr(exg[[2]]$terms,"term.labels")
    slot(obj,"fit")<-predict(obj,xnew=dat,gnew=graph,type="vector")
    attr(dat,"scaled:center")=x.scaling$center
    attr(dat,"scaled:scale")=x.scaling$scale
    
  }
  slot(obj,".fitinfo")$metric=graph$metric
  if(est.only){
    return(slot(obj,"fit"))
  }
  slot(obj,"gmatrix")<-graph
  slot(obj,"xmatrix")<-dat
  return(obj)
})



setMethod("agraph",signature(x="matrix"),function(x,y,...,metric="cosine",est.only=FALSE,control=SemiSupervised.control()){
  cl <- match.call()
  cl[[1]]<-as.name("agraph")
  
  if(is.null(y)|missing(y))stop(paste("Error in y:  must not be NULL or missing"))
  if(missing(y))stop("Error in y:  a response (either continuous or binary) must be provided")
  
  adims=dim(x)
  m=length(y)
  
  if(m==adims[1]){
    L=which(!is.na(y))
    U=which(is.na(y))
    m=length(L)
  }else{
    if(m>adims[1]){
      stop("Error: |y|> dim(graph)[1]")
    }
    L=1:m
    U=(m+1):adims[1]
  }
  
  if(adims[1]<control$U.as.anchor.thresh){
    control$U.as.anchor=FALSE
  }
  if(length(U)==0){
    control$U.as.anchor=FALSE
  }
  x=x.scaleL(x,L)
  
  graph = AnchorGraph(x,metric=metric,control=control)
  
  if(!control$U.as.anchor){
    obj=agraph.default(graph,x,y,...,est.only=est.only,control=control)
    if(est.only){
      return(obj)
    }
  }else{
    tgraph=graph
    tgraph$Z=rbind(graph$Z[L,,drop=FALSE],fitLAE(graph$anchor,graph$anchor,metric=metric,control=control))
    ndat=rbind(x[L,,drop=FALSE],graph$anchor)
    obj=agraph.default(tgraph,ndat,y[L],...,est.only=FALSE,control=control)
  }
  slot(obj,".call")<-cl
  
  obj@.fitinfo$metric=graph$metric
  if(control$U.as.anchor){
    slot(obj,".fitinfo")$dims[1]=adims[1]
    slot(obj,".respinfo")$L=L
    slot(obj,".respinfo")$U=U
    x=as.data.frame(x)
    slot(obj,"fit")<-predict(obj,xnew=x,gnew=graph,type="vector")
  }
  slot(obj,"gmatrix")<-graph
  slot(obj,"xmatrix")<-x
  
  slot(obj,".fitinfo")$as.default=FALSE
  
  if(est.only){
    return(slot(obj,"fit"))
  }
  obj
})
setMethod("agraph",signature(x="data.frame"),function(x,...){
  nms<-names(x)
  x=model.matrix(~.-1,x)
  obj<-agraph(x,...)
  obj@.fitinfo$colnames=nms
  obj
})
setMethod("agraph",signature(x="vector"),function(x,...){
  obj<-agraph(t(t(x)),...)
  obj
})
setMethod("agraph",signature(x="anchor"),function(x,...){
  obj<-agraph.default(graph=x,...)
  gmatrix(obj)=x
  obj
})

## Anchor graph default function. This provides the main functionality for the anchor graph.
## This function is exported but is not expected to be called directly.

agraph.default<-function(graph,x,y,weights,lams,gams,type=c("r","c"),est.only=FALSE,control=SemiSupervised.control()){
  cl <- match.call()
  cl[[1]]<-as.name("agraph")
  lev<-ord<-L<-U<-lev<-NULL
  
  sanity.init(method="agraph")
  create.tunes(method="agraph")
  obj<-run.model(method="agraph")

  obj@.fitinfo$fit[ord]=obj@.fitinfo$fit
  slot(obj,"fit")<-obj@.fitinfo$fit
  if(est.only){
    return(obj@fit)
  }
  obj@.fitinfo$weights[ord]<-obj@.fitinfo$weights
  muy=obj@.fitinfo$resp.info[1]
  ssy=obj@.fitinfo$resp.info[2]
  
  slot(obj,".call")<-cl
  slot(obj,".respinfo")<-list(L=as.vector(L),U=as.vector(U),lev=lev,muy=muy,ssy=ssy)
  obj@.fitinfo$resp.info<-NULL

  slot(obj,"lparm")<-obj@.fitinfo$lam
  slot(obj,"gparm")<-obj@.fitinfo$gam
  xmatrix(obj)=NULL
  obj
}

## Generic functions for class agraph are defined next.

setMethod("show", signature(object = "agraph"),function(object){
  if(class(object)!="agraph"){
    stop("Error: object is not of type agraph")
  }
  
  cl <- object@.call
  ind=grep("agraph",cl)
  if(length(ind)==0){
    cat("Empty agraph object\n")
    return()
  }
  
  
  type=object@type
  adims<-dim(object)
  
  n=adims[1]
  m=adims[2]
  cv.err=object@.cv_str$opt.row[4]
  edf=as.numeric(measures(object)[3])
  df=m-edf
  if(df<0)df=0.0
  
  cat("Anchor Graph Laplacian (agraph) with (n,|L|)=(",n,",",m,") or ",round(m/n*100,0),"% labeled")
  cat("\n\nPerformance Estimates:\nk-CV: ",round(cv.err,3)," GCV: ",round(as.numeric(measures(object))[2],3)," DF: ",round(df,3))
  cat("\n\nFit Estimates:\n")
  cat("Lagrangians: ",round(lparm(object),3)[1]," ",round(gparm(object),3)," ")
  if(object@.fitinfo$xdat){
    cat("Safe-Lagrangian",round(lparm(object),3)[2]," ")
  }
  cat("\n\n")
  invisible(object)
})

setMethod("predict", signature(object = "agraph"),function(object,xnew,gnew,type=c("vector","response","prob"),...){
  if(class(object)!="agraph"){
    stop("Error: object must be of type agraph")
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
    r1<-which(as.numeric(object@.fitinfo$x.vars)>0)
    if(length(r1)>0){
      xnew=x.scaleL(xnew,sanity=TRUE,sanity.only=TRUE)
      if(!is.null(SemiSupervised::xmatrix(object))){
        xnew=scale(xnew,center=attr(slot(object,"xmatrix"),"scaled:center"),scale=attr(slot(object,"xmatrix"),"scaled:scale"))
      }
      p=dim(xnew)[2]
      if(p!=length(object@.fitinfo$x.vars)){stop("Error: supplied xnew does not have correct number of variables")}
    }else{
      xdat=FALSE
    }
  }else{
    if(xdat)stop("Error: Need xnew since model was fit in safe mode")
  }
  if(!missing(gnew)){
    if(class(gnew)!="anchor")stop("Error: gnew supplied but not of class 'anchor'")
  }else{
    
    if(!default.mode){
      if(missing(xnew))stop("Error: Must input xnew or gnew")
      gnew=AnchorGraph(xnew,fit.g=gmatrix(object),metric=object@.fitinfo$metric,control=object@.control)
    }else{
      stop("Error: missing gnew")
    }
  }
  
  if(xdat){
    fit<-as.vector(cbind(gnew$Z,xnew[,r1,drop=FALSE])%*%object@.fitinfo$coef)
  }else{
    fit<-as.vector(gnew$Z%*%object@.fitinfo$coef)
  }
  
  fit<-fit+object@.respinfo$muy
  ctype=object@type
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
