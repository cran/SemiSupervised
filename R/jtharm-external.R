################################################################################
##
##  Programmer: Mark Vere Culp
##
##  Date:        August 26, 2016
##
##  Description: Joint Harmonic Graph functionality. These functions define
##               jtharm graph object. All functions in this file are exported.
##               These are all documented in the R help files.
##
################################################################################


## Definition of jtharm methods. These are all documented.

setGeneric("jtharm", function(x, ...) standardGeneric("jtharm"))
setMethod("jtharm",signature(x="formula"),function(x,data,metric=c("cosine","euclidean"),...,
est.only=FALSE,control=SemiSupervised.control()){
  cl <- match.call()
  cl[[1]]<-as.name("jtharm")
  
  env.parent=parent.frame()
  mf <- match.call(expand.dots = FALSE)
  
  m <- match(c("x", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  eg<-extract.xg(x,env.parent)
  control$dissimilar=eg[[4]]
  if(eg[[1]]>2L){
    mt<-terms(x,data=data)
    attr(mt,"intercept")=0
    mf$formula<-mt
    mf$na.action=na.pass
    eg[[3]]=list(terms=mt,frame=mf)
    mf.x<-eval(mf,env.parent)
    y=model.response(mf.x,"any")
    dat=model.matrix(mt,mf.x)
    if(missing(metric))metric="cosine"
    control$dissimilar=TRUE
  }else{
    dat=NULL
    mt<-terms(eg[[3]])
    attr(mt,"intercept")=0
    mf$formula<-mt
    mf$na.action=na.pass
    eg[[3]]=list(terms=mt,frame=mf)
    mf.g<-eval(mf,env.parent)
    metric=NA
    y=model.response(mf.g,"any")
    control$U.as.anchor=FALSE
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
 if(is.null(control$k)){
   control$k=6L
 }
  if(!control$U.as.anchor){
     if(eg[[1]]<3L){
      graph=as.matrix(model.matrix(mt,mf.g))
    }else{
      dat=x.scaleL(dat,L)
      graph=as.matrix(dG(dat,k=control$k,nok=control$nok,metric=metric))
    }
    obj=jtharm.default(graph,y,...,est.only=est.only,control=control)
    if(est.only){
      return(obj)
    }
    if(eg[[1]]>2L){
      slot(obj,".fitinfo")$as.default=FALSE
    }
  }else{
    dat=x.scaleL(dat,L)
    anchor=getAnchor(dat,list(anchor.seed=control$cv.seed,k=control$U.as.anchor.thresh,iter.max=1000L))
    ndat=rbind(dat[L,,drop=FALSE],anchor)
    tgraph=as.matrix(dG(ndat,k=control$k,nok=control$nok,metric=metric))
    obj=jtharm.default(tgraph,y[L],...,est.only=FALSE,control=control)
    gmatrix(obj)<-tgraph
  }
  slot(obj,".call")<-cl
  eg[[5]]=env.parent
  slot(obj,".terminfo")<-eg
  slot(obj,".fitinfo")$metric=metric
  if(control$U.as.anchor){
    graph=knnGraph(dat,ndat,k=control$k,nok=control$nok,metric=metric)
    fhat<-predict(obj,gnew=graph,type="vector")
    slot(obj,"fit")<-fhat
    slot(obj,".fitinfo")$dims[1]=n
    slot(obj,".respinfo")$L=L
    slot(obj,".respinfo")$U=U
    slot(obj,".fitinfo")$as.default=FALSE
  }
  gmatrix(obj)<-graph
  xmatrix(obj)<-dat
  return(obj)
})

setMethod("jtharm",signature(x="matrix"),function(x,...){
  obj=jtharm.default(graph=x,...)
  gmatrix(obj)=x
  return(obj)
})


## jtharm default function. This provides the main functionality for the jtharm.
## This function is exported but is not expected to be called directly.



"jtharm.default"<-function(graph,y,weights,hs,lams,gams,type=c("r","c"),est.only=FALSE,control=SemiSupervised.control()){
  cl <- match.call()
  cl[[1]]<-as.name("jtharm")
  
  lev<-ord<-muy<-ssy<-L<-U<-lev<-NULL
  sanity.init(method="jtharm")
  create.tunes(method="jtharm")
  control$normalize=FALSE
  obj<-run.model(method="jtharm")
  obj@.fitinfo$fit[ord]=obj@.fitinfo$fit
  slot(obj,"fit")<-obj@.fitinfo$fit
  if(est.only){
    return(obj@fit)
  }
  obj@.fitinfo$fitted_response[ord]=obj@.fitinfo$fitted_response
  obj@.fitinfo$weights[ord]<-obj@.fitinfo$weights
  muy=obj@.fitinfo$resp.info[1]
  ssy=obj@.fitinfo$resp.info[2]
  
  slot(obj,".call")<-cl
  slot(obj,".respinfo")<-list(L=as.vector(L),U=as.vector(U),lev=lev,muy=muy,ssy=ssy)
  obj@.fitinfo$resp.info<-NULL
  
  slot(obj,"hparm")<-obj@.fitinfo$h
  slot(obj,"lparm")<-obj@.fitinfo$lam
  slot(obj,"gparm")<-obj@.fitinfo$gam
  obj
}


setMethod("show", signature(object = "jtharm"),function(object){
  if(class(object)!="jtharm"){
    stop("Error: object is not of type jtharm")
  }
  cl <- object@.call
  ind=grep("jtharm",cl)
  if(length(ind)==0){
    cat("Empty jtharm object\n")
    return()
  }
  
  n=dim(object)[1]
  m=dim(object)[2]
  cv.err=object@.cv_str$opt.row[4]
  dis=object@.control$dissimilar
  edf=measures(object)[3]
  df=m-edf
  if(df<0)df=0.0
  
  cat("Joint Harmonic Fit with (n,|L|)=(",n,",",m,") or ",round(m/n*100,0),"% labeled\n")
  cat("\nPerformance Estimates:\nk-CV: ",round(cv.err,3)," GCV: ",round(measures(object)[2],3)," DF: ",round(df,3))
  cat("\n\nFit Estimates:\n")
  if(dis){
    cat("Graph Kernel h: ",round(hparm(object),3)," ")
  }
  cat("Lagrangian: ",round(lparm(object),3),"    Safe-Lagrangian: ",round(gparm(object),3)," ")
  cat("\n\n")
  invisible(object)
})

setMethod("predict", signature(object = "jtharm"), function(object,xnew,gnew,type=c("vector","response","prob"),pow=1,...){
  if(class(object)!="jtharm"){
    stop("Error: object must be of type jtharm.")
  }
  if(missing(type)){
    type="response"
  }

  default.mode=object@.fitinfo$as.default
  
  form.struct<-object@.terminfo
  if(!is.null(form.struct)){
    if(missing(gnew)){
      if(form.struct[[1]]>1){
        tt=form.struct[[3]]$terms
        Terms <- delete.response(tt)
        m <- model.frame(Terms, xnew, na.action = na.pass)
        xnew <- model.matrix(Terms, m)
      }
    }
  }
  
  if(!missing(xnew)){
    xnew=x.scaleL(xnew,sanity=TRUE,sanity.only=TRUE)
    if(!is.null(xmatrix(object))){
      xnew=scale(xnew,center=attr(xmatrix(object),"scaled:center"),scale=attr(xmatrix(object),"scaled:scale"))
    }
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
    stop("Error: The dims gnew do not match: Expect: (",dim(object)[1],") Inputted: (",np,").")
  }
  ctype=object@type
  yvec=object@fit
  if(object@.control$dissimilar) gnew=exp(-gnew/object@hparm)
  xnew=NULL
  bet=NULL
  fit=inter.predict(yvec,gnew,xnew,bet,pow)
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

