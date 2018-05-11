################################################################################
##
##  Programmer: Mark Vere Culp
##
##  Date:        August 26, 2016
##
##  Description: These are functions defined for virtual class SemiSupervised. These
##               functions are exported graph building and that can be used
##               in common by a subset of approaches in the SemiSupervised library
##               and other libraries either under development or in use, e.g.,
##               the spa library.
##
################################################################################

## control parameters for the SemiSupervised and later package fits.


SemiSupervised.control<-function(normalize=TRUE,stability=NULL,k=NULL,nok=FALSE,
dissimilar=TRUE,l.eps=1e-5,l.thresh=25L,h.thresh=1e-5,U.as.anchor.thresh=600L,
U.as.anchor=TRUE,sig.est=TRUE,sig.frac=0.5,iter.max=1000L,
anchor.seed=100,sfac=5L,cn=4L,LAE.thresh=100L,LAE.eps=1e-4,
cv.fold=3L,cv.seed=100L,cv.cl=TRUE,cv.type="scv",cv.adjust=0.001){
  if(!is.null(stability)){
    if(!is.numeric(stability)||stability<0){
      stop("value of 'stability' must be >= 0")
    }
  }
  if(!is.numeric(h.thresh)||h.thresh<0){
    stop("value of 'h.thresh' must be >= 0")
  }
  if(!is.logical(U.as.anchor)){
    stop("value of 'U.as.anchor' must be either TRUE or FALSE")
  }
  if(U.as.anchor.thresh<0.0){
    U.as.anchor.thresh=floor(abs(U.as.anchor.thresh)+1)
  }
  if(!is.numeric(l.thresh)||l.thresh<0){
    stop("value of 'l.thresh' must be >= 0")
  }
  if(!is.numeric(l.eps)||l.eps<0){
    stop("value of 'l.eps' must be >= 0")
  }
  
  if (!is.numeric(cv.adjust) || cv.adjust< 0){
    stop("value of 'cv.adjust' must be >= 0")
  }
  if(cv.fold<0){
    stop("value of 'cv.fold' must be an integer>= 0")
  }
  cv.fold=floor(cv.fold)
  if(!is.numeric(cv.seed)){
    stop("value of 'cv.seed' must be numeric")
  }
  if(cv.type!="scv"){
    cv.type=1L
  }else{
    cv.type=0L
  }
  iter.max=ceiling(iter.max)
  if(iter.max<1){
    stop("value of 'iter.max' must be >1")
  }

  if(!is.logical(nok)){
    stop("value of 'nok' must be logical")
  }
  if(!is.null(k)){
    if(k<0){
      stop("value of 'k' must be greater than zero")
    }else{
      k=as.integer(floor(k))
    }
  }
  
  if(!is.logical(sig.est)){
    stop("value of 'sig.est' must be TRUE or FALSE")
  }
  if(sig.frac>1 || sig.frac<0){
    stop("value of 'sig.frac' must be in [0,1]")
  }
  if(!is.numeric(sfac)||sfac<0){
    stop("value of 'sfac' must be >= 0")
  }
  if(!is.numeric(LAE.thresh)||LAE.thresh<0){
    stop("value of 'LAE.thresh' must be an integer>= 0")
  }
  LAE.thresh=as.integer(ceiling(LAE.thresh))
  if(!is.numeric(LAE.eps)||LAE.eps<0){
    stop("value of 'LAE.eps' must be >= 0")
  }
  if (!is.numeric(cn) || cn< 0){
    stop("value of 'cn' must be >= 0")
  }
  return(list(normalize=normalize,stability=stability,k=k,nok=nok,
  dissimilar=dissimilar,l.eps=l.eps,l.thresh=l.thresh,h.thresh=h.thresh,
  U.as.anchor.thresh=U.as.anchor.thresh,U.as.anchor=U.as.anchor,
  sig.est=sig.est,sig.frac=sig.frac,iter.max=iter.max,
  anchor.seed=anchor.seed,sfac=sfac,cn=cn,LAE.thresh=LAE.thresh,
  LAE.eps=LAE.eps,cv.fold=cv.fold,cv.seed=cv.seed,cv.cl=cv.cl,
  cv.type=cv.type,cv.adjust=cv.adjust))
}


## This function is used to build the anchor graph.

AnchorGraph <- function(x,metric="cosine",anchor=NULL,fit.g=NULL,control=SemiSupervised.control()){
  if(missing(x)){
    stop("Error: an n X p matrix x of predictors is required but not given")
  }
  x=x.scaleL(x,sanity=TRUE,sanity.only=TRUE)
  n <- nrow(x)
  
  if(!is.null(anchor)){
    p=ncol(x)
    anchor=x.scaleL(anchor,sanity.only=TRUE)
    if(dim(anchor)[2]!=p){
      stop("Provided Anchor points are not the correct dimension")
    }
  }else{
    if(is.null(fit.g)){
      if(is.null(control$k)){
        if(n<2000){
          k=ceiling(n*0.3)
        }else{
          k=600
        }
      }else{
        k=control$k
      }
      control$k=k
      anchor=getAnchor(x,control)
    }else{
      anchor=fit.g$anchor
    }
  }
  Z=fitLAE(x,anchor,metric,control)
  if(is.null(fit.g)){
    T1 = crossprod(Z)
    cZ=colSums(Z)
    ind<-cZ>0
    if(sum(as.numeric(ind))<length(cZ)){
      cZ[!ind]=Inf
    }
    
    rL = T1-T1%*%diag(cZ^{-1})%*%T1
    ##if(is.nan(rL)>0 || is.na(rL)){
    ##  ind<-apply(Z,2,sum)>0
    ##  if(length(ind)>0){
    ##    Z=Z[,ind,drop=FALSE]
    ##    T1 = crossprod(Z)
    ##    rL = T1-T1%*%diag(colSums(Z)^-1)%*%T1
    ##  }
    ##  if(is.nan(rL)>0 || is.na(rL)){
    ##    stop("NA's or NaN's in reduced Laplcacian, try a smaller number of anchor points")
    ##  }
    ##}
    fit.g<-structure(list(Z=Z,rL=rL,anchor=anchor,metric=metric),class="anchor")
  }else{
    fit.g$Z=Z
  }
  return(fit.g)
}

## This function is used to get anchor points

"getAnchor"<-function(x,control){
  set.seed(control$anchor.seed)
  centers<-try(kmeans(x,centers=control$k,iter.max=control$iter.max),silent=TRUE)
  
  if(class(centers)=="try-error"){
    centers=x[sample(1:dim(x)[1],control$k),,drop=FALSE]
  }else{
    centers=centers$centers
  }
  return(centers)
}


## These functions build the k-NN/ epsilon graph. The symm functions are for
## training while the predict functions are for prediction.
## Note: If one were to invoke the predict function for K-NN in training they
##       will get a different answer. This is intentional. When training a
##       symmetric adjustment is used but such an adjustment can not be used
##       for prediction of a new case.


"knnGraph"<-function(x,y,k=6L,nok=FALSE,metric="cosine"){
  n=nrow(x)
  if(metric=="cosine"){
    if(missing(y)){
      x=cosineDist(x)
    }else{
      x=cosineDist(x,y)
    }
  }else{
    if(missing(y)){
      x=euclidDist(x)
    }else{
      x=euclidDist(x,y)
    }
  }
  
  if(!missing(y)){
    x=kgraph.predict(x,k=k,nok=nok)
  }else{
    x=kgraph.symm(x,k=k,nok=nok)
  }
  attr(x,"metric")=metric
  return(x)
}

"kgraph.predict"<-function(x,k=6L,nok=FALSE){
  if(!nok){
    p1=dim(x)[2]
    n1=dim(x)[1]
    for(i in 1:n1){
      x[i,order(x[i,])[(k+1):p1]]=Inf
    }
  }
  attr(x,"distance.graph")=2L
  if(is.null(attr(x,"metric")))attr(x,"metric")=NA
  return(x)
}



"kgraph.symm" <-function(x,k=6L,nok=FALSE){
  if(!nok){
    n<-dim(x)[1]
    m <- matrix(0, n, n)
    for (i in 1:n) {
      comp <- order(x[i, ])[k + 1]
      indx <- x[, i] <= x[i, comp]
      m[i, indx] <- 1
      m[indx, i] <- 1
    }
    m <- as.matrix(x * m)
    ind<-is.nan(m)
    if(sum(ind)>0){
      m[ind]=Inf
    }else{
      m[m == 0] = Inf
    }
    m[x == 0] = 0
  }else{
    m=x
  }
  attr(m,"distance.graph")=2L
  if(is.null(attr(x,"metric")))attr(m,"metric")=NA
  return(m)
}


"epsGraph"<-function(x,y,eps=0.2,noeps=FALSE,metric="cosine"){
  n=nrow(x)
  if(metric=="cosine"){
    if(missing(y)){
      x=cosineDist(x)
    }else{
      x=cosineDist(x,y)
    }
  }else{
    if(missing(y)){
      x=euclidDist(x)
    }else{
      x=euclidDist(x,y)
    }
  }
  if(!missing(y)){
    x=eps.predict(x,eps=eps,noeps=noeps)
  }else{
    x=eps.symm(x,eps=eps,noeps=noeps)
  }
  attr(x,"metric")=metric
  return(x)
}

"eps.predict"<-function(x,eps=0.2,noeps=FALSE){
  if(!noeps){
    p1=dim(x)[2]
    n1=dim(x)[1]
    
    m=x*apply(x<=eps,1,as.numeric)
    n=dim(m)[1]
    for(i in 1:n){
      m[m[,i]==0,i]=Inf
    }
    m[x==0]=0
  }else{
    m=x
  }
  attr(m,"distance.graph")=3L
  if(is.null(attr(x,"metric")))attr(m,"metric")=NA
  return(m)
}
"eps.symm" <-function(x,eps=0.2,noeps=FALSE){
  return(eps.predict(x,eps=eps,noeps=noeps))
}



## This is a cosine/eucldiean distance function.

"cosineDist" <- function(x,y=NULL){
  if(is.null(y)){
    mat=as.matrix(1-tcrossprod(x)/(sqrt(tcrossprod(rowSums(x^2)))))
  }else{
    mat=as.matrix(1-tcrossprod(x,y)/(sqrt(tcrossprod(rowSums(x^2),rowSums(y^2)))))
  }
  attr(mat,"metric")="cosine"
  attr(mat,"distance.graph")=1L
  mat
}

"euclidDist" <- function(x,y=NULL){
  if(is.null(y)){
    mat=as.matrix(dist(x))
  }else{
    x=t(x)
    y=t(y)
    aa <- matrix(colSums (x*x),nrow=1,byrow=T)
    bb <- matrix(colSums (y*y),nrow=1,byrow=T)
    ab <- crossprod(x,y)
    mat=sqrt(abs(repmat(t(aa),1,ncol(bb))+repmat(bb,ncol(aa),1)-2*ab))
  }
  attr(mat,"metric")="eulidean"
  attr(mat,"distance.graph")=1L
  mat
}

## A simple median imputation function.

"impute.median"<-function(object,...){
  handle.na<-function(object){ #plug the mean into missing values
    object[is.na(object)]<-median(object,na.rm=TRUE)
    return(object)
  }
  p<-dim(object)
  if(is.null(p)){ #special case: `object' is a vector
    return(handle.na(object))
  }
  object<-sapply(1:p[2],function(i)as.numeric(object[,i]))
  
  for(i in 1:p[2]){
    if(!is.factor(object[,i]))# if `object' is a factor ignore it
    object[,i]<-handle.na(object[,i])
  }
  return(object)
}

## These next routines are exported only so other packages can use them
## without using the ::: scope operator. This avoids a warning for
## other packages, e.g., R CMD CHECK --as-cran spa issues a warning
## if SemiSupervised:::x.scaleL is used, thus these function were exported.


## Formula internal commands.

"aG"<-function(x){
  return(x)
}


"dG"<-function(x,k=6L,nok=FALSE,metric=NULL){
  if(is.null(metric)){
    metric=attr(x,"metric")
    if(is.null(metric))stop("The original graph call must either use cosineDist or euclidDist to assign the metric")
    
  }
  if(!is.null(attr(x,"distance.graph"))){
    if(!nok){
      return(kgraph.symm(x,k=k))
    }
    return(x)
  }
  return(knnGraph(x,k=k,nok=nok,metric=metric))
}

"sG"<-function(x){
  attr(x,"metric")=NA
  as.matrix(x)
}



## This low level routine scales the data so that it is centered to the labeled means
## with variance one.  Also, the function can perform sanity checks on inputs.

"x.scaleL"<-function(x,L,sanity=TRUE,sanity.only=FALSE){
  n=nrow(x)
  if(sanity){
    if(is.null(n)){
      n=1
    }
    if(n>1){
      x=impute.median(x)
    }else{
      x=matrix(x,1)
      ind<-which(is.na(x[1,]))
      if(length(ind)>0){
        x[,ind]=0.0
      }
    }
    if(is.data.frame(x)){
      x=model.matrix(~.-1,x)
    }else{
      if(is.vector(x))x=t(t(as.numeric(x)))
      x=as.matrix(x)
    }
  }
  if(sanity.only)return(x)
  if(length(attr(x,"scaled:center"))>0){
    return(x)
  }
  if(missing(L)){
    L=1:n
  }
  ms<-as.vector(apply(x[L,,drop=FALSE],2,mean))
  sos=sqrt(as.vector(apply(sweep(x,2,ms)[L,,drop=FALSE]^2,2,sum)))
  
  ind<-sos==0
  if(sum(ind)>0){
    sos[ind]=1.0;
  }
  return(scale(x,center=ms,scale=sos))
}


## This function is a cross validation function routine.
"cv.folds"<-function (n, folds = 3L){
  if(n<3L)stop("Cross Validation requires more than two observations.")
  if(n<4L){
    fls<-list(c(1L,2L),c(1L,3L),c(2L,3L))
  }else{
    while(!(n/folds>2)){
      folds=folds-1
    }
    if(folds<2)folds=2
    fls<-split(sample(1:n), rep(1:folds, length = n))
  }
  fls
}

## This function determines the h.est local kernel parameter.
## This was minorly modified from the kernlab package.

"h.est"<-function (Dist,sym=TRUE,thresh){
  n <- dim(Dist)[1]
  
  dv<-sort(Dist[is.finite(Dist)])
  nz<-length(dv[dv==0])
  dv<-dv[dv>0]
  ldv<-length(dv)
  
  if(ldv>0){
    if(sym){
      dv<-dv[(1:(ldv/2))*2-1]
    }
  }
  rg<-sort(c(0.12,as.vector(quantile(dv,c(0.1,0.5,0.9)))))
  rg[rg>thresh]
}

## These are formula extraction routines.

"extract.a"<-function(form,env){
  g.form=NULL
  x.form=NULL
  a.form=NULL
  dg=1L
  frm=paste(form)
  rhs=strsplit(frm[[3]],"\\+")[[1]]
  ind<-grep("aG\\(",rhs)
  if(length(ind)>0L){
    g.form=as.formula(paste0(frm[2],frm[1],rhs[ind],collapse="+"),env)
    a.form=g.form[[3]]
    g.form=as.formula(paste(frm[2],frm[1],"1",sep=""),env)
  }else{
    dg=2L
  }
  if(length(ind)<1L){
    x.form=as.formula(paste0(frm[2],frm[1],paste0(rhs,collapse="+")),env)
  }
  rhs=rhs[-ind]
  if(length(rhs)>0L){
    x.form=as.formula(paste0(frm[2],frm[1],paste0(rhs,collapse="+")),env)
  }else{
    x.form=NULL
  }
  return(list(dg,x.form,g.form,TRUE,a.form))
}

"extract.xg"<-function(form,env){
  ctrl=TRUE
  g.form=NULL
  x.form=NULL
  dg=1L
  frm=paste(form)
  rhs=strsplit(frm[[3]],"\\+")[[1]]
  ind<-grep("dG\\(",rhs)
  if(length(ind)<1L){
    ind<-grep("sG\\(",rhs)
    dg=2
    if(length(ind)>0L)ctrl=FALSE
  }
  if(length(ind)>0L){
    g.form=as.formula(paste0(frm[2],frm[1],rhs[ind],collapse="+"),env)
  }else{
    dg=3
  }
  if(length(ind)<1L){
    x.form=as.formula(paste0(frm[2],frm[1],paste0(rhs,collapse="+")),env)
  }
  rhs=rhs[-ind]
  if(length(rhs)>0L){
    x.form=as.formula(paste0(frm[2],frm[1],paste0(rhs,collapse="+")),env)
  }else{
    x.form=NULL
  }
  return(list(dg,x.form,g.form,ctrl))
}

"extract.g"<-function(form,env){
  ctrl=TRUE
  g.form=NULL
  dg=1L
  frm=paste(form)
  rhs=strsplit(frm[[3]],"\\+")[[1]]
  ind<-grep("dG\\(",rhs)
  if(length(ind)<1L){
    ind<-grep("sG\\(",rhs)
    dg=2
    ctrl=FALSE
  }
  if(length(ind)>0L){
    g.form=as.formula(paste0(frm[2],frm[1],rhs[ind],collapse="+"),env)
  }else{
    dg=3
  }
  return(list(dg,NULL,g.form,ctrl))
}

## A interpolation prediction routine.

"inter.predict" <-function(y,W,x=NULL,bet=NULL,pow,fhat=FALSE){
  np=dim(W)[1]
  if(!is.null(x)){
    n=dim(x)[1]
    if(n==1){
      f1=sum(x*bet)
    }else{
      f1=as.vector(x%*%bet)
    }
  }else{
    f1=rep(0.0,np)
  }
  f<-sapply(1:np,function(i){a=abs(W[i,]/(1-W[i,]));a=a^(pow);sum(a*y)/sum(a)})
  
  ind<-which(is.nan(f))
  if(length(ind)>0){
    for(i in ind){
      val<-as.numeric(W[i,]==1)
      f[i]=sum(val*y)/sum(val)
    }
  }
  ind<-which(is.nan(f)|is.na(f))
  if(length(ind)>0)f[ind]=0.0;
  fit=as.numeric(f1+f)
  if(fhat){
    return(data.frame(fit=fit,f1=f,f2=f1))
  }
  return(fit)
}
