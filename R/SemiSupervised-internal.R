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

## These are internal graph functions.

"lap"<-function(g,stab=1e-3,normalize=FALSE){
  if (!is.numeric(stab) || stab< 0){
    stop("value of 'stab' must be >= 0")
  }
  x=list()
  g=as.matrix(g)
  m=dim(g)
  if(m[1]!=m[2]){
    stop("Graph must have same number of rows and columns (i.e. g must be an adjaceny matrix)\n")
  }
  obj<-new("lapGraph")
  
  x[[1]]=g
  x[[2]]=as.double(stab)
  x[[3]]=as.numeric(normalize)
  x[[4]]=m[1]
  slot(obj,"graph")=x
  obj
}

"g.scaleL"<-function(g){
  g@graph[[1]]=as.numeric(g@graph[[1]])
  g
}

## A repmat function that is useful.

"repmat" <- function(X,m,n){
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

## The low level LAE fitting routine.

"fitLAE"<-function(x,anchor,metric,control=SemiSupervised.control()){
  n=nrow(x)
  if(is.null(n)){
    x=x.scaleL(x,sanity=TRUE,sanity.only=TRUE)
    n=1
  }
  d <- dim(anchor)[2]
  m <- dim(anchor)[1]
  
  sfac=control$sfac
  if(sfac>m){
    sfac=m
  }
  cn=control$cn
  thresh=control$l.thresh
  eps=control$l.eps
  
  if(metric=="cosine"){
    Dis<-cosineDist(x,anchor)
  }else{
    Dis <- euclidDist(x,anchor)
  }
  pos = t(sapply(1:n,function(i)order(Dis[i,])[1:sfac]))
  if(sfac==1)pos=t(pos)
    
  storage.mode(pos)<-"integer"

  ind <- (pos-1)*n+repmat(matrix(c(1:n),ncol=1),1,sfac) ##CHECK
  storage.mode(pos)<-"integer"
   val<-.Call("LAE",t(x),t(anchor),pos,as.integer(cn),as.integer(thresh),as.numeric(eps),PACKAGE="SemiSupervised")
  Z1 <- rep(0.0,m*n)
  Z1 [as.numeric((pos-1)*n)+rep(1:n,sfac)] <- as.numeric(t(matrix(val,sfac,n)))
  return(matrix(Z1,n,m))
}

## Check sanity on initialization
sanity.init<-function(envir=sys.frame(-1),method="s4pm"){
  if(any(envir$graph=="",na.rm=TRUE)){
    stop(paste("Error in graph:  must be supplied"))
  }
  if(!is.na(match(method,c("s4pm","jtharm","mex","mreg")))){
    if(class(envir$graph)=="lapGraph"){
      n=envir$graph[[4]]
      graph=matrix(envir$graph[[1]],n,n)
    }else{
      graph=x.scaleL(envir$graph,sanity.only=TRUE)
      n=nrow(graph)
      n1=ncol(graph)
      if(n1!=n){
        stop("Error: Graph must be symmetric.")
      }
    }
  }
  if(!is.na(match(method,c("amex","agraph")))){
    if(class(envir$graph)!="anchor")stop("Error: The graph must be created as an Anchor Graph of class `anchor'\n")
    graph=envir$graph
    n=dim(graph$Z)[1]
  }
  assign("n",n,envir=envir)
  if(any(envir$y=="",na.rm=TRUE)){
    stop("Error in y: a response (either continuous or binary) must be provided")
  }
  y<-envir$y
  m=length(y)
  if(m==n){
    L=as.numeric(which(!is.na(y)))
    U=as.numeric(which(is.na(y)))
    m=length(L)
    y=y[L]
  }else{
    if(m>n){
      stop("Error: |y|> dim(graph)[1]")
    }
    L=1:m
    U=(m+1):n
  }
  ord=c(L,U)
  assign("m",m,envir=envir)
  assign("ord",ord,envir=envir)
  assign("L",L,envir=envir)
  assign("U",U,envir=envir)
  
  type=envir$type
  if(length(type)>1){
    type="r"
    if(is.factor(y))type="c"
  }
  lev=NA
  if(type=="r"){
    y=as.numeric(y)
    if(is.null(envir$control$stability))envir$control$stability=1e-3
    type=0L
  }else{
    lev<-paste(c(-1,1))
    if(is.factor(y))lev=levels(y)
    if(nlevels(y)>2){
      stop(paste("Error in y: Currently the response must be binary"))
    }
    y=as.numeric(as.factor(y))
    y=c(-1.0,1.0)[y]
    if(is.null(envir$control$stability))envir$control$stability=1e-5
    type=1L
  }
  assign("type",type,envir=envir)
  assign("lev",lev,envir=envir)
  assign("y",y,envir=envir)
  
  if(!any(envir$weights=="",na.rm=TRUE)){
    weights=envir$weights
    if(is.vector(weights) & length(weights)==n){
      weights=envir$weights[ord]
    }else{
      stop("Error: weights must be a vector of length ",n);
    }
  }else{
    weights=rep(1.0,n)
  }
  if(!is.na(match(method,c("s4pm","jtharm","mex","mreg")))){
    graph=graph[ord,ord,drop=FALSE]
  }
  if(!is.na(match(method,c("amex","agraph")))){
    graph$Z=graph$Z[ord,,drop=FALSE]
  }
  assign("graph",graph,envir=envir)
  assign("weights",weights,envir=envir)
  xdat=FALSE
  if(!is.na(match(method,c("s4pm","agraph")))){
    check<-any(envir$x=="",na.rm=TRUE)
    xdat=FALSE
    if(!check){
      if(!is.null(envir$x)){
        x=x.scaleL(envir$x,sanity.only=TRUE)
        x=x[ord,,drop=FALSE]
        assign("x",x,envir=envir)
        xdat=TRUE
      }
    }
  }
  assign("xdat",xdat,envir=envir)
  is.kern=FALSE
  if(!is.na(match(method,c("mreg","mex","amex")))){
    if(any(envir$K=="",na.rm=TRUE)){
      if(method=="amex")stop("Error in K: Kernel gramm matrix must be provided. Use `agraph' if you want to fit the anchor graph only model.")
      stop("Error in K: Kernel gramm matrix must be provided. Use `s4pm' if you want to fit the graph only model.")
     }
    K=x.scaleL(envir$K,sanity.only=TRUE)
    if(nrow(K)!=n){
      stop("Error in K: K must be a kernel marix with ",n," rows")
    }
    if(ncol(K)!=n){
      stop("Error in K: K must be a kernel marix with ",n," columns")
    }
    K=K[ord,ord,drop=FALSE]
    assign("K",K,envir=envir)
    is.kern=TRUE
  }
  assign("is.kern",is.kern,envir=envir)
  is.anchor=FALSE
  if(!is.na(match(method,c("agraph","amex")))){
    is.anchor=TRUE
  }
  assign("is.anchor",is.anchor,envir=envir)
  
}
## Create tuning parameters

create.tunes<-function(envir=sys.frame(-1),method="s4pm"){
  if(!is.na(match(method,c("s4pm","jtharm","mex","mreg")))){
    if(envir$control$dissimilar){
      if(any(envir$hs=="",na.rm=TRUE)){
        hs=as.vector(h.est(envir$graph,thresh=envir$control$h.thresh))
      }else{
        hs=as.numeric(na.omit(envir$hs))
        hs=abs(hs)
      }
    }else{
      hs=NA
    }
    lhs=length(hs)
    assign("hs",hs,envir=envir)
  }else{
    lhs=1
    hs=NA
  }
 
  if(any(envir$gams=="",na.rm=TRUE)){
    gams=c(0.0,0.001,0.01)
    if(method=="jtharm"){
         gams=c(0.001,0.01,0.1)
    }
  }else{
    gams=as.numeric(na.omit(envir$gams))
    gams=abs(gams)
  }
  
  if(any(envir$lams=="",na.rm=TRUE)){
    lam1=c(0.01,0.1,1.0,2.0,10.0)
    xdat=TRUE
    if(!is.na(match(method,c("s4pm","agraph")))){
      xdat=envir$xdat
    }
    lam2=c(0.0,0.1,1.0,2.0,10.0)
    if(!xdat){
      lam2=NA
    }
    if(envir$is.kern){
      lam2=c(0.1,1.0,2.0,10.0)  ##require an ambient penatly for a kernel
    }
    
    l1=length(lam1)
    l2=length(lam2)
    lams=matrix(0,l1*l2,2)
    k3=1
    for(i1 in 1:l1){
      for(j1 in 1:l2){
        lams[k3,]=c(lam1[i1],lam2[j1])
        k3=k3+1
      }
    }
  }else{
    lams=na.omit(envir$lams)
    if(is.vector(lams)){
      lams=t(as.matrix(lams))
    }
    if(dim(lams)[2]<2){
      val<-2-dim(lams)[2]
      for(i in 1:val){
        lams=cbind(lams,0.1)
      }
    }
  }
  llam=dim(lams)[1]
  lgam=length(gams)
  lhs=length(hs)
  assign("lgam",lgam,envir=envir)
  assign("lhs",lhs,envir=envir)
  assign("gams",gams,envir=envir)

  tunes<-lapply(hs,function(i){
    a=1
    if(lhs==1){
      a=0
    }
    tunes<-list()
    k3=1
    for(j1 in 1:llam){
      for(k1 in 1:lgam){
        tvec=NULL
        if(!envir$is.anchor){
          tvec<-c(tvec,i)
        }
        if(method=="jtharm"){
          tvec<-c(tvec,lams[j1,1])
        }else{
          tvec<-c(tvec,lams[j1,])
        }
        tvec<-c(tvec,gams[k1])
        tunes[[k3]]<-tvec
        k3<-k3+1
      }
    }
    do.call("rbind",tunes)
  })
  assign("tunes",tunes,envir=envir)
}

## Generic run the model function
run.model<-function(envir=sys.frame(-1),method="s4pm",call.method=NULL){
  control=envir$control
  m=envir$m
  n=envir$n
  
  if(is.null(call.method)){
    call.method=method
  }
  fnc<-paste(call.method,".cv",sep="")

  inpts<-list()
  inpts$type=envir$type
  inpts$control=envir$control
  inpts$y=envir$y
  
  if(envir$is.anchor){
    inpts$W=envir$graph$rL
    inpts$tune=envir$tunes[[1]]
    x.ind=1:n
    if((control$cv.type==0L) &&(m<n)){
      x.ind=1:m
    }
   inpts$z=envir$graph$Z[x.ind,,drop=FALSE]
   inpts$weights<-envir$weights[x.ind]
    if(!envir$is.kern && method!="jtharm"){
      if(envir$xdat){
        inpts$x<-envir$x[x.ind,,drop=FALSE]
      }
    }
    if(envir$is.kern){
      inpts$K<-envir$K[x.ind,x.ind,drop=FALSE]
    }
    g1<-do.call(fnc,inpts)
    fin.tunes<-g1$tunes
    fin.folds<-g1$folds
  }else{
    fin.g1<-lapply(envir$tunes,function(tunes.use){
      hs=tunes.use[1,1]
      if(!control$dissimilar){
        lw=g.scaleL(lap(envir$graph,control$cv.adjust,control$normalize))@graph
      }else{
        if(is.na(hs)){
          stop("Error: hs variable is NA.")
        }
        lw=g.scaleL(lap(exp(-envir$graph/hs),control$stability,control$normalize))@graph
      }
      x.ind=1:n
      if(control$cv.type==0L && (m<n)){
        gs1=.Call("lgraph",list(lw),as.integer(m),as.integer(n),1L,PACKAGE="SemiSupervised")
        if(gs1[[2]][1]==0L){
          lw[[1]]=gs1[[1]]
          lw[[4]]=as.integer(m)
          x.ind=1:m
        }
      }
      inpts$weights<-envir$weights[x.ind]
      if(!envir$is.kern && method!="jtharm"){
        if(envir$xdat){
          inpts$x<-envir$x[x.ind,,drop=FALSE]
        }
      }
      if(envir$is.kern){
        inpts$K<-envir$K[x.ind,x.ind,drop=FALSE]
      }
      inpts$graph=lw
      inpts$tune=tunes.use
      do.call(fnc,inpts)
    })
    fin.tunes<-do.call("rbind",lapply(fin.g1,function(i)i$tunes))
    fin.folds<-fin.g1[[1]]$folds
  }


  nms<-dimnames(fin.tunes)[[2]]
  ind<-match("cv.errs",nms)
  flag=FALSE
  mn=NULL
  if(is.na(ind)){
    flag=TRUE
  }else{
    mn=which.min(fin.tunes[,ind])
  }
  if(is.null(mn)){
    flag=TRUE
  }else{
    tvec=as.vector(fin.tunes[mn,])
    tind<-list()
    tind[[1]]<-match("hs",nms)
    tind[[2]]<-match("lam",nms)
    tind[[3]]<-match("lam1",nms)
    tind[[4]]<-match("lam2",nms)
    tind<-do.call("c",tind)
    tind<-na.omit(tind)
    tvec<-tvec[tind]
  }
  if(!flag){
     inpts<-list()
     inpts$weights<-envir$weights
     inpts$type=envir$type
     inpts$control=envir$control
     inpts$y=envir$y
     
     if(!envir$is.kern && method!="jtharm"){
       if(envir$xdat){
         inpts$x<-envir$x
       }
     }
     if(envir$is.kern){
       inpts$K<-envir$K
     }
     if(envir$is.anchor){
       inpts$z=envir$graph$Z
       inpts$W=envir$graph$rL
     }else{
       if(!control$dissimilar){
         lw=g.scaleL(lap(envir$graph,control$stability,control$normalize))@graph
       }else{
         if(is.na(tvec[1])){
           stop("Error: hs variable is NA.")
         }
         lw=g.scaleL(lap(exp(-envir$graph/tvec[1]),control$stability,control$normalize))@graph
       }
       inpts$graph=lw
     }
     if(control$cv.type==0L && (m<n)){
       mat<-list()
       for(i in 1:envir$lgam){
         mat[[i]]<-c(tvec,envir$gams[i])
       }
       inpts$tunes=do.call("rbind",mat)
       g1<-do.call(fnc,inpts)
       mn=which.min(g1$tunes[,ind])
       
       if(length(mn)==0){
         flag=TRUE
       }else{
         tvec=g1$tunes[mn,,drop=FALSE]
         cv.str=list(stagewise=fin.tunes,classic=g1$tunes,folds=g1$folds,opt.row=NA)
       }
     }else{
       tvec=fin.tunes[mn,,drop=FALSE]
       cv.str=list(stagewise=NULL,classic=fin.tunes,folds=fin.folds,opt.row=NA)
     }
  }
  fnc<-paste(call.method,".fit",sep="")
  
  if(!flag){
    tvec=as.vector(tvec)
    cv.str$opt.row=tvec
    
    inpts$tunes=tvec
    model<-do.call(fnc,inpts)
    if(method=="jtharm"){
      model$lam=tvec[2]
      model$gam=tvec[3]
    }
    if(envir$is.anchor){
      model$lam=tvec[1:2]
      model$gam=tvec[3]
    }
    if((!envir$is.anchor)&&(method!="jtharm")){
      model$lam=tvec[2:3]
      model$gam=tvec[4]
    }
    model$h=tvec[1]
  }else{
    model=list(fit=rep(0.0,n))
  }
  if(envir$est.only){
    return(model$fit)
  }
  obj<-new(method)
  if(flag)return(obj)
  slot(obj,".cv_str")<-cv.str
  slot(obj,".control")<-control
  model$metric=NA
  model$colnames=NULL
  model$as.default=TRUE
  slot(obj,".fitinfo")<-model
  slot(obj,"type")<-c("r","c")[envir$type+1]
  ymatrix(obj)=envir$y[1:m]
  gmatrix(obj)=NULL
  ssy=model$resp.info[2]
  
  if(envir$type==0L){
    err=mean((envir$y[1:m]-model$fit[1:m])^2)/ssy^2
    gcv=err/(1-model$edf/m)^2
    conf=c(sqrt(err),gcv,model$edf)
  }else{
    err=mean(log(1+exp(-2*((envir$y[1:m])*(model$fit[1:m])))))
    gcv=err/(1-model$edf/m)^2
    conf=c(sqrt(err),gcv,model$edf)
  }
  slot(obj,"measures")<-conf
  obj
}


