fvec<-c(NA,n,m,p,NA,NA,NA,NA,NA)
ctrl.s4pm1<-s4pm.control(normalize=TRUE)
ctrl.agraph<-agraph.control()
if(n<1001){
  ctrl.s4pm1$U.as.anchor.thresh=n+1
  ctrl.agraph$U.as.anchor.thresh=n+1
}
for(j1 in 1:nruns){
  set.seed(j1)
  L<-sort(sample(1:n,m))
  U<-setdiff(1:n,L)
  
  yval<-y
  yval[U]<-NA
  mu<-mean(y[L])
  msy<-sqrt(mean( (y[L]-mu)^2))
  ind<-which(apply(x[L,,drop=FALSE],2,var)>0)
  dat<-data.frame(y=yval,x)[,ind]
  
  fvec[1]<-j1
  fvec[8]<-system.time(g1<-s4pm(y~.,data=dat,type=type,control=ctrl.s4pm1))[3]
  if(type=="r"){
    fvec[5]<-sqrt(mean((y[U]-g1@fit[U])^2))/msy
  }else{
    f<-factor(sign(g1@fit),levels=c(-1,1))
    tab<-table(f[U],y[U])
    fvec[5]<-1-sum(diag(tab))/sum(tab)
  }
  
  fvec[9]<-system.time(g1<-agraph(y~.,data=dat,type=type,control=ctrl.agraph))[3]
  if(type=="r"){
    fvec[6]<-sqrt(mean((y[U]-g1@fit[U])^2))/msy
  }else{
    f<-factor(sign(g1@fit),levels=c(-1,1))
    tab<-table(f[U],y[U])
    fvec[6]<-1-sum(diag(tab))/sum(tab)
  }
  
  cv.errs<-rep(0.0,lalp)
  t.errs<-rep(0.0,lalp)
  if(type=="r"){
    for(i in 1:lalp){
      g1<-cv.glmnet(as.matrix(x[L,ind,drop=FALSE]), yval[L],alpha=alpha[i])
      f<-predict(g1,newx=x[U,ind,drop=FALSE])
      t.errs[i]<-sqrt(mean((y[U]-f)^2))/msy
      cv.errs[i]<-min(g1$cvm)
    }
  }
  if(type=="c"){
    for(i in 1:lalp){
      g1<-cv.glmnet(as.matrix(x[L,ind,drop=FALSE]), (yval[L]+1)/2,alpha=alpha[i],family="binomial")
      f<-as.vector(predict(g1,newx=x[U,ind,drop=FALSE],type="class"))
      tab<-table(f,(y[U]))
      t.errs[i]<-1-sum(diag(tab))/sum(tab)
      cv.errs[i]<-min(g1$cvm)
    }
  }
  fvec[7]<-t.errs[which.min(cv.errs)]
  
  
  cat("fn: ",nms," iteration=",j1," time=",(proc.time()-t1)[3]/60,"\n")
  final[iter,-1]<-fvec
  final[iter,1]<-nms
  iter<-iter+1
  write.csv(final,filenm,quote=FALSE,row.names=FALSE)
}
