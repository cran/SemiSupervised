#include"Anchor_Reg.h"

SEXP AREG_CV(SEXP folds,SEXP omits,SEXP _Z,SEXP _X,SEXP W,SEXP y,SEXP case_wts,SEXP tunes,SEXP type,SEXP thresh,SEXP eps,SEXP _stab){
  SEXP adims=PROTECT(getAttrib(_Z, R_DimSymbol));
  int i=0,j=0,k=0,n=INTEGER(adims)[0],pz=INTEGER(adims)[1];
  int K=LENGTH(folds),i1=0,m1=0,px=0,qt=0,m=LENGTH(y);
  int ntunes,prot=0,lfold=0,v,v1;
  double tmp;
  
  if(!Rf_isNull(tunes)){
    ntunes=INTEGER(getAttrib(tunes, R_DimSymbol))[0];
  }else{
    ntunes=0;
  }
  if(!Rf_isNull(_X)){
    qt=1;
    px=INTEGER(getAttrib(_X, R_DimSymbol))[1];
  }
  SEXP X=PROTECT(duplicate(_X));
  SEXP Z=PROTECT(duplicate(_Z));
  SEXP yval=PROTECT(duplicate(y));
  SEXP wts=PROTECT(duplicate(case_wts));
  SEXP indx=PROTECT(allocVector(INTSXP,n));
  SEXP fit=PROTECT(allocVector(REALSXP,n));
  SEXP _m=PROTECT(allocVector(INTSXP,1));
  SEXP convs=PROTECT(allocVector(INTSXP,ntunes));
  SEXP errs=PROTECT(allocVector(INTSXP,ntunes));
  
  SEXP cv_errs=PROTECT(allocVector(REALSXP,ntunes));
  SEXP tvec=PROTECT(allocVector(REALSXP,3));
  SEXP cv_mode=PROTECT(allocVector(INTSXP,1));
  SEXP mod=PROTECT(allocVector(VECSXP,2));
  SEXP final=PROTECT(allocVector(VECSXP,3));
  SEXP resp_info=PROTECT(allocVector(REALSXP,2));
  
  prot=16;
  for(i=0;i<LENGTH(cv_errs);i++){
    REAL(cv_errs)[i]=0.0;
    INTEGER(convs)[i]=0;
    INTEGER(errs)[i]=0;
  }

  for(i=0;i<n;i++){
    REAL(fit)[i]=0.0;
  }
  INTEGER(cv_mode)[0]=1;
  for(k=0;k<K;k++){
    m1=0;
    lfold=LENGTH(VECTOR_ELT(folds,k));
    INTEGER(_m)[0]=lfold;

    for(i=0;i<n;i++){
      if(i<lfold){
        INTEGER(indx)[i]=INTEGER(VECTOR_ELT(folds,k))[i]-1;
      }else{
        INTEGER(indx)[i]=INTEGER(VECTOR_ELT(omits,k))[i-lfold]-1;
      }
    }
    for(i=0;i<m;i++){
      REAL(yval)[i]=REAL(y)[INTEGER(indx)[i]];
      REAL(wts)[i]=REAL(case_wts)[INTEGER(indx)[i]];
    }
    if(INTEGER(type)[0]==0L){
      scaleL(REAL(yval),REAL(resp_info),lfold,m,1);
    }
    
    INTEGER(_m)[0]=lfold;
    if(qt>0){
      for(i=0;i<px;i++){
        for(j=0;j<m;j++){
          v=INTEGER(indx)[j]+i*n;
          v1=j+i*n;
          REAL(X)[v1]=REAL(_X)[v];
        }
      }
    }
    for(i=0;i<pz;i++){
      for(j=0;j<m;j++){
        v=INTEGER(indx)[j]+i*n;
        v1=j+i*n;
        REAL(Z)[v1]=REAL(_Z)[v];
      }
    }
    for(i=0;i<ntunes;i++){
      for(i1=0;i1<n;i1++){
        REAL(fit)[i1]=0.0;
      }
      REAL(tvec)[0]=REAL(tunes)[i];
      REAL(tvec)[1]=REAL(tunes)[ntunes+i];
      REAL(tvec)[2]=REAL(tunes)[2*ntunes+i];
      
      
      SET_VECTOR_ELT(mod,0,PROTECT(ARIDGE(Z,X,W,yval,wts,_m,tvec,type,thresh,eps,cv_mode,_stab,fit)));
      INTEGER(convs)[m1]+=INTEGER(VECTOR_ELT(mod,0))[1];
      INTEGER(errs)[m1]+=INTEGER(VECTOR_ELT(mod,0))[0];
      
      if(INTEGER(type)[0]==0L){
        for(i1=lfold;i1<m;i1++){
          tmp=REAL(yval)[i1]-REAL(fit)[i1];
          if(tmp!=tmp){
            tmp=100.0;
          }
          REAL(cv_errs)[m1]+=tmp*tmp;
        }
      }else{
        for(i1=lfold;i1<m;i1++){
          tmp=log(1+exp(-2.0*REAL(yval)[i1]*REAL(fit)[i1]));
          if(tmp!=tmp){
            tmp=100.0;
          }
          REAL(cv_errs)[m1]+=tmp;
        }
      }
      m1++;
      
      UNPROTECT(1);
    }
  }
  SET_VECTOR_ELT(final,0,PROTECT(duplicate(cv_errs)));
  SET_VECTOR_ELT(final,1,PROTECT(duplicate(convs)));
  SET_VECTOR_ELT(final,2,PROTECT(duplicate(errs)));
  prot+=3;
  UNPROTECT(prot);
return final;
}

SEXP ARIDGE(SEXP Z,SEXP _X,SEXP _W,SEXP _y,SEXP case_wts,SEXP _m,SEXP tunes,SEXP type,SEXP thresh,SEXP eps,SEXP cv_mode,SEXP _stab,SEXP eta){
  SEXP adims1 = PROTECT(getAttrib(Z, R_DimSymbol));
  int n=INTEGER(adims1)[0], pz= INTEGER(adims1)[1],px=0,mpa;
  int m=asInteger(_m),i,prot,qt=0;
  double lam1=REAL(tunes)[0],lam2=REAL(tunes)[1],gam=REAL(tunes)[2];
  SEXP y=PROTECT(duplicate(_y));
  SEXP W=PROTECT(duplicate(_W));
  SEXP X=PROTECT(duplicate(_X));
  SEXP wts=PROTECT(duplicate(case_wts));
  if(!Rf_isNull(X)){
    qt=1;
    px=INTEGER(getAttrib(X, R_DimSymbol))[1];
  }
  SEXP varsx=PROTECT(allocVector(INTSXP,px+1-qt));
  for(i=0;i<(px+1-qt);i++){
    INTEGER(varsx)[i]=0;
  }
  SEXP resp_info=PROTECT(allocVector(REALSXP,2));
  if(qt==1){
    px=rmv_zero_var(REAL(X),m,n,px,INTEGER(varsx));
    scaleL(REAL(X),REAL(resp_info),m,n,px);
    if(px==0){
      qt=0;
    }
  }
  mpa=px;
  if(pz>px){
    mpa=pz;
  }
  SEXP bwork=PROTECT(allocVector(REALSXP,px*px+1));
  SEXP beta=PROTECT(allocVector(REALSXP,pz+px));
  SEXP awork=PROTECT(allocVector(REALSXP,mpa*mpa));
  SEXP zwork=PROTECT(allocVector(REALSXP,(mpa+n)*mpa));
  SEXP lvec=PROTECT(allocVector(VECSXP,7));
  SEXP adjust=PROTECT(allocVector(REALSXP,2));
  SEXP conv=PROTECT(allocVector(INTSXP,2));
  SEXP meas=PROTECT(allocVector(REALSXP,1));
  SEXP dims=PROTECT(allocVector(INTSXP,4));
  
  prot=16;
  REAL(adjust)[0]=0.0;
  
  for(i=0;i<LENGTH(Z);i++){
    REAL(zwork)[i]=REAL(Z)[i];
  }
  for(i=m;i<n;i++){
    REAL(wts)[i]*=gam/(REAL(wts)[i]+gam);
  }
  REAL(meas)[0]=0.0;

  if(asInteger(type)==0L){
    INTEGER(conv)[0]=agraph_regress(REAL(Z),X,REAL(W),REAL(y),REAL(wts), REAL(resp_info),REAL(awork),REAL(bwork),REAL(zwork),pz,px,n,m,lam1,lam2,asReal(_stab),REAL(beta),&REAL(beta)[pz],REAL(eta));
    
    INTEGER(conv)[1]=0;
  }else{
    INTEGER(conv)[0]=agraph_lscore(REAL(Z),X,REAL(W),REAL(y),REAL(wts),REAL(awork),REAL(bwork),REAL(zwork),pz,px,n,m,lam1,lam2,asReal(_stab),asInteger(thresh),asReal(eps),INTEGER(conv),REAL(resp_info),REAL(beta),REAL(eta));
  }
  //Scale back
  for(i=0;i<n;i++){
    REAL(eta)[i]=REAL(eta)[i]*REAL(resp_info)[1]+REAL(resp_info)[0];
  }
  for(i=0;i<(pz+px*qt);i++){
    REAL(beta)[i]*=REAL(resp_info)[1];
  }
  if(INTEGER(cv_mode)[0]==1L){
    UNPROTECT(prot);
    return conv;
  }
  
  if(INTEGER(conv)[0]==0L){
    SEXP zwork1=PROTECT(duplicate(Z));
    prot++;
   REAL(meas)[0]=get_DF(REAL(Z),X,REAL(W),REAL(wts),REAL(awork),REAL(bwork),REAL(zwork),REAL(zwork1),n, m,pz,px);

  }else{
     REAL(meas)[0]=(double)m;
   }
  INTEGER(dims)[0]=n;
  INTEGER(dims)[1]=m;
  INTEGER(dims)[2]=px*qt;
  INTEGER(dims)[3]=pz;
    
  SET_VECTOR_ELT(lvec,0,PROTECT(duplicate(beta)));
  SET_VECTOR_ELT(lvec,1,PROTECT(duplicate(varsx)));
  SET_VECTOR_ELT(lvec,2,PROTECT(duplicate(conv)));
  SET_VECTOR_ELT(lvec,3,PROTECT(duplicate(meas)));
  SET_VECTOR_ELT(lvec,4,PROTECT(duplicate(dims)));
  SET_VECTOR_ELT(lvec,5,PROTECT(duplicate(resp_info)));
  SET_VECTOR_ELT(lvec,6,PROTECT(duplicate(wts)));
  prot+=7;
  UNPROTECT(prot);
  return(lvec);
}

int agraph_lscore(double *Z,SEXP x,double *rL, double *y,double *wts,double *awork,double *bwork,double *zwork,int a,int p,int n,int m,double lam1,double lam2,double stab,int thresh,double ep,int *conv,double *adjust,double *beta,double *eta){
  int i=0,k=0,prot=0,err,a2=a*a,na=n*a,pa=p+a;
  double vec=0.0,sv=0.0;
  SEXP peta=PROTECT(allocVector(REALSXP,n));
  SEXP iwts=PROTECT(allocVector(REALSXP,n));
  SEXP etan=PROTECT(allocVector(REALSXP,n));
  SEXP z=PROTECT(allocVector(REALSXP,m));
  SEXP W=PROTECT(allocVector(REALSXP,a2));
  for(i=0;i<n;i++){
    REAL(etan)[i]=eta[i];
    if(i<m){
      REAL(z)[i]=y[i];
    }
    REAL(iwts)[i]=wts[i];
  }
  prot=5;
  do{
    for(i=0;i<n;i++){
      REAL(peta)[i]=exp(eta[i]);
      REAL(peta)[i]/=(1+REAL(peta)[i]);
      REAL(iwts)[i]=REAL(peta)[i]*(1-REAL(peta)[i])*wts[i];
      if(i<m){
        REAL(z)[i]=eta[i]+((y[i]+1.0)/2.0-REAL(peta)[i])/REAL(iwts)[i];
      }
    }
    for(i=0;i<a2;i++){
      REAL(W)[i]=rL[i];
    }
    for(i=0;i<na;i++){
      zwork[i]=Z[i];
    }
    err=agraph_regress(Z,x,REAL(W),REAL(z),REAL(iwts),adjust,awork,bwork,zwork,a,p,n,m,lam1,lam2,stab,beta,&beta[a],REAL(etan));
    if(err!=0L){
      break;
    }
    sv=0.0;
    vec=0.0;

    for(i=0;i<n;i++){
      vec+=wts[i]*(REAL(etan)[i]-eta[i])*(REAL(etan)[i]-eta[i]);
      sv+=wts[i]*eta[i]*eta[i];
      eta[i]=REAL(etan)[i];
    }
    vec/=sv;
    if(vec!=vec){
      err=-1;
      for(i=0;i<pa;i++)
        beta[i]=0.0;
      break;
    }
    if(vec<ep)
      break;
    k++;
  } while(k<thresh);
  for(i=0;i<n;i++){
    wts[i]=REAL(iwts)[i];
  }
  conv[1]=k;
  UNPROTECT(prot);
  return(err);
}
int agraph_regress(double *z,SEXP x, double *rL,double *y, double *wts,double *adjust,double *awork, double *bwork,double *zwork,
                   int a, int p, int n,int m,double lam1,double lam2,double stab,double *f1, double *f2, double *eta){
  int i,j,one=1,info=0,np=n*p;
  double alpha,beta;
  
  scaleL(y,adjust,m,m,1);
  
  dlmult(zwork,wts,n,a);
  alpha=1.0;
  beta=0.0;
  F77_CALL(dgemv)("T",&m,&a,&alpha,zwork,&n,y,&one,&beta,f1,&one);
  
  alpha=1.0;
  beta=lam1;
  F77_CALL(dgemm)("T","N",&a,&a,&n,&alpha,zwork,&n,z,&n,&beta,rL,&a);
  for(i=0;i<a;i++){
    rL[i+a*i]+=stab;
    
  }
  F77_CALL(dposv)( "L", &a, &one, rL, &a,f1, &a, &info );
  if(info!=0L){
    info=1L;
    for(i=0;i<n;i++){
      eta[i]=0.0;
    }
    return(info);
  }
  alpha=1.0;
  beta=0.0;
  F77_CALL(dgemv)("N",&n,&a,&alpha,z,&n,f1,&one,&beta,eta,&one);//eta=Hyz
  if(p==0){
    return(info);
  }
  
  for(i=0;i<n;i++){
    if(i<m){
      zwork[i]=(y[i]-eta[i])*wts[i];
    }else{
      zwork[i]=-eta[i]*wts[i];
    }
  }
  
  alpha=1.0;
  beta=0.0;
  F77_CALL(dgemv)("T",&n,&p,&alpha,REAL(x),&n,zwork,&one,&beta,f2,&one);//f2=x^TW(I-H)yz
 
  
  for(i=0;i<np;i++){
    zwork[i]=REAL(x)[i];
  }
  dlmult(zwork,wts,n,p);
  alpha=1.0;
  beta=0.0;
  F77_CALL(dgemm)("T","N",&p,&p,&n,&alpha,zwork,&n,REAL(x),&n,&beta,bwork,&p);// X^T WX
  
  F77_CALL(dgemm)("T","N",&p,&a,&n,&alpha,zwork,&n,z,&n,&beta,awork,&p);// X^T WZ
  
  for(j=0;j<p;j++){
    for(i=0;i<a;i++){
      zwork[i+j*a]=awork[j+i*p]; //Z^TWX
    }
  }
  
  F77_CALL(dtrtrs)("L","N","N", &a,&p,rL,&a,zwork,&a,&info);
  F77_CALL(dtrtrs)("L","T","N", &a,&p,rL,&a,zwork,&a,&info);
  
  alpha=-1.0;
  beta=1.0;
  F77_CALL(dgemm)("N","N",&p,&p,&a,&alpha,awork,&p,zwork,&a,&beta,bwork,&p);
  for(i=0;i<p;i++){
    bwork[i+p*i]+=lam2+stab;
  }
  F77_CALL(dposv)( "L", &p, &one,bwork, &p,f2, &p, &info );//(x^TW(I-H)x+lam2*I)^{-1}x^TW(I-H)yz
  if(info!=0L){
    info=10L;
    for(i=0;i<n;i++){
      eta[i]=0.0;
    }
    return(info);
  }
  for(i=0;i<a;i++){
    awork[i]=f1[i];
  }
  alpha=1.0;
  beta=0.0;
  F77_CALL(dgemv)("N",&a,&p,&alpha,zwork,&a,f2,&one,&beta,f1,&one);
  for(i=0;i<a;i++){
    awork[i]-=f1[i];
  }
  alpha=-1.0;
  beta=1.0;
  F77_CALL(dgemv)("N",&n,&a,&alpha,z,&n,f1,&one,&beta,eta,&one);//eta=Hx(x^TW(I-H)x+lam2*I)^{-1}x^TW(I-H)yz

  for(i=0;i<a;i++){
    f1[i]=awork[i];
  }
  alpha=1.0;
  beta=1.0;
  F77_CALL(dgemv)("N",&n,&p,&alpha,REAL(x),&n,f2,&one,&beta,eta,&one);//eta=x(x^TW(I-H)x+lam2*I)^{-1}x^TW(I-H)yz+Hyz

  return(info);
  
  
}

double get_DF(double *z, SEXP X,double *rL,double *wts,double *awork,double *bwork,double *zwork,double *zwork1,int n, int m, int a, int p){
  int info=0L,i,prot=0,j;
  double alpha,beta,df=0.0;
  
  dlmult(zwork1,wts,n,a);
  alpha=1.0;
  beta=0.0;
  F77_CALL(dgemm)("T","N",&a,&a,&m,&alpha,z,&n,zwork1,&n,&beta,awork,&a);
  F77_CALL(dtrtrs)("L","N","N", &a,&a,rL,&a,awork,&a,&info);
  F77_CALL(dtrtrs)("L","T","N", &a,&a,rL,&a,awork,&a,&info);
  if(info!=0L){
    df=0.0;
  }else{
    for(i=0;i<a;i++){
      df+=awork[i+a*i];
    }
  }
  if(p>0){
    alpha=-1.0;
    beta=0.0;
    F77_CALL(dgemm)("T","N",&p,&a,&n,&alpha,REAL(X),&n,zwork1,&n,&beta,awork,&p);
    SEXP xwork=PROTECT(allocVector(REALSXP,m*p));
    prot++;
    for(i=0;i<m;i++){
      for(j=0;j<p;j++){
        REAL(xwork)[i+j*m]=REAL(X)[i+j*n];
      }
    }
    alpha=-1.0;
    beta=1.0;
    F77_CALL(dgemm)("N","N",&m,&p,&a,&alpha,z,&n,zwork,&a,&beta,REAL(xwork),&m);
    alpha=1.0;
    F77_CALL(dtrsm)("R","L","T","N",&m,&p,&alpha,bwork,&p,REAL(xwork),&m);
    F77_CALL(dtrsm)("R","L","N","N",&m,&p,&alpha,bwork,&p,REAL(xwork),&m);
    alpha=1.0;
    beta=0.0;
    F77_CALL(dgemm)("T","N",&a,&p,&m,&alpha,zwork1,&n,REAL(xwork),&m,&beta,zwork,&a);
    F77_CALL(dtrtrs)("L","N","N", &a,&p,rL,&a,zwork,&a,&info);
    F77_CALL(dtrtrs)("L","T","N", &a,&p,rL,&a,zwork,&a,&info);
    if(info!=0L){
      df=0.0;
    }else{
      dlmult(REAL(xwork),wts,m,p);
      alpha=1.0;
      beta=0.0;
      F77_CALL(dgemm)("T","N",&p,&p,&m,&alpha,REAL(X),&n,REAL(xwork),&m,&beta,bwork,&p);
      
      alpha=1.0;
      beta=1.0;
      
      F77_CALL(dgemm)("N","N",&p,&p,&a,&alpha,awork,&p,zwork,&a,&beta,bwork,&p);
      for(i=0;i<p;i++){
        df+=bwork[i+p*i];
      }
    }
    UNPROTECT(prot);
  }
  return(df);  
}

