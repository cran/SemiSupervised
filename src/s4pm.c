 #include"s4pm.h"

SEXP cv_s4pm_fit(SEXP folds,SEXP omits,SEXP _A,SEXP _x,SEXP y,SEXP case_wts,SEXP tunes,SEXP l_thresh,SEXP l_eps,SEXP _type){
  int i=0,j=0,k=0,qt=0,n=asInteger(VECTOR_ELT(_A,3)),K=LENGTH(folds),m=LENGTH(y),m1=0,i1=0,k3=0;
  int prot=0,lfold=0,v=0,v1=0,p=0,n2=n*n;
  int ntunes;

  double tmp,last,next,lratio;
  
  SEXP A=PROTECT(duplicate(_A));
  SEXP x=PROTECT(duplicate(_x));
  
  if(!Rf_isNull(tunes)){
    ntunes=INTEGER(getAttrib(tunes, R_DimSymbol))[0];
  }else{
    ntunes=0;
  }
  if(!Rf_isNull(x)){
    qt=1;
    p=INTEGER(getAttrib(x, R_DimSymbol))[1];
  }
  SEXP yval=PROTECT(duplicate(y));
  SEXP indx=PROTECT(allocVector(INTSXP,m));
  SEXP cv_errs=PROTECT(allocVector(REALSXP,ntunes));
  SEXP _m=PROTECT(allocVector(INTSXP,1));
  SEXP mat=PROTECT(allocVector(REALSXP,n*n));
  SEXP row_sum=PROTECT(allocVector(REALSXP,n));
  SEXP wts=PROTECT(allocVector(REALSXP,n));
  SEXP fit=PROTECT(allocVector(REALSXP,n));
  
  SEXP final=PROTECT(allocVector(VECSXP,3));
  SEXP tvec=PROTECT(allocVector(REALSXP,4));
  SEXP convs=PROTECT(allocVector(INTSXP,ntunes));
  SEXP errs=PROTECT(allocVector(INTSXP,ntunes));
  
  SEXP cvmode=PROTECT(allocVector(INTSXP,1));
  SEXP resp_info=PROTECT(allocVector(REALSXP,2));
  SEXP mod=PROTECT(allocVector(VECSXP,2));
  
  prot=17;
  
  INTEGER(cvmode)[0]=1L;
  
  last=REAL(tunes)[ntunes];
  lap(REAL(VECTOR_ELT(A,0)),REAL(row_sum),n,asReal(VECTOR_ELT(A,1)),last);
  if(asInteger(VECTOR_ELT(A,2))==1){
    dlrmult(REAL(VECTOR_ELT(A,0)),REAL(row_sum),n,n);
  }
  
  for(i=0;i<n2;i++){
    REAL(mat)[i]=REAL(VECTOR_ELT(A,0))[i];
  }
  
  for(i=0;i<n;i++){
    REAL(wts)[i]=REAL(case_wts)[i];
    REAL(fit)[i]=0.0;
  }
  for(i=0;i<LENGTH(cv_errs);i++){
    REAL(cv_errs)[i]=0.0;
    INTEGER(convs)[i]=0;
    INTEGER(errs)[i]=0;
  }
  
  for(k=0;k<K;k++){
    m1=0;
    lfold=LENGTH(VECTOR_ELT(folds,k));
    for(i=0;i<m;i++){
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
    if(INTEGER(_type)[0]==0L){
      scaleL(REAL(yval),REAL(resp_info),lfold,m,1);
    }
    
    for(i=0;i<m;i++){
      for(j=0;j<m;j++){
        v=INTEGER(indx)[j]+INTEGER(indx)[i]*n;
        v1=j+i*n;
        REAL(VECTOR_ELT(A,0))[v1]=REAL(mat)[v];
      }
    }
    for(j=0;j<m;j++){
      for(i=m;i<n;i++){
        v=INTEGER(indx)[j]+i*n;
        v1=i+j*n;
        REAL(VECTOR_ELT(A,0))[v1]=REAL(mat)[v];
        v1=j+i*n;
        REAL(VECTOR_ELT(A,0))[v1]=REAL(mat)[v];
      }
    }
    if(qt>0){
      for(i=0;i<p;i++){
        for(j=0;j<m;j++){
          v=INTEGER(indx)[j]+i*n;
          v1=j+i*n;
          REAL(x)[v1]=REAL(_x)[v];
        }
      }
    }
    last=REAL(tunes)[ntunes];
    for(i=0;i<ntunes;i++){
      next=REAL(tunes)[ntunes+i];
      lratio=next/last;
       if(fabs(lratio-1.0)>1e-6){
         for(i1=0;i1<n2;i1++){
           REAL(VECTOR_ELT(A,0))[i1]*=lratio;
         }
       }
      last=next;
      INTEGER(_m)[0]=lfold;
      for(i1=0;i1<n;i1++){
        REAL(fit)[i1]=0.0;
      }
      REAL(tvec)[0]=REAL(tunes)[i];
      REAL(tvec)[1]=next;
      REAL(tvec)[2]=REAL(tunes)[i+2*ntunes];
      REAL(tvec)[3]=REAL(tunes)[i+3*ntunes];
      SET_VECTOR_ELT(mod,0,PROTECT(s4pm_fit(A,x,yval,wts,_m,tvec,cvmode,_type,l_thresh,l_eps,fit)));
      INTEGER(convs)[m1]+=INTEGER(VECTOR_ELT(VECTOR_ELT(mod,0),1))[0];
      INTEGER(errs)[m1]+=INTEGER(VECTOR_ELT(VECTOR_ELT(mod,0),0))[0];
      k3++;
      if(INTEGER(_type)[0]==0L){
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
           //tmp=0.0;
           //if(REAL(yval)[i1]*REAL(fit)[i1]<0.0)tmp=1.0;
           REAL(cv_errs)[m1]+=tmp;
         }
         //REAL(cv_errs)[m1]/=(double)(m-lfold);//Graph Kernel h:  0.12  Lagrangians:  0.01   0  Safe-Lagrangian 0.1
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

SEXP s4pm_fit(SEXP _A,SEXP _x,SEXP _y,SEXP case_wts,SEXP _m,SEXP tunes,SEXP CV_Mode,
                   SEXP _type,SEXP l_thresh,SEXP l_eps,SEXP eta){
  int n=asInteger(VECTOR_ELT(_A,3)),m = asInteger(_m),i=0,p=1,qt=0,prot=0,pstar=0;
  int type=INTEGER(_type)[0];
  double lam=0.0;
  
  SEXP y=PROTECT(duplicate(_y));
  SEXP A = PROTECT(duplicate(_A));
  SEXP row_sum = PROTECT(allocVector(REALSXP,n));
  SEXP x = PROTECT(duplicate(_x));
  SEXP lvec=PROTECT(allocVector(VECSXP,9));
  
  SEXP wts=PROTECT(allocVector(REALSXP,n));
  SEXP resp_info=PROTECT(allocVector(REALSXP,2));
  lam=REAL(tunes)[2];
  if(!Rf_isNull(x)){
    qt=1;
    p=INTEGER(getAttrib(x, R_DimSymbol))[1];
    pstar=p;
  }
  SEXP vars=PROTECT(allocVector(INTSXP,p));
  for(i=0;i<p;i++){
    INTEGER(vars)[i]=0;
  }
  if(qt>0){
    p=rmv_zero_var(REAL(x),m,n,p,INTEGER(vars));
    scaleL(REAL(x),REAL(resp_info),m,n,p);
    if(p==0){
      qt=0;
    }
  }
  SEXP fn=PROTECT(allocVector(REALSXP,n+n*qt));
  SEXP bet=PROTECT(allocVector(REALSXP,p));
  SEXP ppwork=PROTECT(allocVector(REALSXP,p*p));
  SEXP npwork=PROTECT(duplicate(x));
  SEXP dP=PROTECT(allocVector(REALSXP,n*p));
  SEXP conv=PROTECT(allocVector(INTSXP,1));
  SEXP err=PROTECT(allocVector(INTSXP,1));
  SEXP df=PROTECT(allocVector(REALSXP,1));
  SEXP cv_out=PROTECT(allocVector(VECSXP,3));
  SEXP dims=PROTECT(allocVector(INTSXP,3));
  
  prot=18;
  for(i=0;i<LENGTH(fn);i++){
    REAL(fn)[i]=0.0;
  }
  for(i=0;i<LENGTH(bet);i++){
    REAL(bet)[i]=0.0;
  }
  for(i=0;i<LENGTH(ppwork);i++){
    REAL(ppwork)[i]=0.0;
  }
  for(i=0;i<m;i++){
    REAL(wts)[i]=REAL(case_wts)[i];
  }
  for(i=m;i<n;i++){
    REAL(wts)[i]=REAL(case_wts)[i]-REAL(case_wts)[i]*REAL(case_wts)[i]/(REAL(case_wts)[i]+REAL(tunes)[3]);
  }
  
   if(asInteger(CV_Mode)!=1L){
     lap(REAL(VECTOR_ELT(A,0)),REAL(row_sum),n,asReal(VECTOR_ELT(A,1)),REAL(tunes)[1]);
     if(asInteger(VECTOR_ELT(A,2))==1){
      dlrmult(REAL(VECTOR_ELT(A,0)),REAL(row_sum),n,n);
     }
   }
  
  for(i=0;i<n;i++){
    REAL(dP)[i]=REAL(VECTOR_ELT(A,0))[i+n*i];
  }
  
  INTEGER(conv)[0]=0;
   if(type==0L){
    INTEGER(err)[0]=regress(REAL(VECTOR_ELT(A,0)),x,REAL(y),REAL(wts),REAL(eta),REAL(fn),REAL(bet),
                            REAL(ppwork),REAL(resp_info),npwork,lam,REAL(tunes)[3],m,n,qt,p);
  }else{
    INTEGER(err)[0]=lscore(REAL(VECTOR_ELT(A,0)),x,REAL(y),REAL(wts),REAL(eta),REAL(fn),REAL(bet),
                           REAL(ppwork),REAL(dP),npwork,INTEGER(conv),REAL(resp_info),lam,REAL(tunes)[3],REAL(l_eps)[0],
                           INTEGER(l_thresh)[0],m,n,qt,p);
  }
  
  //Scale back
  for(i=0;i<n;i++){
    REAL(eta)[i]=REAL(eta)[i]*REAL(resp_info)[1]+REAL(resp_info)[0];
    REAL(fn)[i]=REAL(fn)[i]*REAL(resp_info)[1];
  }
  if(qt>0){
    for(i=0;i<p;i++){
      REAL(bet)[i]*=REAL(resp_info)[1];
    }
    for(i=0;i<n;i++){
      REAL(fn)[i+n]=REAL(fn)[i+n]*REAL(resp_info)[1];
    }
  }
  SEXP tmp=PROTECT(allocVector(REALSXP,2));
  prot++;
  mean_center(REAL(fn),REAL(tmp),m,n,1);

  if(asInteger(CV_Mode)>0L){
    SET_VECTOR_ELT(cv_out,0,PROTECT(duplicate(err)));
    SET_VECTOR_ELT(cv_out,1,PROTECT(duplicate(conv)));
    prot+=2;
    UNPROTECT(prot);
    return cv_out;
  }
  if(INTEGER(err)[0]==0L){
    if(REAL(tunes)[3]>0.0){
      REAL(df)[0]=get_df_s4pm(REAL(VECTOR_ELT(A,0)),REAL(wts),npwork,REAL(ppwork),m,n,n,p*qt);
    }else{
      REAL(df)[0]=get_df_s4pm(REAL(VECTOR_ELT(A,0)),REAL(wts),npwork,REAL(ppwork),m,m,n,p*qt);
    }
  }else{
    REAL(df)[0]=0.0;
  }
  INTEGER(dims)[0]=n;
  INTEGER(dims)[1]=m;
  INTEGER(dims)[2]=p*qt;

  
  SET_VECTOR_ELT(lvec,0,PROTECT(duplicate(fn)));
  SET_VECTOR_ELT(lvec,1,PROTECT(duplicate(wts)));
  SET_VECTOR_ELT(lvec,2,PROTECT(duplicate(bet)));
  SET_VECTOR_ELT(lvec,3,PROTECT(duplicate(conv)));
  SET_VECTOR_ELT(lvec,4,PROTECT(duplicate(err)));
  SET_VECTOR_ELT(lvec,5,PROTECT(duplicate(df)));
  SET_VECTOR_ELT(lvec,6,PROTECT(duplicate(dims)));
  SET_VECTOR_ELT(lvec,7,PROTECT(duplicate(resp_info)));
  SET_VECTOR_ELT(lvec,8,PROTECT(duplicate(vars)));
  prot+=9;

  UNPROTECT(prot);
  return lvec;
}

int lscore(double* M,SEXP x,double *y,double *wts,double* eta,double *fn,double *bet,
           double *ppwork,double* dP,SEXP npwork,int *conv,double *adjust,double lam,double gam,
           double ep,int thresh,int m, int n,int qt, int p){
  int i=0,j=0,k=0,prot=0,err,np=n*p;
  double vec=0.0,sv=0.0;
  
  SEXP peta=PROTECT(allocVector(REALSXP,n));
  SEXP iwts=PROTECT(allocVector(REALSXP,n));
  SEXP etan=PROTECT(allocVector(REALSXP,n));
  SEXP z=PROTECT(allocVector(REALSXP,m));
  
  prot=4;
  do{
    for(i=0;i<n;i++){
      REAL(peta)[i]=exp(eta[i]);
      REAL(peta)[i]/=(1+REAL(peta)[i]);
      REAL(iwts)[i]=REAL(peta)[i]*(1-REAL(peta)[i])*wts[i];
      if(i<m){
        REAL(z)[i]=eta[i]+((y[i]+1.0)/2.0-REAL(peta)[i])/REAL(iwts)[i];
      }
    }
    if(k>0){
      for(i=0;i<n;i++){
        for(j=i+1;j<n;j++){
          M[j+n*i]=M[n*j+i];
        }
      }
      for(i=0;i<n;i++){
        M[i+n*i]=dP[i];
      }
      if(qt>0){
        for(i=0;i<np;i++){
          REAL(npwork)[i]=REAL(x)[i];
        }
      }
    }
    err=regress(M,x,REAL(z),REAL(iwts),REAL(etan),fn,bet,ppwork,adjust,npwork,lam,gam,m,n,qt,p);
    if(err>0){
      break;
    }
    sv=0.0;
    vec=0.0;
    for(i=0;i<n;i++){
      vec+=(REAL(etan)[i]-eta[i])*(REAL(etan)[i]-eta[i]);
      sv+=eta[i]*eta[i];
      eta[i]=REAL(etan)[i];
    }
    vec/=sv;
    if(vec<ep)
      break;
    k++;
  } while(k<thresh);
  for(i=0;i<n;i++){
    wts[i]=REAL(iwts)[i];
  }
  
  conv[0]=k;
  UNPROTECT(prot);
  return(err);
}
    
int regress(double* M,SEXP x,double *y,double *wts,double* eta,double *fn,double *bet,
            double *ppwork,double *adjust,SEXP npwork,double lam,double gam,int m, int n,int qt, int p){
  int i=0,one=1,info=0,err=0,err1=0,type=qt;
  double alpha=1.0,beta=0.0;
  
  scaleL(y,adjust,m,m,1);
  
  
  for(i=0;i<n;i++){
    M[i+i*n]+=wts[i];
  }
  
  switch(type){
    case 0:
    for(i=0;i<m;i++){
      eta[i]=y[i]*wts[i];
    }
    for(i=m;i<n;i++){
      eta[i]=0.0;
    }
      alpha=1.0;
    beta=0.0;
    F77_CALL(dposv)("L", &n, &one, M, &n,eta, &n, &info);
    if(info!=0L){
      if(info<0L){
        err=1L;
        }
      if(info>0L){
        err=2L;
      }
      type=3L;
      }else{
        for(i=0;i<n;i++){
          fn[i]=eta[i];
        }
      }
    break;
    case 1:
    if(gam>0){
      err1=get_coef(bet,M,REAL(x),y,wts,ppwork,REAL(npwork),lam,m,n,n,p);
    }else{
      err1=get_coef(bet,M,REAL(x),y,wts,ppwork,REAL(npwork),lam,m,m,n,p);
    }
    err+=err1;
    if(err1>0L){
      type=3L;
    }else{
      F77_CALL(dgemv)("N", &n,&p,&alpha,REAL(x),&n,bet,&one,&beta,fn,&one);
      for(i=0;i<m;i++){
        eta[i]=wts[i]*(y[i]-fn[i]);
        fn[i+n]=fn[i];
      }
      for(i=m;i<n;i++){
        eta[i]=-wts[i]*fn[i];
        fn[i+n]=fn[i];
      }
      F77_CALL(dtrtrs)("L","N","N", &n,&one,M,&n,eta,&n,&info);
      F77_CALL(dtrtrs)("L","T","N", &n,&one,M,&n,eta,&n,&info);
      for(i=0;i<n;i++){
        fn[i]=eta[i];
        eta[i]+=fn[i+n];
      }
    }
    break;
    default:
    type=3L;
    break;
  }
  if(type==3L){
    for(i=0;i<n;i++){
      eta[i]=0.0;
      fn[i]=0.0;
      if(qt>1){
        fn[i+n]=0.0;
      }
    }
    if(qt>0){
      for(i=0;i<p;i++){
        bet[i]=0.0;
      }
    }
  }
  
  return err;
}

int get_coef(double *bet,double *M,double *x,double *y,double *wts,double *Vpp,double *Z,
             double lam,int m,int n,int fn,int p){
  int i,one=1,err=0L,np=fn*p;
  double alpha=1.0,beta=0.0;
  
  dlmult(Z,wts,fn,p);
  F77_CALL(dposv)("L", &fn, &p, M, &fn,Z, &fn, &i);
  if(i!=0L){
    if(i<0L){
      err=1L;
    }
    if(i>0L){
      err=2L;
    }
    for(i=0;i<p;i++){
      bet[i]=0.0;
    }
  }else{
    for(i=0;i<np;i++){
      Z[i]=x[i]-Z[i];
    }
    dlmult(Z,wts,fn,p);
    
    alpha=1.0;
    beta=0.0;
    F77_CALL(dgemv)("T", &m,&p,&alpha,Z,&fn,y,&one,&beta,bet,&one);

    
    alpha=1.0;
    beta=0.0;
    F77_CALL(dgemm)("T","N", &p,&p,&n,&alpha,Z,&fn,x,&fn,&beta,Vpp,&p);
    for(i=0;i<p;i++){
      Vpp[i+p*i]+=lam;
    }
    
    
    F77_CALL(dposv)("L", &p, &one, Vpp, &p,bet, &p, &i);
    if(i!=0L){
      if(i<0L){
        err=10L;
      }
      if(i>0L){
        err=20L;
      }
      for(i=0;i<p;i++){
        bet[i]=0.0;
      }
    }
  }
  return err;
}
double get_df_s4pm(double *M,double* wts,SEXP x,double *Vpp,int m,int n,int fn,int p){
  int i,j,info;
  double tr=0.0,alpha,beta;
  
  F77_CALL(dtrtri)("L","N", &fn,M,&fn,&info);
  for(i=0;i<m;i++){
    for(j=i;j<fn;j++){
      tr+=M[j+fn*i]*M[j+fn*i]*wts[i];
    }
  }
  if(p>0){
    SEXP Q=PROTECT(duplicate(x));
    SEXP A = PROTECT(allocVector(REALSXP,p*p));
    for(i=0;i<LENGTH(A);i++){
      REAL(A)[i]=0.0;
    }
    for(i=0;i<n;i++){
      wts[i]=1.0/wts[i];
    }
    dlmult(REAL(Q),wts,fn,p);
    alpha=1.0;
    beta=0.0;
    F77_CALL(dgemm)("T","N", &p,&p,&m,&alpha,REAL(x),&fn,REAL(Q),&fn,&beta,REAL(A),&p);
    for(i=0;i<n;i++){
      wts[i]=1.0/wts[i];
    }
    F77_CALL(dpotri)("L",&p,Vpp,&p,&info);
    if(info==0L){
      for(j=0;j<p;j++){
        for(i=0;i<j;i++){
          tr+=REAL(A)[i*p+j]*Vpp[i*p+j];
          tr+=REAL(A)[j*p+i]*Vpp[i*p+j];
        }
        tr+=REAL(A)[j+p*j]*Vpp[j+p*j];
      }
    }
    UNPROTECT(2);
  }
  return tr;
}

