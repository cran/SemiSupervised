#include"joint_harmonic.h"

SEXP cv_jtharm_fit(SEXP folds,SEXP omits,SEXP _A,SEXP y,SEXP case_weights,SEXP tunes,SEXP _stab,SEXP _type){
  int i=0,j=0,k=0,K=LENGTH(folds),n=LENGTH(case_weights),m=LENGTH(y),i1=0,ntunes;
  int prot=0,lfold=0,v=0,v1=0,info,one=1;
  double tmp=0.0,last,next,gam;
  double alpha=1.0,beta=0.0;
  
  if(!Rf_isNull(tunes)){
    ntunes=INTEGER(getAttrib(tunes, R_DimSymbol))[0];
  }else{
    ntunes=0;
  }
 
  
  SEXP S=PROTECT(duplicate(VECTOR_ELT(_A,0)));
  SEXP P=PROTECT(duplicate(S));
  SEXP PS=PROTECT(duplicate(S));
  SEXP ref=PROTECT(duplicate(S));
  SEXP mat=PROTECT(duplicate(S));
  SEXP _S=PROTECT(duplicate(S));
  
  SEXP yval=PROTECT(duplicate(y));
  
  SEXP row_sum=PROTECT(allocVector(REALSXP,n));
  SEXP yyu=PROTECT(allocVector(REALSXP,n));
  SEXP eta=PROTECT(allocVector(REALSXP,n));
  
  SEXP final=PROTECT(allocVector(VECSXP,2));
  SEXP cv_errs=PROTECT(allocVector(REALSXP,ntunes));
  
  SEXP nnwork=PROTECT(allocVector(REALSXP,n*n));
  SEXP indx=PROTECT(allocVector(INTSXP,n));
  SEXP ip=PROTECT(allocVector(INTSXP,n));
  
  SEXP resp_info=PROTECT(allocVector(REALSXP,2));
  SEXP wts=PROTECT(duplicate(case_weights));
  SEXP errs=PROTECT(allocVector(INTSXP,ntunes));
 
  prot=18;
 
  for(i=0;i<LENGTH(cv_errs);i++){
    REAL(cv_errs)[i]=0.0;
    INTEGER(errs)[i]=0L;
  }
  for(i=0;i<LENGTH(yyu);i++){
    REAL(yyu)[i]=0.0;
  }
  for(i=0;i<LENGTH(indx);i++){
    INTEGER(indx)[i]=i;
  }
  for(i=0;i<LENGTH(nnwork);i++){
    REAL(nnwork)[i]=0.0;
  }
  for(i=0;i<LENGTH(wts);i++){
    REAL(wts)[i]=sqrt(REAL(wts)[i]);
  }
  lap(REAL(mat),REAL(row_sum),n,asReal(VECTOR_ELT(_A,1)),1.0);
  if(asInteger(VECTOR_ELT(_A,2))==1){
    dlrmult(REAL(mat),REAL(row_sum),n,n);
  }
  
  for(i1=0;i1<ntunes;i1++){
    next=REAL(tunes)[ntunes+i1];
    gam=REAL(tunes)[ntunes*2+i1];
    if(fabs(next-last)>CV_THRESH){
      for(i=0;i<LENGTH(P);i++){
        REAL(P)[i]=REAL(mat)[i]*next;
      }
       for(i=0;i<LENGTH(S);i++){
        REAL(_S)[i]=REAL(VECTOR_ELT(_A,0))[i];
      }
      dlrmult(REAL(_S),REAL(wts),n,n);
      for(i=0;i<LENGTH(PS);i++){
        REAL(PS)[i]=REAL(_S)[i]+REAL(P)[i];
      }
      for(i=0;i<n;i++){
        REAL(PS)[i+n*i]+=asReal(_stab);
      }
      
      F77_CALL(dgesv)(&n, &n,REAL(PS),&n,INTEGER(ip),REAL(_S),&n,&info);
      if(info!=0L){
        info=1L;
        REAL(cv_errs)[i1]=100000000000.0;
        REAL(errs)[i1]=info;
      }else{
        alpha=1.0;
        beta=0.0;
        
        F77_CALL(dgemm)("N","N", &n,&n,&n,&alpha,REAL(P),&n,REAL(_S),&n,&beta,REAL(ref),&n);
        for(i=0;i<LENGTH(S);i++){
          REAL(S)[i]=REAL(_S)[i];
        }
      }
      last=next;
    }
  
    if(info!=0L){
      REAL(cv_errs)[i1]=100000000000.0;
      INTEGER(errs)[i1]=info;
      if(info!=0L)continue;
    }
    for(k=0;k<K;k++){
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
      }
      scaleL(REAL(yval),REAL(resp_info),lfold,LENGTH(yval),1);

      for(i=0;i<m;i++){
        for(j=0;j<m;j++){
          v=INTEGER(indx)[j]+INTEGER(indx)[i]*n;
          v1=j+i*n;
          REAL(PS)[v1]=REAL(ref)[v];
          REAL(S)[v1]=REAL(_S)[v];
        }
      }
     
      for(i=m;i<n;i++){
        for(j=0;j<m;j++){
          v=INTEGER(indx)[j]*n+i;
          v1=j*n+i;
          REAL(PS)[v1]=REAL(ref)[v];
          REAL(S)[v1]=REAL(_S)[v];
        }
      }
      
      for(j=0;j<m;j++){
        for(i=m;i<n;i++){
          v=INTEGER(indx)[j]+i*n;
          v1=j+i*n;
          REAL(PS)[v1]=REAL(ref)[v];
          REAL(S)[v1]=REAL(_S)[v];
        }
      }
   
      for(i=m;i<n;i++){
        for(j=m;j<n;j++){
          v=i+j*n;
          REAL(PS)[v]=REAL(ref)[v];
          REAL(S)[v]=REAL(_S)[v];
        }
      }
      
      for(i=lfold;i<n;i++){
        REAL(PS)[i+n*i]+=gam;
      }
      
      
      for(i=0;i<lfold;i++){
        REAL(yyu)[i]=REAL(yval)[i];
      }
      INTEGER(errs)[i1]+=get_yu(REAL(yyu),REAL(PS),REAL(yval),lfold,n,INTEGER(ip));
      if(INTEGER(errs)[i1]!=0L){
        REAL(cv_errs)[i1]=100000000000.0;
      }else{
        alpha=1.0;
        beta=0.0;
        F77_CALL(dgemv)("N", &n,&n,&alpha,REAL(S),&n,REAL(yyu),&one,&beta,REAL(eta),&one);
        
        if(INTEGER(_type)[0]==0L){
          for(i=lfold;i<m;i++){
            tmp=REAL(yval)[i]-REAL(eta)[i];
            if(tmp!=tmp){
              tmp=100.0;
            }
            REAL(cv_errs)[i1]+=tmp*tmp;
          }
        }else{
          for(i=lfold;i<m;i++){
            tmp=log(1+exp(-2.0*REAL(yval)[i]*REAL(eta)[i]));
            if(tmp!=tmp){
              tmp=100.0;
            }
            REAL(cv_errs)[i1]+=tmp;
          }
        }
      }
       
    }
  }
  SET_VECTOR_ELT(final,0,PROTECT(duplicate(cv_errs)));
  SET_VECTOR_ELT(final,1,PROTECT(duplicate(errs)));
  prot+=2;
  
  UNPROTECT(prot);
  return final;
}


SEXP jt_harm_fit(SEXP _A,SEXP _y,SEXP case_weights,SEXP lam,SEXP gam,SEXP _stab,SEXP f_only,SEXP f){
  int n = LENGTH(case_weights),m=LENGTH(_y);
  int prot=0,info,i,one=1;
  double alpha=1.0,beta=0.0;

  SEXP y=PROTECT(duplicate(_y));
  SEXP S = PROTECT(duplicate(VECTOR_ELT(_A,0)));
  SEXP P = PROTECT(duplicate(S));
  SEXP PS = PROTECT(duplicate(S));
  SEXP row_sum = PROTECT(allocVector(REALSXP,n));
  SEXP yyu=PROTECT(allocVector(REALSXP,n));
  SEXP lvec = PROTECT(allocVector(VECSXP, 5));
  SEXP ip=PROTECT(allocVector(INTSXP,n));
  SEXP df=PROTECT(allocVector(REALSXP,1));
  SEXP dims=PROTECT(allocVector(INTSXP,2));
  SEXP errs=PROTECT(allocVector(INTSXP,1));
  SEXP wts=PROTECT(duplicate(case_weights));
  SEXP resp_info=PROTECT(allocVector(REALSXP,2));
  prot=13;
  INTEGER(errs)[0]=0L;
  
  scaleL(REAL(y),REAL(resp_info),m,m,1);
  
  for(i=0;i<m;i++){
    REAL(yyu)[i]=REAL(y)[i];
  }
  for(i=m;i<n;i++){
    REAL(yyu)[i]=0.0;
  }
  for(i=0;i<LENGTH(wts);i++){
    REAL(wts)[i]=sqrt(REAL(wts)[i]);
  }
  lap(REAL(P),REAL(row_sum),n,asReal(VECTOR_ELT(_A,1)),REAL(lam)[0]);
  if(asInteger(VECTOR_ELT(_A,2))==1){
    dlrmult(REAL(P),REAL(row_sum),n,n);
  }
  
  dlrmult(REAL(S),REAL(wts),n,n);
  for(i=0;i<LENGTH(PS);i++){
    REAL(PS)[i]=REAL(S)[i]+REAL(P)[i];
  }
  for(i=0;i<n;i++){
    REAL(PS)[i+n*i]+=asReal(_stab);
  }
  
  
  F77_CALL(dgesv)(&n, &n,REAL(PS),&n,INTEGER(ip),REAL(S),&n,&info);
  if(info!=0){
    for(i=0;i<n;i++){
      REAL(f)[i]=0.0;
      if(i>m){
        REAL(yyu)[i]=0.0;
      }
    }
    INTEGER(errs)[0]=1L;
  }else{
    F77_CALL(dgemm)("N","N", &n,&n,&n,&alpha,REAL(P),&n,REAL(S),&n,&beta,REAL(PS),&n);
    if(m<n){
      for(i=m;i<n;i++){
        REAL(PS)[i+n*i]+=REAL(gam)[0];
      }
     INTEGER(errs)[0]=get_yu(REAL(yyu),REAL(PS),REAL(y),m,n,INTEGER(ip));
    }
    alpha=1.0;
    F77_CALL(dgemv)("N", &n,&n,&alpha,REAL(S),&n,REAL(yyu),&one,&beta,REAL(f),&one);
  }
  if(INTEGER(f_only)[0]==1L){
    UNPROTECT(prot);
    return yyu;
  }
  for(i=0;i<LENGTH(yyu);i++){
    REAL(yyu)[i]=REAL(yyu)[i]*REAL(resp_info)[1]+REAL(resp_info)[0];
    REAL(f)[i]=REAL(f)[i]*REAL(resp_info)[1]+REAL(resp_info)[0];
  }
  
  INTEGER(dims)[0]=n;
  INTEGER(dims)[1]=m;
  if(info!=0){
    REAL(df)[0]=m;
  }else{
    REAL(df)[0]=get_df_jt(REAL(PS),REAL(S),m,n,INTEGER(ip));
  }
  SET_VECTOR_ELT(lvec,0,PROTECT(duplicate(yyu)));
  SET_VECTOR_ELT(lvec,1,PROTECT(duplicate(df)));
  SET_VECTOR_ELT(lvec,2,PROTECT(duplicate(dims)));
  SET_VECTOR_ELT(lvec,3,PROTECT(duplicate(errs)));
  SET_VECTOR_ELT(lvec,4,PROTECT(duplicate(resp_info)));
  prot+=5;
  UNPROTECT(prot);
  return lvec;
}
int get_yu(double *yu,double *PS,double *y,int m, int n, int *ip){
  int i,info,one=1,u=n-m;
  double alpha=-1.0,beta=0.0;
  F77_CALL(dgemv)("N", &u,&m,&alpha,&PS[m],&n,y,&one,&beta,&yu[m],&one);
  F77_CALL(dgesv)(&u, &one,&PS[m*n+m],&n,ip,&yu[m],&n,&info);
  if(info!=0L){
    for(i=0;i<u;i++){
      (&yu[m])[i]=0.0;
    }
    info=10L;
  }
  return info;
}
double get_df_jt(double *PS,double *S,int m,int n,int *ip){
  int i,j,info,u=n-m;
  double df=0.0;
  
  for(i=0;i<m;i++){
    df+=S[i+n*i];
  }
  if(m<n){
    F77_CALL(dgetrs)("N",&u,&m,&PS[m*n+m],&n,ip,&PS[m],&n,&info);
    if(info==0L){
      for(i=m;i<n;i++){
        for(j=0;j<m;j++){
          df-=S[i*n+j]*PS[i+j*n];
        }
      }
    }
  }
  return df;
}

