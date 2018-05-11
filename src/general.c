#include "general.h"

SEXP lgraph(SEXP _A,SEXP _m,SEXP _n,SEXP _q){
  int q=asInteger(_q);
  int n=asInteger(_n);
  int m=asInteger(_m);
  int ishift,prot;
  int u=n-m,i,j,info;
  double alpha=1.0,beta=0.0;
  
  SEXP B=PROTECT(duplicate(_A));
  SEXP r_sum = PROTECT(allocVector(REALSXP,n));
  SEXP llwork=PROTECT(allocVector(REALSXP,m*m));
  SEXP ans = PROTECT(allocVector(VECSXP, q+1));
  SEXP errs = PROTECT(allocVector(INTSXP, q));
  
  prot=5;
  for(ishift=0;ishift<q;ishift++){
    lap(REAL(VECTOR_ELT(VECTOR_ELT(B,ishift),0)),REAL(r_sum),n,0.0,1.0);
    
    F77_CALL(dposv)("L", &u, &m, &REAL(VECTOR_ELT(VECTOR_ELT(B,ishift),0))[n*m+m], &n,&REAL(VECTOR_ELT(VECTOR_ELT(B,ishift),0))[m], &n, &info);
    INTEGER(errs)[ishift]=info;
    if(info!=0){
      for(i=0;i<LENGTH(llwork);i++){
        REAL(llwork)[i]=0.0;
      }
    }else{
      F77_CALL(dgemm)("N","N", &m,&m,&u,&alpha,&REAL(VECTOR_ELT(VECTOR_ELT(B,ishift),0))[n*m],&n,&REAL(VECTOR_ELT(VECTOR_ELT(B,ishift),0))[m],&n,&beta,REAL(llwork),&m);

      for(i=0;i<m;i++){
        for(j=0;j<m;j++){
          REAL(llwork)[i+j*m]+=REAL(VECTOR_ELT(VECTOR_ELT(_A,ishift),0))[i+j*n];
        }
      }
    }

    SET_VECTOR_ELT(ans,ishift,PROTECT(duplicate(llwork)));
    prot++;
  }
  SET_VECTOR_ELT(ans,q,PROTECT(duplicate(errs)));
  prot++;
  UNPROTECT(prot);
  return ans;
}


int ginv(double *P,double *alp,double *work,int *rank,int n,int get_df){
  int i=-1L,one=1L,info,prot=0;
  SEXP kwork=PROTECT(allocVector(REALSXP,WORK_BUFFER)); //not sure how big to make this thing in general
  prot++;
  double rcond=-1;
  if(get_df==1L){
    F77_CALL(dgelss)(&n, &n, &one,P,&n,alp,&n,work,&rcond,&rank[1],REAL(kwork),&i,&info);
    i=(int)REAL(kwork)[0];
    if(i<1)i=1L;
    SEXP new_work=PROTECT(allocVector(REALSXP,i));
    prot++;
    F77_CALL(dgelss)(&n, &n, &one,P,&n,alp,&n,work,&rcond,&rank[1],REAL(new_work),&i,&info);
  }else{
    SEXP iwork=PROTECT(allocVector(INTSXP,n*n+11*n));
    prot++;
    F77_CALL(dgelsd)(&n, &n, &one,P,&n,alp,&n,work,&rcond,&rank[1],REAL(kwork),&i,INTEGER(iwork),&info);
    i=(int)REAL(kwork)[0];
    if(i<1)i=1;
    SEXP new_work=PROTECT(allocVector(REALSXP,i));
    prot++;
    F77_CALL(dgelsd)(&n, &n, &one,P,&n,alp,&n,work,&rcond,&rank[1],REAL(new_work),&i,INTEGER(iwork),&info);
  }
  UNPROTECT(prot);
  return info;
}

void lap(double* A,double *r_sum,int n,double stab,double lam){
  int i,j;
  for(j=0;j<n;j++){
    r_sum[j]=0.0;
    for(i=0;i<n;i++){
      r_sum[j]+=A[i+j*n]+stab;
    }
    A[j+j*n]=lam*(r_sum[j]-A[j+j*n]-stab);
    for(i=0;i<n;i++){
      if(i!=j)
        A[i+j*n]=-lam*A[i+j*n]-lam*stab;
    }
  }
  for(i=0;i<n;i++){
    r_sum[i]=1/sqrt(r_sum[i]);
  }
}


void dlmult(double *A,double *d,int n,int p){
  int i=0,j=0;
  for(;j<n;j++){
    for(i=0;i<p;i++)
      A[j+i*n]=d[j]*A[j+i*n];
  }
}

void drmult(double *A,double *d,int n,int p){
  int i=0,j=0;
  for(;j<p;j++){
    for(i=0;i<n;i++)
      A[i+j*n]=d[j]*A[i+j*n];
  }
}

void dlrmult(double *A,double *d,int n,int fn){
  int i=0,j=0;
  for(;j<n;j++){
    for(i=0;i<n;i++)
      A[i+j*fn]=d[j]*A[i+j*fn]*d[i];
  }
}

int rmv_zero_var(double *x, int m, int n,int p, int *vars){
  int i=0,j=0,skip=0;
  double sos=0.0,mu=0.0,tmp;
  
  for(i=0;i<p;i++){
    mu=0.0;
    for(j=0;j<m;j++){
      mu+=x[j+i*n];
    }
    mu/=(double)m;
    vars[i]=1;
    sos=0.0;
    for(j=0;j<m;j++){
      tmp=x[j+i*n]-mu;
      sos+=tmp*tmp;
    }
    if(sos>0.0){
      for(j=0;j<n;j++){
        x[j+skip*n]=x[j+i*n];
      }
      skip++;
    }else{
      vars[i]=0;
    }
  }
  return skip;
}

void scaleL(double *x,double *moments,int m, int n,int p){
  mean_center(x,moments,m,n,p);
  sos(x,moments,m,n,p);
}
void sum_vec(double *svec,double *mat,int start,int end,int p,double cst){
  int i,j;
  for(j=0;j<p;j++){
    for(i=start;i<end;i++){
      svec[j]+=mat[i+j*end];
    }
   svec[j]=svec[j]*cst;
  }
}


int colzero(double *x,int m,int n,int p,int* vars){
  int i=0,j=0,skip=0;
  double rs=0.0;
  for(i=0;i<p;i++){
    vars[i]=1;
    rs=0.0;
    for(j=0;j<m;j++){
      rs+=x[j+i*n]*x[j+i*n];
    }
    if(rs>0.0){
      for(j=0;j<n;j++){
        x[j+skip*n]=x[j+i*n];
      }
      skip++;
    }else{
      vars[i]=0;
    }
  }
  return(skip);
}


void mean_center(double* x,double *mu,int m,int n,int p){
  int i=0,j=0;
  double fac=1.0/(double)m;
  

  for(i=0;i<p;i++){
    mu[0]=0.0;
    for(j=0;j<m;j++){
      mu[0]+=x[j+i*n];
    }
    mu[0]=mu[0]*fac;
    for(j=0;j<n;j++){
      x[j+i*n]-=mu[0];
     }
  }
}
void sos(double *x,double *sos,int m, int n,int p){
  int i=0,j=0;
  for(i=0;i<p;i++){
    sos[1]=0.0;
    for(j=0;j<m;j++){
      sos[1]+=x[j+i*n]*x[j+i*n];
    }
    sos[1]=sqrt(fabs(sos[1]));
    if(sos[1]>0.0){
      for(j=0;j<n;j++){
        x[j+i*n]=x[j+i*n]/sos[1];
      }
    }
  }
}

double mean(double *v,int n){
  int i;
  double s=0.0;
  for (i = 0; i < n; i++) s += v[i];
  s /= n;
  return s;
}

void printit(double *a,int n){
  int i;
  for(i=0;i<n;i++){
    Rprintf(" %f",a[i]);
  }
  Rprintf("\n");
}
double sign(double x){
  if(x>0.0){
    return(1.0);
  }else{
    return(-1.0);
  }
}
double max(double a,double b){
  if(a<b) return(b);
  return(a);
}



void q_sort(double *a,int n){
  
  qsort((void*)a, n, sizeof(double), compare);
}
int compare(const void *a, const void *b)
{
  double x=*(double *)a;
  double y=*(double *)b;
  
  if(x>y) return -1;
  if(x<y) return 1;
  return 0;
  
}
