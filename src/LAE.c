#include"LAE.h"

SEXP LAE(SEXP x,SEXP anchor,SEXP pos,SEXP cn, SEXP thresh, SEXP eps ){
  SEXP adims = PROTECT(getAttrib(x, R_DimSymbol));
  if (TYPEOF(adims) != INTSXP) error("non-integer dims");
  int n = INTEGER(adims)[1], p = INTEGER(adims)[0];
  int n_u=INTEGER(getAttrib(anchor, R_DimSymbol))[0];
  int p_pos=INTEGER(getAttrib(pos, R_DimSymbol))[1];
  
  SEXP U=PROTECT(allocVector(REALSXP,p_pos*n_u));
  SEXP xvec=PROTECT(allocVector(REALSXP,p));
  SEXP z=PROTECT(allocVector(REALSXP,p_pos*n));
  int i,j,k,indx,k3,prot=4;

   for(i=0;i<n;i++){
    k3=0;
    for(j=0;j<p;j++){
      REAL(xvec)[j]=REAL(x)[j+i*p];
    }
    for(k=0;k<p_pos;k++){
      indx=INTEGER(pos)[i+k*n]-1;
      for(j=0;j<n_u;j++){
        REAL(U)[k3]=REAL(anchor)[j+n_u*indx];
        k3++;
      }
    }
    LAE_FIT(REAL(xvec),REAL(U),&REAL(z)[p_pos*i],REAL(PROTECT(duplicate(xvec))),asInteger(cn),
            asInteger(thresh),asReal(eps),n_u,p_pos,i);
    UNPROTECT(1);
  }
  UNPROTECT(prot);
  return(z);
}

void LAE_FIT(double* x,double* U,double* z1,double* dif,int cn,int thresh,double eps,int n,int p_u,int s12){
  int prot=1,i,j,one=1,k;
  SEXP z0=PROTECT(allocVector(REALSXP,p_u));
  SEXP v=PROTECT(allocVector(REALSXP,p_u));
  SEXP dgv=PROTECT(allocVector(REALSXP,p_u));
  SEXP work=PROTECT(allocVector(REALSXP,p_u));
  SEXP swork=PROTECT(allocVector(REALSXP,p_u));
 
  prot=5;
  double si=1.0/(double)p_u,dt0=0.0,dt1=1.0,alpha;
  double gv,gz,tmp,tmp7,gvz,tmp2,check;
  double alp=1.0,bet=0.0;
  double twos=1.0,b,bt0=1.0,bt1=0.0;
  
  for(i=0;i<n;i++){
    dif[i]=x[i];
  }
  for(i=0;i<p_u;i++){
    REAL(v)[i]=0.0;
    REAL(z0)[i]=si;
    z1[i]=si;
    REAL(dgv)[i]=0.0;
    REAL(work)[i]=0.0;
    REAL(swork)[i]=0.0;
  }
  for(i=0;i<cn;i++){
    alpha=(dt0-1.0)/dt1;
      for(j=0;j<p_u;j++){
      REAL(v)[j]=z1[j]+alpha*(z1[j]-REAL(z0)[j]);
    }
    alp=-1.0;
    bet=1.0;
    F77_CALL(dgemv)("N",&n,&p_u,&alp,U,&n,REAL(v),&one,&bet,dif,&one);
    gv=0.0;
    for(j=0;j<n;j++){
      gv+=dif[j]*dif[j];
    }
    gv*=0.5;
    bet=0.0;
    F77_CALL(dgemv)("T",&n,&p_u,&alp,U,&n,dif,&one,&bet,REAL(dgv),&one);
    twos=0.5;
    for(k=0;k<thresh;k++){
      twos=2*twos;
      b=twos*bt0;
      tmp7=1.0/b;
      for(j=0;j<p_u;j++){
        REAL(work)[j]=REAL(v)[j]-REAL(dgv)[j]*tmp7;
        REAL(swork)[j]=REAL(work)[j];
      }
      simplex_ov(REAL(work),REAL(swork),p_u);
      alp=-1.0;
      bet=0.0;
      F77_CALL(dgemv)("N",&n,&p_u,&alp,U,&n,REAL(work),&one,&bet,dif,&one);
      gz=0.0;
      for(j=0;j<n;j++){
        dif[j]+=x[j];
        gz+=dif[j]*dif[j];
      }
      gz*=0.5;
      gvz=gv;
      tmp=0.0;
      for(j=0;j<p_u;j++){
        tmp2=REAL(work)[j]-REAL(v)[j];
        gvz+=tmp2*REAL(dgv)[j];
        tmp+=tmp2*tmp2;
      }
      tmp*=0.5*b;
      gvz+=tmp;
      if(gz <= gvz){
        bt0=bt1;
        bt1=b;
        for(j=0;j<p_u;j++){
          REAL(z0)[j]=z1[j];
          z1[j]=REAL(work)[j];
        }
        break;
      }
    }
    
    if(bt1 == 0L){
      bt0=bt1;
      bt1=b;
      for(j=0;j<p_u;j++){
        REAL(z0)[j]=z1[j];
        z1[j]=REAL(work)[j];
      }
    }
    
    dt0=dt1;
    bt0=bt1;
    dt1=0.5*(1.0+sqrt(1.0+4.0*dt1*dt1));
    check=0.0;
    for(j=0;j<p_u;j++){
      check+=fabs(REAL(z0)[j]-z1[j]);
    }
    if(check<eps){
        break;
    }
    for(j=0;j<n;j++){
      dif[j]=x[j];
    }
  }
  UNPROTECT(prot);
}
SEXP simplex_check(SEXP _x,SEXP _n,SEXP _p){
  SEXP x=PROTECT(duplicate(_x));
  SEXP sx=PROTECT(duplicate(_x));
  int n=asInteger(_n);
  int p=asInteger(_p);
  simplex(REAL(x),REAL(sx),n,p);
  
  UNPROTECT(2);
  return(x);
  
}

void order_to_s(double *x,double *z, int *ord,int n,int s){
  int i,k;
  double mnvec;
  
  
  for(i=0;i<n;i++){
    z[i]=x[i];
  }
  
  q_sort(z,n);
  
  mnvec=z[n-s];
 
  k=0;
  for(i=0;i<n;i++){
    if(x[i]<=mnvec){
      ord[k]=i;
      k++;
    }
  }
}


void simplex( double *x,double *work,int n, int p){
  int i;
  
  for(i=0;i<p;i++){
    simplex_ov(&x[i*n],&work[i*n],n);
  }
}
void simplex_ov( double *x, double *sx,int n){
  double s1=0.0,v=0.0,b=0.0;
  int j=1,kk=0;

  q_sort(sx,n);
  for(j=0;j<n;j++){
    b=sx[j];
    v=b*(double)(j)-s1;
    s1+=b;
    if(!(v>=-1.0)){
      kk=j;
      break;
    }
  }
  if(kk==0){
    kk=n;
    s1=(s1-1.0)/(double)kk;
  }else{
    s1=(s1-b-1.0)/(double)kk;
  }
  for(j=0;j<n;j++){
    x[j]=max(x[j]-s1,0.0);
  }
}
