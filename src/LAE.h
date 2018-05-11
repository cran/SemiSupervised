#include"general.h"
SEXP LAE(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP simplex_check(SEXP,SEXP,SEXP);

void order_to_s(double*,double*,int*,int,int);
void LAE_FIT(double*,double*,double*,double*,int,int,double,int,int,int);
void simplex( double *,double*, int,int);
void simplex_ov( double *,double *,int);
