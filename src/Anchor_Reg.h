#include"general.h"

SEXP AREG_CV(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP ARIDGE(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
int agraph_lscore(double*,SEXP,double*,double*,double*,double*,double*,double*,int,int,int,int,double,double,double,int,double,int*,double*,double*,double*);
double get_DF(double*,SEXP,double*,double*,double*,double*,double*,double*,int,int,int,int);
int agraph_regress(double*,SEXP,double*,double*,double*,double*, double*,double*,double*,int,int,int,int,double,double,double,double*,double*,double*);
