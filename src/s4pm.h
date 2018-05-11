#include "general.h"

SEXP s4pm_fit(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP cv_s4pm_fit(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
int regress(double*,SEXP,double*,double*,double*,double*,double*,
            double*,double*,SEXP,double,double,int,int,int,int);
int lscore(double*,SEXP,double*,double*,double*,double*,double*,
           double*,double*,SEXP,int*,double*,double,double,double,int,int,int,int,int);
int get_coef(double*,double*,double*,double*,double*,double*,double*,
              double,int,int,int,int);
double get_df_s4pm(double*,double*,SEXP,double*,int,int,int,int);
