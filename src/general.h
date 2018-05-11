#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "R_ext/Rdynload.h"
#include <R_ext/Complex.h>
#include <R_ext/RS.h>
#include <math.h>  

#define La_extern
#define BLAS_extern

#include <R_ext/Lapack.h>

#undef La_extern
#undef BLAS_extern

#define STABILITY_THRESH             1E-15
#define WORK_BUFFER                  10
#define CV_THRESH                    1E-15

SEXP lgraph(SEXP,SEXP,SEXP,SEXP);
void lap(double*,double*,int,double,double);
void scaleL(double*,double*,int,int,int);
void mean_center(double*,double*,int,int,int);
void sum_vec(double*,double*,int,int,int,double);
void dlmult(double*,double*,int,int);
void drmult(double*,double*,int,int);
void dlrmult(double*,double*,int,int);
int rmv_zero_var(double*,int,int,int,int*);
void sos(double*,double*,int,int,int);
double mean(double*,int);
double sign(double);
void printit(double*,int);
void q_sort(double*,int);
int compare(const void *, const void *);
double max(double,double);
int colzero(double*,int,int,int,int*);
int ginv(double*,double*,double*,int*,int,int);

