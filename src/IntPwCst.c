/************************************************/
/* Routine for estimating the cumulative hazard */
/* (piecewise constant hazard)                  */
/* Author: H. Charvat                           */
/* Last modified: 2016/05/10                    */
/* Part of the mexhaz package                   */
/************************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>

SEXP IntPwCst(SEXP nph, SEXP param, SEXP leint, SEXP lerem, SEXP whint)
{
  SEXP result;
  int lx = length(lerem);
  int lnph = length(nph);
  int lp = length(param);
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(param = coerceVector(param,REALSXP));
  PROTECT(leint = coerceVector(leint,REALSXP));
  PROTECT(lerem = coerceVector(lerem,REALSXP));
  PROTECT(whint = coerceVector(whint,INTSXP));
  PROTECT(result = allocVector(REALSXP,(lx*2)));
  double *Param = REAL(param);
  double *Nph = REAL(nph);
  double *Leint = REAL(leint);
  double *Lerem = REAL(lerem);
  int *Whint = INTEGER(whint);
  double *Result = REAL(result);
  int nfct = lnph/lx;
  int npar = lp/nfct;
  int i, k, z, t2;
  double Temp;
  for (z=0; z<lx; z++){
    Result[z] = 0;
    Result[z+lx] = 0;
    t2 = z*nfct;
    for (i=0; i<nfct; i++){
      Result[z] += Param[npar*i+Whint[z]]*Nph[i+t2];
    }
    Result[z+lx] += exp(Result[z])*Lerem[z];
    if (Whint[z]>0){
      for (k=Whint[z]; k>0; k--){
	Temp = 0;
	for (i=0; i<nfct; i++){
	  Temp += Param[npar*i+(k-1)]*Nph[i+t2];
	}
	Result[z+lx] += exp(Temp)*Leint[k-1];
      }
    }
    Result[z+lx] = log(Result[z+lx]);
  }
  UNPROTECT(6);
  return result;
}
