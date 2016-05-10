/************************************************/
/* Routine for estimating the cumulative hazard */
/* (log-hazard described by linear B-spline)    */
/* Author: H. Charvat                           */
/* Last modified: 2016/05/10                    */
/* Part of the mexhaz package                   */
/************************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>

SEXP IntBs1(SEXP x, SEXP nph, SEXP param, SEXP leint, SEXP whint, SEXP knots)
{
  SEXP result;
  int lx = length(x);
  int lnph = length(nph);
  int lk = length(knots);
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(param = coerceVector(param,REALSXP));
  PROTECT(leint = coerceVector(leint,REALSXP));
  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(whint = coerceVector(whint,INTSXP));
  PROTECT(knots = coerceVector(knots,REALSXP));
  PROTECT(result = allocVector(REALSXP,(lx*2)));
  double *Param = REAL(param);
  double *Nph = REAL(nph);
  double *Leint = REAL(leint);
  double *X = REAL(x);
  int *Whint = INTEGER(whint);
  double *Knots = REAL(knots);
  double *Result = REAL(result);
  int nnph = lnph/lx;
  int lsdk = lk-1;
  int i, ii, k, z, Which, t2;
  double Temp, Beta1, Beta2;
  double *MyParam = malloc((lsdk+1)*sizeof(double));
  
  MyParam[0] = 0;

  for (z=0; z<lx; z++){

    Which = Whint[z];
    Result[z] = 0;
    Result[z+lx] = 0;

    t2 = z*nnph;
    for (i=0; i<lsdk; i++){
      MyParam[i+1] = Param[i];
      ii = 1;
      while (ii<nnph){
	MyParam[i+1] += Param[ii*lsdk+i]*Nph[ii+t2];
	ii++;
      }
    }

    Beta1 = MyParam[Which];
    Beta2 = MyParam[Which+1];
    Temp = Beta2-Beta1;
    if (Temp!=0){
      Result[z] += (1/Leint[Which])*(Beta1*(Knots[Which+1]-X[z])+Beta2*(X[z]-Knots[Which]));
      Result[z+lx] += (Leint[Which]/Temp)*(exp(Result[z])-exp(Beta1));
      }
    else {
      Result[z] += Beta1;
      Result[z+lx] += (X[z]-Knots[Which])*exp(Beta1);
    }

    if (Which>0){
      for (k=Which; k>0; k--){
	Beta1 = MyParam[k-1];
	Beta2 = MyParam[k];
	Temp = Beta2-Beta1;
	if (Temp!=0) {
	  Result[z+lx] += (Leint[k-1]/Temp)*(exp(Beta2)-exp(Beta1));
	}
	else {
	  Result[z+lx] += Leint[k-1]*exp(Beta1);
	}
      }
    }

    Result[z+lx] = log(Result[z+lx]);
  }

  UNPROTECT(7);
  free(MyParam);
  return result;
}
