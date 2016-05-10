/*********************************************/
/* Routine for estimating the variance of    */
/* the log-hazard and log-cumulative hazard  */
/* by the Delta Method                       */
/* (piecewise constant hazard)               */
/* Author: H. Charvat                        */
/* Last modified: 2016/05/10                 */
/* Part of the mexhaz package                */
/*********************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>

SEXP DeltaPwCst(SEXP nph, SEXP param, SEXP others, SEXP leint, SEXP lerem, SEXP whint, SEXP varcov)
{
  SEXP result;
  int lx = length(lerem);
  int lnph = length(nph);
  int loth = length(others);
  int npar = length(param);
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(param = coerceVector(param,REALSXP));
  PROTECT(others = coerceVector(others,REALSXP));
  PROTECT(leint = coerceVector(leint,REALSXP));
  PROTECT(lerem = coerceVector(lerem,REALSXP));
  PROTECT(whint = coerceVector(whint,INTSXP));
  PROTECT(varcov = coerceVector(varcov,REALSXP));
  PROTECT(result = allocVector(REALSXP,(2*lx)));
  double *Nph = REAL(nph);
  double *Param = REAL(param);
  double *Others = REAL(others);
  double *Leint = REAL(leint);
  double *Lerem = REAL(lerem);
  int *Whint = INTEGER(whint);
  double *Varcov = REAL(varcov);
  double *Result = REAL(result);
  int nnph = lnph/lx;
  int nind = loth/lx;
  int nbase = (npar-nind)/nnph;
  double tempL, InvtempL;
  double *MyGradL = malloc(npar*sizeof(double));
  double *MyGradS = malloc(npar*sizeof(double));
  double *Res = malloc(nbase*sizeof(double));
  double *MyParam = malloc(nbase*sizeof(double));
  double *tempLvec = malloc(nbase*sizeof(double));
  int i, ii, j, k, z, t2, t3, Which;

  for (z=0; z<lx; z++){

    t3 = nind*z;
    for (i=0; i<nind; i++){
      MyGradL[i] = Others[i+t3];
      MyGradS[i] = Others[i+t3];
    }

    Which = Whint[z];
    Result[z] = 0;
    Result[lx+z] = 0;
    tempL = 0;

    t2 = z*nnph;
    for (i=0; i<nbase; i++){
      MyParam[i] = Param[i+nind];
      Res[i] = 0;
      tempLvec[i] = 0;
      ii = 1;
      while (ii<nnph){
	MyParam[i] += Param[ii*nbase+i+nind]*Nph[ii+t2];
	ii++;
      }
    }

    // Calculation of lambda, Lambda and necessary integrals //
    Res[Which] = 1;
    tempL = exp(MyParam[Which])*Lerem[z];
    tempLvec[Which] = tempL;
    if (Which>0){
      for (k=(Which-1); k>=0; k--){
	tempL += exp(MyParam[k])*Leint[k];
	tempLvec[k] = exp(MyParam[k])*Leint[k];
      }
    }
    InvtempL = 1/tempL;

    ii = 0;
    while (ii<nnph){
      for (j=0; j<nbase; j++){
	MyGradL[nind + ii*nbase+j] = Res[j]*Nph[ii+t2];
	MyGradS[nind + ii*nbase+j] = tempLvec[j]*Nph[ii+t2]*InvtempL;
      }
      ii++;
    }

    for (i=0; i<npar; i++){
      for (j=0; j<npar; j++){
	Result[z] += MyGradL[i]*Varcov[j+npar*i]*MyGradL[j];
	Result[lx+z] += MyGradS[i]*Varcov[j+npar*i]*MyGradS[j];
      }
    }

  }

  UNPROTECT(8);
  free(Res);
  free(MyGradL);
  free(MyGradS);
  free(MyParam);
  free(tempLvec);
  return result;
}
