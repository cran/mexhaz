/*********************************************/
/* Routine for estimating the variance of    */
/* the log-hazard and log-cumulative hazard  */
/* by the Delta Method                       */
/* (Weibull hazard)                          */
/* Author: H. Charvat                        */
/* Last modified: 2017/12/08                 */
/* Part of the mexhaz 1.4 package            */
/*********************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include "SplineFunc.h"

SEXP DeltaWeibR(SEXP x, SEXP nph, SEXP fixobs, SEXP paramt, SEXP varcov, SEXP grad)
{
  SEXP varlhaz, varlcum, gradlhaz, gradlcum, rlist, rlistnames;
  int lx = length(x);
  int lnph = length(nph);
  int lfix = length(fixobs);
  int npar = length(paramt);

  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(fixobs = coerceVector(fixobs,REALSXP));
  PROTECT(paramt = coerceVector(paramt,REALSXP));
  PROTECT(varcov = coerceVector(varcov,REALSXP));
  PROTECT(grad = coerceVector(grad,INTSXP));
  PROTECT(varlhaz = allocVector(REALSXP,lx));
  PROTECT(varlcum = allocVector(REALSXP,lx));
  int isGrad = INTEGER(grad)[0];
  int A1 = 1;
  int A2 = 1;
  if (isGrad){
    A1 = lx;
    A2 = npar;
  }
  PROTECT(gradlhaz = allocVector(REALSXP,A1*A2));
  PROTECT(gradlcum = allocVector(REALSXP,A1*A2));
  int nprotect = 10;

  double *X = REAL(x);
  double *Nph = REAL(nph);
  double *FixObs = REAL(fixobs);
  double *ParamT = REAL(paramt);
  double *Varcov = REAL(varcov);
  double *VarLHaz = REAL(varlhaz);
  double *VarLCum = REAL(varlcum);

  int nnph = lnph/lx;
  int nfix = lfix/lx;
  int i, j, z, t2, t3;
  double TempH;

  double *MyGradLH = (double *)R_alloc(npar,sizeof(double));
  double *MyGradLC = (double *)R_alloc(npar,sizeof(double));

  double **GradLHaz = dmatrix(REAL(gradlhaz), A1, A2);
  double **GradLCum = dmatrix(REAL(gradlcum), A1, A2);

  for (z=0; z<lx; z++){

    t3 = nfix*z;
    MyGradLH[0] = 1;
    MyGradLC[0] = 1;
    for (i=0; i<nfix; i++){
      MyGradLH[i+1] = FixObs[i+t3];
      MyGradLC[i+1] = FixObs[i+t3];
    }

    VarLHaz[z] = 0;
    VarLCum[z] = 0;

    t2 = z*nnph;
    TempH = ParamT[nfix+1];
    for (i=0; i<nnph; i++){
      TempH += ParamT[nfix+2+i]*Nph[i+t2];
    }
    TempH = X[z]*exp(TempH); 

    MyGradLH[nfix+1] = 1 + TempH;
    MyGradLC[nfix+1] = TempH;
    for (i=0; i<nnph; i++){
      MyGradLH[nfix+2+i] = Nph[i+t2]*(1+TempH);
      MyGradLC[nfix+2+i] = Nph[i+t2]*TempH;
    }

    for (i=0; i<npar; i++){
      for (j=0; j<npar; j++){
	VarLHaz[z] += MyGradLH[i]*Varcov[j+npar*i]*MyGradLH[j];
	VarLCum[z] += MyGradLC[i]*Varcov[j+npar*i]*MyGradLC[j];
      }
      if (isGrad){
	GradLHaz[i][z] = MyGradLH[i];
	GradLCum[i][z] = MyGradLC[i];
      }
    }

  }

  if (isGrad){
    /* assemble the return objects as a list */
    PROTECT(rlist= allocVector(VECSXP, 4));
    SET_VECTOR_ELT(rlist, 0, varlhaz);
    SET_VECTOR_ELT(rlist, 1, varlcum);
    SET_VECTOR_ELT(rlist, 2, gradlhaz);
    SET_VECTOR_ELT(rlist, 3, gradlcum);

    /* add names to the list elements */
    PROTECT(rlistnames = allocVector(STRSXP, 4));
    SET_STRING_ELT(rlistnames, 0, mkChar("VarLogHaz"));
    SET_STRING_ELT(rlistnames, 1, mkChar("VarLogCum"));
    SET_STRING_ELT(rlistnames, 2, mkChar("GradLogHaz"));
    SET_STRING_ELT(rlistnames, 3, mkChar("GradLogCum"));
  }
  else {
    /* assemble the return objects as a list */
    PROTECT(rlist= allocVector(VECSXP, 2));
    SET_VECTOR_ELT(rlist, 0, varlhaz);
    SET_VECTOR_ELT(rlist, 1, varlcum);

    /* add names to the list elements */
    PROTECT(rlistnames = allocVector(STRSXP, 2));
    SET_STRING_ELT(rlistnames, 0, mkChar("VarLogHaz"));
    SET_STRING_ELT(rlistnames, 1, mkChar("VarLogCum"));
  }
  setAttrib(rlist, R_NamesSymbol, rlistnames);

  UNPROTECT(nprotect+2);
  return rlist;
}
