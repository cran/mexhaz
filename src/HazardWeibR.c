/************************************************/
/* Routine for estimating the cumulative hazard */
/* (piecewise constant hazard)                  */
/* Author: H. Charvat                           */
/* Last modified: 2020/10/02                    */
/* Part of the mexhaz 1.9 package               */
/************************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>

SEXP HazardWeibR(SEXP x, SEXP nph, SEXP fixobs, SEXP param, SEXP paramf)
{
  SEXP loghaz, hazcum0, hazcum, test, rlist, rlistnames;
  int lx = length(x);
  int lnph = length(nph);
  int lfix = length(fixobs);

  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(fixobs = coerceVector(fixobs,REALSXP));
  PROTECT(param = coerceVector(param,REALSXP));
  PROTECT(paramf = coerceVector(paramf,REALSXP));
  PROTECT(loghaz = allocVector(REALSXP,lx));
  PROTECT(hazcum0 = allocVector(REALSXP,1));
  PROTECT(hazcum = allocVector(REALSXP,lx));
  PROTECT(test = allocVector(LGLSXP,1));
  int nprotect = 9;

  double *X = REAL(x);
  double *Nph = REAL(nph);
  double *FixObs = REAL(fixobs);
  double *Param = REAL(param);
  double *ParamF = REAL(paramf);
  double *LogHaz = REAL(loghaz);
  double *HazCum = REAL(hazcum);

  int nnph = lnph/lx;
  int nfix = lfix/lx;
  int i, z, t1, t2;
  double tempF, tempN;
  double Total = 0;

  for (z=0; z<lx; z++){
    tempF = ParamF[0];
    t1 = z*nfix;
    for (i=0; i<nfix; i++){
      tempF += FixObs[i+t1]*ParamF[i+1];
    }
    t2 = z*nnph;
    tempN = Param[0];
    for (i=0; i<nnph; i++){
      tempN += Param[i+1]*Nph[i+t2];
    }
    LogHaz[z] = log(X[z])*(exp(tempN)-1)+tempN;
    HazCum[z] = pow(X[z],exp(tempN));
    Total += LogHaz[z] + HazCum[z] + tempF;
    LogHaz[z] += tempF;
    HazCum[z] *= exp(tempF);
  }
  REAL(hazcum0)[0] = 0;
  LOGICAL(test)[0] = (isinf(fabs(Total)) || isnan(Total));

  /* assemble the return objects as a list */
  PROTECT(rlist= allocVector(VECSXP, 4));
  SET_VECTOR_ELT(rlist, 0, loghaz);
  SET_VECTOR_ELT(rlist, 1, hazcum0);
  SET_VECTOR_ELT(rlist, 2, hazcum);
  SET_VECTOR_ELT(rlist, 3, test);

  /* add names to the list elements */
  PROTECT(rlistnames = allocVector(STRSXP, 4));
  SET_STRING_ELT(rlistnames, 0, mkChar("LogHaz"));
  SET_STRING_ELT(rlistnames, 1, mkChar("HazCum0"));
  SET_STRING_ELT(rlistnames, 2, mkChar("HazCum"));
  SET_STRING_ELT(rlistnames, 3, mkChar("Test"));
  setAttrib(rlist, R_NamesSymbol, rlistnames);

  UNPROTECT(nprotect+2);
  return rlist;
}
