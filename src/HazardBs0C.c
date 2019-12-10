/************************************************/
/* Routine for estimating the cumulative hazard */
/* (piecewise constant hazard)                  */
/* For counting process data (t0, t1)           */
/* Author: H. Charvat                           */
/* Last modified: 2019/11/26                    */
/* Part of the mexhaz 1.7 package               */
/************************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>

SEXP HazardBs0C(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP matk)
{
  SEXP loghaz, hazcum0, hazcum, test, rlist, rlistnames;
  int lx = length(x);
  int lnph = length(nph);
  int lfix = length(fixobs);
  int nptd = length(param);

  PROTECT(x0 = coerceVector(x0,REALSXP));
  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(timecat0 = coerceVector(timecat0,INTSXP));
  PROTECT(timecat = coerceVector(timecat,INTSXP));
  PROTECT(fixobs = coerceVector(fixobs,REALSXP));
  PROTECT(param = coerceVector(param,REALSXP));
  PROTECT(paramf = coerceVector(paramf,REALSXP));
  PROTECT(matk = coerceVector(matk,REALSXP));
  PROTECT(loghaz = allocVector(REALSXP,lx));
  PROTECT(hazcum0 = allocVector(REALSXP,1));
  PROTECT(hazcum = allocVector(REALSXP,lx));
  PROTECT(test = allocVector(LGLSXP,1));
  int nprotect = 13;

  double *X0 = REAL(x0);
  double *X = REAL(x);
  double *Nph = REAL(nph);
  int *TimeCat0 = INTEGER(timecat0);
  int *TimeCat = INTEGER(timecat);
  double *FixObs = REAL(fixobs);
  double *Param = REAL(param);
  double *ParamF = REAL(paramf);
  double *MatK = REAL(matk);
  double *LogHaz = REAL(loghaz);
  double *HazCum = REAL(hazcum);

  int nnph = lnph/lx;
  int nfix = lfix/lx;
  int nbase = nptd/nnph;
  int i, k, z, t1, t2, tcz0, tcz;
  double Temp, tempF;
  double Total = 0;

  for (z=0; z<lx; z++){
    tempF = 0;
    t1 = z*nfix;
    for (i=0; i<nfix; i++){
      tempF += FixObs[i+t1]*ParamF[i];
    }
    tcz0 = TimeCat0[z];
    tcz = TimeCat[z];
    t2 = z*nnph;
    LogHaz[z] = 0;
    for (i=0; i<nnph; i++){
      LogHaz[z] += Param[nbase*i+tcz]*Nph[i+t2];
    }
    HazCum[z] = exp(LogHaz[z])*X[z];

    for (k=tcz; k>tcz0; k--){
      Temp = 0;
      for (i=0; i<nnph; i++){
	Temp += Param[nbase*i+(k-1)]*Nph[i+t2];
      }
      HazCum[z] += exp(Temp)*MatK[k-1];
    }

    Temp = 0;
    for (i=0; i<nnph; i++){
      Temp += Param[nbase*i+tcz0]*Nph[i+t2];
    }
    HazCum[z] -= exp(Temp)*X0[z];

    Total += LogHaz[z] + log(HazCum[z]) + tempF;
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
