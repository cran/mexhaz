/*********************************************/
/* Routine for estimating the gradient and   */
/* the Hessian of the log-hazard and of the  */
/* cumulative hazard                         */
/* (Weibull hazard,                          */ 
/* 1 time, no expected hazard)               */
/* Author: H. Charvat                        */
/* Last modified: 2020/09/05                 */
/* Part of the mexhaz 2.0 package            */
/*********************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include "SplineFunc.h"

SEXP HGH_WbR(SEXP x, SEXP nph, SEXP fixobs, SEXP param, SEXP paramf)
{
  SEXP loghaz, hazcum, test, gradlhaz, gradcum, hesslhaz, hesscum, rlist, rlistnames;
  int lx = length(x);
  int lnph = length(nph);
  int lfix = length(fixobs);
  int npar = length(param)+length(paramf);
  int nhess = 0.5*npar*(npar+1);

  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(fixobs = coerceVector(fixobs,REALSXP));
  PROTECT(param = coerceVector(param,REALSXP));
  PROTECT(paramf = coerceVector(paramf,REALSXP));
  PROTECT(loghaz = allocVector(REALSXP,lx));
  PROTECT(hazcum = allocVector(REALSXP,lx));
  PROTECT(gradlhaz = allocVector(REALSXP,lx*npar));
  PROTECT(gradcum = allocVector(REALSXP,lx*npar));
  PROTECT(hesslhaz = allocVector(REALSXP,lx*nhess));
  PROTECT(hesscum = allocVector(REALSXP,lx*nhess));
  PROTECT(test = allocVector(LGLSXP,1));
  int nprotect = 12;

  double *X = REAL(x);
  double *Nph = REAL(nph);
  double *FixObs = REAL(fixobs);
  double *Param = REAL(param);
  double *ParamF = REAL(paramf);
  double *LogHaz = REAL(loghaz);
  double *HazCum = REAL(hazcum);

  int nnph = lnph/lx;
  int nfix = lfix/lx;
  int i, j, z, t2, t3, cc;
  double tempF, tempH, tempLC, tempC, expTH, expTF;
  double tempCLC, tempLCLC;
  double Total = 0;

  double **GradLHaz = dmatrix(REAL(gradlhaz), lx, npar);
  double **GradCum = dmatrix(REAL(gradcum), lx, npar);
  double **HessLHaz = dmatrix(REAL(hesslhaz), lx, nhess);
  double **HessCum = dmatrix(REAL(hesscum), lx, nhess);

  for (z=0; z<lx; z++){

    tempF = 0;
    t3 = nfix*z;
    for (i=0; i<nfix; i++){
      tempF += FixObs[i+t3]*ParamF[i];
    }
    expTF = exp(tempF);

    t2 = z*nnph;
    tempH = 0;
    for (i=0; i<nnph; i++){
      tempH += Param[i]*Nph[i+t2];
    }
    expTH = exp(tempH);
    tempLC = log(X[z])*expTH;
    
    LogHaz[z] = log(X[z])*(expTH-1)+tempH;
    tempC = pow(X[z],expTH);
    Total += LogHaz[z] + tempC + tempF;
    LogHaz[z] += tempF;
    tempC *=expTF;
    HazCum[z] = tempC;
    tempCLC = tempLC*tempC;
    tempLCLC = (1+tempLC)*tempCLC;
    
    for (i=0; i<nfix; i++){
      GradLHaz[i][z] = FixObs[i+t3];
      GradCum[i][z] = FixObs[i+t3]*tempC;
    }

    for (i=0; i<nnph; i++){
      GradLHaz[nfix+i][z] = Nph[i+t2]*(1+tempLC);
      GradCum[nfix+i][z] = Nph[i+t2]*tempCLC;
    }

    cc = 0;
    for (i=0; i<nfix; i++){
      for (j=i; j<npar; j++){
	HessLHaz[cc][z] = 0;
	HessCum[cc][z] = FixObs[i+t3]*GradCum[j][z];
	cc++;
      }
    }
    
    for (i=0; i<nnph; i++){
      for (j=i; j<nnph; j++){
	HessLHaz[cc][z] = Nph[i+t2]*Nph[j+t2]*tempLC;
	HessCum[cc][z] = Nph[i+t2]*Nph[j+t2]*tempLCLC;
	cc++;
      }
    }
    
  }
  LOGICAL(test)[0] = (isinf(fabs(Total)) || isnan(Total));

  /* assemble the return objects as a list */
  PROTECT(rlist= allocVector(VECSXP, 7));
  SET_VECTOR_ELT(rlist, 0, loghaz);
  SET_VECTOR_ELT(rlist, 1, hazcum);
  SET_VECTOR_ELT(rlist, 2, test);
  SET_VECTOR_ELT(rlist, 3, gradlhaz);
  SET_VECTOR_ELT(rlist, 4, gradcum);
  SET_VECTOR_ELT(rlist, 5, hesslhaz);
  SET_VECTOR_ELT(rlist, 6, hesscum);
  
  /* add names to the list elements */
  PROTECT(rlistnames = allocVector(STRSXP, 7));
  SET_STRING_ELT(rlistnames, 0, mkChar("LogHaz"));
  SET_STRING_ELT(rlistnames, 1, mkChar("HazCum"));
  SET_STRING_ELT(rlistnames, 2, mkChar("Test"));
  SET_STRING_ELT(rlistnames, 3, mkChar("GradLogHaz"));
  SET_STRING_ELT(rlistnames, 4, mkChar("GradCum"));
  SET_STRING_ELT(rlistnames, 5, mkChar("HessLHaz"));
  SET_STRING_ELT(rlistnames, 6, mkChar("HessCum"));
  setAttrib(rlist, R_NamesSymbol, rlistnames);

  UNPROTECT(nprotect+2);
  return rlist;
}
