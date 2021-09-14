/*********************************************/
/* Routine for estimating the gradient and   */
/* the Hessian of the log-hazard and of the  */
/* cumulative hazard (aggregated)            */
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

SEXP HGHAggr_WbR(SEXP x, SEXP nph, SEXP fixobs, SEXP statobs, SEXP nbyclust, SEXP param, SEXP paramf)
{
  SEXP loghaz, hazcum, test, gradlhaz, gradcum, hesslhaz, hesscum, rlist, rlistnames;
  int lx = length(x);
  int lnph = length(nph);
  int lfix = length(fixobs);
  int nclust = length(nbyclust);
  int npar = length(param)+length(paramf);
  int nhess = 0.5*npar*(npar+1);

  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(fixobs = coerceVector(fixobs,REALSXP));
  PROTECT(statobs = coerceVector(statobs,INTSXP));
  PROTECT(nbyclust = coerceVector(nbyclust,INTSXP));
  PROTECT(param = coerceVector(param,REALSXP));
  PROTECT(paramf = coerceVector(paramf,REALSXP));
  PROTECT(loghaz = allocVector(REALSXP,nclust));
  PROTECT(hazcum = allocVector(REALSXP,nclust));
  PROTECT(gradlhaz = allocVector(REALSXP,nclust*npar));
  PROTECT(gradcum = allocVector(REALSXP,nclust*npar));
  PROTECT(hesslhaz = allocVector(REALSXP,nclust*nhess));
  PROTECT(hesscum = allocVector(REALSXP,nclust*nhess));
  PROTECT(test = allocVector(LGLSXP,1));
  int nprotect = 14;

  double *X = REAL(x);
  double *Nph = REAL(nph);
  double *FixObs = REAL(fixobs);
  int *StatObs = INTEGER(statobs);
  int *NByClust = INTEGER(nbyclust);
  double *Param = REAL(param);
  double *ParamF = REAL(paramf);
  double *LogHaz = REAL(loghaz);
  double *HazCum = REAL(hazcum);

  int nnph = lnph/lx;
  int nfix = lfix/lx;
  int Cst3 = nfix*npar-0.5*nfix*(nfix-1);
  int i, j, nc, m, t2, t3, cc;
  double tempF, tempLH, tempH, tempLC, tempC, expTH, expTF;
  double tempCLC, tempLCLC;
  double Total = 0;
  int z = 0;

  double *tempGC = (double *)R_alloc(npar,sizeof(double));

  double **GradLHaz = dmatrix(REAL(gradlhaz), nclust, npar);
  double **GradCum = dmatrix(REAL(gradcum), nclust, npar);
  double **HessLHaz = dmatrix(REAL(hesslhaz), nclust, nhess);
  double **HessCum = dmatrix(REAL(hesscum), nclust, nhess);
  
  for (nc=0; nc<nclust; nc++){
    
    LogHaz[nc] = 0;
    HazCum[nc] = 0;
    cc = 0;
    for (i=0; i<npar; i++){
      GradLHaz[i][nc] = 0;
      GradCum[i][nc] = 0;
      for (j=i; j<npar; j++){
	HessLHaz[cc][nc] = 0;
	HessCum[cc][nc] = 0;
	cc++;
      }
    }
    
    for (m=0; m<NByClust[nc]; m++){

      t3 = nfix*z;
      tempF = 0;
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
      
      tempLH = log(X[z])*(expTH-1)+tempH;
      tempC = pow(X[z],expTH);
      
      Total = tempLH + tempC + tempF;
      tempLH += tempF;
      LogHaz[nc] += StatObs[z]*tempLH;
      tempC *=expTF;
      HazCum[nc] += tempC;
      tempCLC = tempLC*tempC;
      tempLCLC = (1+tempLC)*tempCLC;
      
      for (i=0; i<nfix; i++){
	GradLHaz[i][nc] += StatObs[z]*FixObs[i+t3];
	tempGC[i] = FixObs[i+t3]*tempC;
	GradCum[i][nc] += tempGC[i];
      }
      
      cc = Cst3;
      for (i=0; i<nnph; i++){
	GradLHaz[nfix+i][nc] += Nph[i+t2]*(1+tempLC)*StatObs[z];
	tempGC[nfix+i] = Nph[i+t2]*tempCLC;
	GradCum[nfix+i][nc] += tempGC[nfix+i];
	for (j=i; j<nnph; j++){
	  HessLHaz[cc][nc] += Nph[i+t2]*Nph[j+t2]*tempLC*StatObs[z];
	  HessCum[cc][nc] += Nph[i+t2]*Nph[j+t2]*tempLCLC;
	  cc++;
	}
      }
      
      cc = 0;
      for (i=0; i<nfix; i++){
	for (j=i; j<npar; j++){
	  HessLHaz[cc][nc] += 0;
	  HessCum[cc][nc] += FixObs[i+t3]*tempGC[j];
	  cc++;
	}
      }
      
      z++;

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
