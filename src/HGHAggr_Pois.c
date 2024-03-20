/*********************************************/
/* Routine for estimating the gradient and   */
/* the Hessian of the log-hazard and of the  */
/* cumulative hazard (aggregated)            */
/* (Poisson rate, aggregated data)           */
/* Author: H. Charvat                        */
/* Last modified: 2024/03/20                 */
/* Part of the mexhaz 2.5 package            */
/*********************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include "SplineFunc.h"

SEXP HGHAggr_Pois(SEXP fixobs, SEXP statobs, SEXP offobs, SEXP nbyclust, SEXP paramf)
{
  SEXP loghaz, hazcum, test, gradlhaz, gradcum, hesslhaz, hesscum, rlist, rlistnames;
  int lx = length(statobs);
  int lfix = length(fixobs);
  int nclust = length(nbyclust);
  int npar = length(paramf);
  int nhess = 0.5*npar*(npar+1);
  
  PROTECT(fixobs = coerceVector(fixobs,REALSXP));
  PROTECT(statobs = coerceVector(statobs,INTSXP));
  PROTECT(offobs = coerceVector(offobs,REALSXP));
  PROTECT(nbyclust = coerceVector(nbyclust,INTSXP));
  PROTECT(paramf = coerceVector(paramf,REALSXP));
  PROTECT(loghaz = allocVector(REALSXP,nclust));
  PROTECT(hazcum = allocVector(REALSXP,nclust));
  PROTECT(gradlhaz = allocVector(REALSXP,nclust*npar));
  PROTECT(gradcum = allocVector(REALSXP,nclust*npar));
  PROTECT(hesslhaz = allocVector(REALSXP,1));
  PROTECT(hesscum = allocVector(REALSXP,nclust*nhess));
  PROTECT(test = allocVector(LGLSXP,1));
  int nprotect = 12;

  double *FixObs = REAL(fixobs);
  int *StatObs = INTEGER(statobs);
  double *OffObs = REAL(offobs);
  int *NByClust = INTEGER(nbyclust);
  double *ParamF = REAL(paramf);
  double *LogHaz = REAL(loghaz);
  double *HazCum = REAL(hazcum);

  int nfix = lfix/lx;
  int i, j, nc, m, t3, cc;
  double tempF, expTF;
  double Total = 0;
  int z = 0;

  double *tempGC = (double *)R_alloc(npar,sizeof(double));
  double **GradLHaz = dmatrix(REAL(gradlhaz), nclust, npar);
  double **GradCum = dmatrix(REAL(gradcum), nclust, npar);
  double **HessCum = dmatrix(REAL(hesscum), nclust, nhess);

  for (nc=0; nc<nclust; nc++){

    LogHaz[nc] = 0;
    HazCum[nc] = 0;
    cc = 0;
    for (i=0; i<npar; i++){
      GradLHaz[i][nc] = 0;
      GradCum[i][nc] = 0;
      for (j=i; j<npar; j++){
	HessCum[cc][nc] = 0;
	cc++;
      }
    }

    for (m=0; m<NByClust[nc]; m++){
    
      tempF = 0;
      t3 = nfix*z;
      for (i=0; i<nfix; i++){
	tempF += FixObs[i+t3]*ParamF[i];
      }
      //tempF += OffObs[z];
      expTF = exp(tempF + OffObs[z]);
            
      // Calculation of lambda, Lambda //
      
      Total += tempF + OffObs[z];
      LogHaz[nc] += StatObs[z]*tempF;
      HazCum[nc] += expTF;
      
      for (i=0; i<nfix; i++){
	GradLHaz[i][nc] += StatObs[z]*FixObs[i+t3];
	tempGC[i] = FixObs[i+t3]*expTF;
	GradCum[i][nc] += tempGC[i];
      }
      
      cc = 0;
      for (i=0; i<nfix; i++){
	for (j=i; j<npar; j++){
	  HessCum[cc][nc] += FixObs[i+t3]*tempGC[j];
	  cc++;
	}
      }

      z++;

    } 
  }
  LOGICAL(test)[0] = (isinf(fabs(Total)) || isnan(Total));
  REAL(hesslhaz)[0] = 0;

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
