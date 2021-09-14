/*********************************************/
/* Routine for estimating the gradient and   */
/* the Hessian of the log-hazard and of the  */
/* cumulative hazard (aggregated)            */
/* (log-hazard described by a B-spline,      */
/* 2 times, no expected hazard)              */
/* Author: H. Charvat                        */
/* Last modified: 2020/09/05                 */
/* Part of the mexhaz 2.0 package            */
/*********************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include "SplineFunc.h"

SEXP HGHAggr_BsL(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk)
{
  SEXP loghaz, hazcum0, hazcum, test, gradlhaz, gradcum0, gradcum, hesslhaz, hesscum0, hesscum, rlist, rlistnames;
  int lx = length(x);
  int lnph = length(nph);
  int lfix = length(fixobs);
  int ltotk = length(totk);
  int lleg = length(n);
  int nclust = length(nbyclust);
  int npar = length(param)+length(paramf);
  int nhess = 0.5*npar*(npar+1);
  
  PROTECT(x0 = coerceVector(x0,REALSXP));
  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(timecat0 = coerceVector(timecat0,INTSXP));
  PROTECT(timecat = coerceVector(timecat,INTSXP));
  PROTECT(fixobs = coerceVector(fixobs,REALSXP));
  PROTECT(statobs = coerceVector(statobs,INTSXP));
  PROTECT(nbyclust = coerceVector(nbyclust,INTSXP));
  PROTECT(param = coerceVector(param,REALSXP));
  PROTECT(paramf = coerceVector(paramf,REALSXP));
  PROTECT(deg = coerceVector(deg,INTSXP));
  PROTECT(n = coerceVector(n,REALSXP));
  PROTECT(lw = coerceVector(lw,REALSXP));
  PROTECT(matk = coerceVector(matk,REALSXP));
  PROTECT(totk = coerceVector(totk,REALSXP));
  PROTECT(loghaz = allocVector(REALSXP,nclust));
  PROTECT(hazcum = allocVector(REALSXP,nclust));
  PROTECT(hazcum0 = allocVector(REALSXP,nclust));
  PROTECT(gradlhaz = allocVector(REALSXP,nclust*npar));
  PROTECT(gradcum0 = allocVector(REALSXP,nclust*npar));
  PROTECT(gradcum = allocVector(REALSXP,nclust*npar));
  PROTECT(hesslhaz = allocVector(REALSXP,1));
  PROTECT(hesscum0 = allocVector(REALSXP,nclust*nhess));
  PROTECT(hesscum = allocVector(REALSXP,nclust*nhess));
  PROTECT(test = allocVector(LGLSXP,1));
  int nprotect = 25;

  double *X0 = REAL(x0);
  double *X = REAL(x);
  double *Nph = REAL(nph);
  int *TimeCat0 = INTEGER(timecat0);
  int *TimeCat = INTEGER(timecat);
  double *FixObs = REAL(fixobs);
  int *StatObs = INTEGER(statobs);
  int *NByClust = INTEGER(nbyclust);
  double *Param = REAL(param);
  double *ParamF = REAL(paramf);
  int Deg = INTEGER(deg)[0];
  double *N = REAL(n);
  double *lW = REAL(lw);
  double *MatK = REAL(matk);
  double *TotK = REAL(totk);
  double *LogHaz = REAL(loghaz);
  double *HazCum0 = REAL(hazcum0);
  double *HazCum = REAL(hazcum);

  int nnph = lnph/lx;
  int nfix = lfix/lx;
  int nbase = ltotk-Deg;
  int Cst1 = 2*Deg;
  int Cst2;
  int Cst3 = nfix*npar-0.5*nfix*(nfix-1);
  int i, j, k, l, nc, m, t2, t3, tcz, tcz0, cc;
  double tempF, tempH, tempL, tempL0, expTF, NexpTF, NiNeTF, NlNeTF;
  double Total = 0;
  int z = 0;
  double *TotKPos = &TotK[Deg];

  double *MyParam = (double *)R_alloc(nbase+1,sizeof(double));
  double *TempD = (double *)R_alloc(Cst1,sizeof(double));
  double *Res = (double *)R_alloc(nbase+1,sizeof(double));
  double *tempLvec = (double *)R_alloc(nbase,sizeof(double));
  double *tempHess = (double *)R_alloc(nbase*nbase,sizeof(double));
  double *Res0 = (double *)R_alloc(nbase+1,sizeof(double));
  double *tempLvec0 = (double *)R_alloc(nbase,sizeof(double));
  double *tempHess0 = (double *)R_alloc(nbase*nbase,sizeof(double));
  double *tempGC = (double *)R_alloc(npar,sizeof(double));
  double *tempGC0 = (double *)R_alloc(npar,sizeof(double));
  
  double **GradLHaz = dmatrix(REAL(gradlhaz), nclust, npar);
  double **GradCum0 = dmatrix(REAL(gradcum0), nclust, npar);
  double **GradCum = dmatrix(REAL(gradcum), nclust, npar);
  double **HessCum0 = dmatrix(REAL(hesscum0), nclust, nhess);
  double **HessCum = dmatrix(REAL(hesscum), nclust, nhess);

  double (*Fpt)(double, double*, double*, double*, double*, int, int, double*);
  if (Deg==1){
    Fpt = &DeltaSpline1;
    Cst2 = 1;
  }
  else if (Deg==2){
    Fpt = &DeltaSpline2;
    Cst2 = 2;
  }
  else {
    Fpt = &DeltaSpline3;
    Cst2 = 4;
  }

  MyParam[0] = 0;

  for (nc=0; nc<nclust; nc++){
    
    LogHaz[nc] = 0;
    HazCum[nc] = 0;
    HazCum0[nc] = 0;
    cc = 0;
    for (i=0; i<npar; i++){
      GradLHaz[i][nc] = 0;
      GradCum[i][nc] = 0;
      GradCum0[i][nc] = 0;
      for (j=i; j<npar; j++){
	HessCum[cc][nc] = 0;
	HessCum0[cc][nc] = 0;
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
      
      tcz0 = TimeCat0[z];
      tcz = TimeCat[z];
      tempL0 = 0;
      tempL = 0;
      
      t2 = z*nnph;
      for (i=0; i<nbase; i++){
	MyParam[i+1] = Param[i];
	tempLvec0[i] = 0;
	tempLvec[i] = 0;
	for (j=1; j<nnph; j++){
	  MyParam[i+1] += Param[j*nbase+i]*Nph[j+t2];
	}
	for (k=0; k<nbase; k++){
	  tempHess0[k+i*nbase] = 0;
	  tempHess[k+i*nbase] = 0;
	}
      }
      
      for (i=0; i<tcz; i++){
	tempL += IntDSpline23H((*Fpt), TotKPos[i-1], TotKPos[i], &TotK[i], &MatK[Cst2*i], TempD, &MyParam[i], N, lW, lleg, nbase, i, tempLvec, tempHess, Res);
      }
      tempL += IntDSpline23H((*Fpt), TotKPos[tcz-1], X[z], &TotK[tcz], &MatK[Cst2*tcz], TempD, &MyParam[tcz], N, lW, lleg, nbase, tcz, tempLvec, tempHess, Res);
      tempH = (*Fpt)(X[z], &TotK[tcz], &MatK[Cst2*tcz], TempD, &MyParam[tcz], nbase, tcz, Res);
      
      for (i=0; i<tcz0; i++){
	tempL0 += IntDSpline23H((*Fpt), TotKPos[i-1], TotKPos[i], &TotK[i], &MatK[Cst2*i], TempD, &MyParam[i], N, lW, lleg, nbase, i, tempLvec0, tempHess0, Res0);
      }
      tempL0 += IntDSpline23H((*Fpt), TotKPos[tcz0-1], X0[z], &TotK[tcz0], &MatK[Cst2*tcz0], TempD, &MyParam[tcz0], N, lW, lleg, nbase, tcz0, tempLvec0, tempHess0, Res0);
      (*Fpt)(X0[z], &TotK[tcz0], &MatK[Cst2*tcz0], TempD, &MyParam[tcz0], nbase, tcz0, Res0);
      
      Total += tempL + tempH + tempF;
      LogHaz[nc] += StatObs[z]*(tempH + tempF);
      HazCum[nc] += tempL*expTF;
      HazCum0[nc] += tempL0*expTF;
      
      for (i=0; i<nfix; i++){
	GradLHaz[i][nc] += StatObs[z]*FixObs[i+t3];
	tempGC[i] = FixObs[i+t3]*tempL*expTF;
	tempGC0[i] = FixObs[i+t3]*tempL0*expTF;
	GradCum[i][nc] += tempGC[i];
	GradCum0[i][nc] += tempGC0[i];
      }

      cc = Cst3;
      for (i=0; i<nnph; i++){
	NexpTF = Nph[i+t2]*expTF;
	NiNeTF = Nph[i+t2]*NexpTF;
	for (j=0; j<nbase; j++){
	  GradLHaz[nfix + i*nbase+j][nc] += Res[j+1]*Nph[i+t2]*StatObs[z];
	  tempGC[nfix + i*nbase+j] = tempLvec[j]*NexpTF;
	  tempGC0[nfix + i*nbase+j] = tempLvec0[j]*NexpTF;
	  GradCum[nfix + i*nbase+j][nc] += tempGC[nfix + i*nbase+j];
	  GradCum0[nfix + i*nbase+j][nc] += tempGC0[nfix + i*nbase+j];
	  for (k=j; k<nbase; k++){
	    HessCum[cc][nc] += tempHess[k+j*nbase]*NiNeTF;
	    HessCum0[cc][nc] += tempHess0[k+j*nbase]*NiNeTF;
	    cc++;
	  }
	  for (l=i+1; l<nnph; l++){
	    NlNeTF = Nph[l+t2]*NexpTF;
	    for (k=0; k<nbase; k++){
	      HessCum[cc][nc] += tempHess[k+j*nbase]*NlNeTF;
	      HessCum0[cc][nc] += tempHess0[k+j*nbase]*NlNeTF;
	      cc++;
	    }
	  }	  
	}
      }
            
      cc = 0;
      for (i=0; i<nfix; i++){
	for (j=i; j<npar; j++){
	  HessCum[cc][nc] += FixObs[i+t3]*tempGC[j];
	  HessCum0[cc][nc] += FixObs[i+t3]*tempGC0[j];
	  cc++;
	}
      }
      
      z++;
      
    }
  }
  LOGICAL(test)[0] = (isinf(fabs(Total)) || isnan(Total));
  REAL(hesslhaz)[0] = 0;

  /* assemble the return objects as a list */
  PROTECT(rlist= allocVector(VECSXP, 10));
  SET_VECTOR_ELT(rlist, 0, loghaz);
  SET_VECTOR_ELT(rlist, 1, hazcum0);
  SET_VECTOR_ELT(rlist, 2, hazcum);
  SET_VECTOR_ELT(rlist, 3, test);
  SET_VECTOR_ELT(rlist, 4, gradlhaz);
  SET_VECTOR_ELT(rlist, 5, gradcum0);
  SET_VECTOR_ELT(rlist, 6, gradcum);
  SET_VECTOR_ELT(rlist, 7, hesslhaz);
  SET_VECTOR_ELT(rlist, 8, hesscum0);
  SET_VECTOR_ELT(rlist, 9, hesscum);
  
  /* add names to the list elements */
  PROTECT(rlistnames = allocVector(STRSXP, 10));
  SET_STRING_ELT(rlistnames, 0, mkChar("LogHaz"));
  SET_STRING_ELT(rlistnames, 1, mkChar("HazCum0"));
  SET_STRING_ELT(rlistnames, 2, mkChar("HazCum"));
  SET_STRING_ELT(rlistnames, 3, mkChar("Test"));
  SET_STRING_ELT(rlistnames, 4, mkChar("GradLogHaz"));
  SET_STRING_ELT(rlistnames, 5, mkChar("GradCum0"));
  SET_STRING_ELT(rlistnames, 6, mkChar("GradCum"));
  SET_STRING_ELT(rlistnames, 7, mkChar("HessLHaz"));
  SET_STRING_ELT(rlistnames, 8, mkChar("HessCum0"));
  SET_STRING_ELT(rlistnames, 9, mkChar("HessCum"));
  setAttrib(rlist, R_NamesSymbol, rlistnames);

  UNPROTECT(nprotect+2);
  return rlist;
}
