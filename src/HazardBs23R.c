/************************************************/
/* Routine for estimating the cumulative hazard */
/* (when the log-hazard is described            */
/* by a quadratic or cubic B-spline)            */
/* Author: H. Charvat                           */
/* Last modified: 2019/11/26                    */
/* Part of the mexhaz 1.7 package               */
/************************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include "SplineFunc.h"

SEXP HazardBs23R(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk)
{
  SEXP loghaz, hazcum0, hazcum, test, rlist, rlistnames;
  int lx = length(x);
  int lnph = length(nph);
  int lleg = length(n);
  int ltotk = length(totk);
  int lfix = length(fixobs);

  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(timecat = coerceVector(timecat,INTSXP));
  PROTECT(fixobs = coerceVector(fixobs,REALSXP));
  PROTECT(param = coerceVector(param,REALSXP));
  PROTECT(paramf = coerceVector(paramf,REALSXP));
  PROTECT(deg = coerceVector(deg,INTSXP));
  PROTECT(n = coerceVector(n,REALSXP));
  PROTECT(lw = coerceVector(lw,REALSXP));
  PROTECT(matk = coerceVector(matk,REALSXP));
  PROTECT(totk = coerceVector(totk,REALSXP));
  PROTECT(loghaz = allocVector(REALSXP,lx));
  PROTECT(hazcum0 = allocVector(REALSXP,1));
  PROTECT(hazcum = allocVector(REALSXP,lx));
  PROTECT(test = allocVector(LGLSXP,1));
  int nprotect = 15;

  double *X = REAL(x);
  double *Nph = REAL(nph);
  int *TimeCat = INTEGER(timecat);
  double *FixObs = REAL(fixobs);
  double *Param = REAL(param);
  double *ParamF = REAL(paramf);
  int Deg = INTEGER(deg)[0];
  double *N = REAL(n);
  double *lW = REAL(lw);
  double *MatK = REAL(matk);
  double *TotK = REAL(totk);
  double *LogHaz = REAL(loghaz);
  double *HazCum = REAL(hazcum);

  int nnph = lnph/lx;
  int nfix = lfix/lx;
  int nbase = ltotk-Deg;
  int Cst1 = 2*Deg;
  int Cst2 = 2*Deg-2;
  int i, j, z, tcz, t1;
  double tempL, tempH, tempF;
  double Total = 0;
  double *TotKPos = &TotK[Deg];

  double *MyParam = (double *)R_alloc((nbase+1),sizeof(double));
  double *TempD = (double *)R_alloc(Cst1,sizeof(double));

  double (*Fpt)(double, double*, double*, double*, double*);
  if (Deg==2){
    Fpt = &Spline2;
  }
  else {
    Fpt = &Spline3;
  }

  MyParam[0] = 0;

  if (nnph==1){
    int lk = nbase-Deg+1;
    double *CumVecSpl = (double *)R_alloc((lk+1),sizeof(double));
    double tempV;
    for (i=0; i<nbase; i++){
      MyParam[i+1] = Param[i];
    }
    CumVecSpl[0] = 0;
    for (i=0; i<lk; i++){
      tempV = IntSpline23((*Fpt), TotKPos[i-1], TotKPos[i], &TotK[i], &MatK[Cst2*i], TempD, &MyParam[i], N, lW, lleg);
      CumVecSpl[i+1] = CumVecSpl[i] + tempV;
    }
    for (z=0; z<lx; z++){
      tempF = 0;
      t1 = z*nfix;
      for (i=0; i<nfix; i++){
	tempF += FixObs[i+t1]*ParamF[i];
      }
      tcz = TimeCat[z];
      tempL = CumVecSpl[tcz];
      tempL += IntSpline23((*Fpt), TotKPos[tcz-1], X[z], &TotK[tcz], &MatK[Cst2*tcz], TempD, &MyParam[tcz], N, lW, lleg);
      tempH = (*Fpt)(X[z], &TotK[tcz], &MatK[Cst2*tcz], TempD, &MyParam[tcz]);
      Total += tempL + tempH + tempF;
      LogHaz[z] = tempH + tempF;
      HazCum[z] = tempL*exp(tempF);
    }
  }
  else {
    int t2;
    for (z=0; z<lx; z++){
      tempF = 0;
      t1 = z*nfix;
      for (i=0; i<nfix; i++){
	tempF += FixObs[i+t1]*ParamF[i];
      }
      t2 = z*nnph;
      for (i=0; i<nbase; i++){
	MyParam[i+1] = Param[i];
	for (j=1; j<nnph; j++){
	  MyParam[i+1] += Param[j*nbase+i]*Nph[j+t2];
	}
      }
      tcz = TimeCat[z];
      tempL = 0;
      for (i=0; i<tcz; i++){
        tempL += IntSpline23((*Fpt), TotKPos[i-1], TotKPos[i], &TotK[i], &MatK[Cst2*i], TempD, &MyParam[i], N, lW, lleg);
      }
      tempL += IntSpline23((*Fpt), TotKPos[tcz-1], X[z], &TotK[tcz], &MatK[Cst2*tcz], TempD, &MyParam[tcz], N, lW, lleg);
      tempH = (*Fpt)(X[z], &TotK[tcz], &MatK[Cst2*tcz], TempD, &MyParam[tcz]);
      Total += tempL + tempH + tempF;
      LogHaz[z] = tempH + tempF;
      HazCum[z] = tempL*exp(tempF);
    }
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
