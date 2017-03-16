/************************************************/
/* Routine for estimating the cumulative hazard */
/* (when the log-hazard is described            */
/* by a restricted cubic B-spline)              */
/* For counting process data (t0, t1)           */
/* Author: H. Charvat                           */
/* Last modified: 2017/03/16                    */
/* Part of the mexhaz 1.3 package               */
/************************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include "SplineFunc.h"

SEXP HazardNsC(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP intk, SEXP nsadj1, SEXP nsadj2)
{
  SEXP loghaz, logcum, test, rlist, rlistnames;
  int lx = length(x);
  int lnph = length(nph);
  int lleg = length(n);
  int lfix = length(fixobs);

  PROTECT(x0 = coerceVector(x0,REALSXP));
  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(timecat0 = coerceVector(timecat0,INTSXP));
  PROTECT(timecat = coerceVector(timecat,INTSXP));
  PROTECT(fixobs = coerceVector(fixobs,REALSXP));
  PROTECT(param = coerceVector(param,REALSXP));
  PROTECT(paramf = coerceVector(paramf,REALSXP));
  PROTECT(deg = coerceVector(deg,INTSXP));
  PROTECT(n = coerceVector(n,REALSXP));
  PROTECT(lw = coerceVector(lw,REALSXP));
  PROTECT(matk = coerceVector(matk,REALSXP));
  PROTECT(totk = coerceVector(totk,REALSXP));
  PROTECT(intk = coerceVector(intk,REALSXP));
  PROTECT(nsadj1 = coerceVector(nsadj1,REALSXP));
  PROTECT(nsadj2 = coerceVector(nsadj2,REALSXP));
  PROTECT(loghaz = allocVector(REALSXP,lx));
  PROTECT(logcum = allocVector(REALSXP,lx));
  PROTECT(test = allocVector(LGLSXP,1));
  int nprotect = 19;

  double *X0 = REAL(x0);
  double *X = REAL(x);
  double *Nph = REAL(nph);
  int *TimeCat0 = INTEGER(timecat0);
  int *TimeCat = INTEGER(timecat);
  double *FixObs = REAL(fixobs);
  double *Param = REAL(param);
  double *ParamF = REAL(paramf);
  int *Deg = INTEGER(deg);
  double *N = REAL(n);
  double *lW = REAL(lw);
  double *MatK = REAL(matk);
  double *TotK = REAL(totk);
  double *IntK = REAL(intk);
  double *NsAdj1 = REAL(nsadj1);
  double *NsAdj2 = REAL(nsadj2);
  double *LogHaz = REAL(loghaz);
  double *LogCum = REAL(logcum);

  int nnph = lnph/lx;
  int nfix = lfix/lx;
  int nbase = Deg[1]-5;
  int leB = Deg[1]-1;
  int firstK = Deg[2];
  int i, j, z, tcz0, tcz, t1;
  double tempL, tempH, tempF;
  double Total = 0;
  double TempD[6];

  double *MyBasisB = (double *)R_alloc(leB,sizeof(double));

  if (nnph==1){
    for (z=0; z<lx; z++){
      tempF = 0;
      t1 = z*nfix;
      for (i=0; i<nfix; i++){
	tempF += FixObs[i+t1]*ParamF[i];
      }
      tcz0 = TimeCat0[z];
      tcz = TimeCat[z];
      tempL = 0;
      for (i=tcz0; i<tcz; i++){
	tempL += IntNSpl(IntK[i], IntK[i+1], &TotK[i], &MatK[4*i], NsAdj1, NsAdj2, MyBasisB, TempD, Param, N, lW, lleg, leB, nbase, (i+firstK));
      }
      tempL += IntNSpl(IntK[tcz], X[z], &TotK[tcz], &MatK[4*tcz], NsAdj1, NsAdj2, MyBasisB, TempD, Param, N, lW, lleg, leB, nbase, (tcz+firstK));
      tempL -= IntNSpl(IntK[tcz0], X0[z], &TotK[tcz0], &MatK[4*tcz0], NsAdj1, NsAdj2, MyBasisB, TempD, Param, N, lW, lleg, leB, nbase, (tcz0+firstK));
      tempH = NSpl(X[z], &TotK[tcz], &MatK[4*tcz], NsAdj1, NsAdj2, MyBasisB, TempD, Param, leB, nbase, (tcz+firstK));
      tempL = log(tempL);
      Total += tempH + tempL + tempF;
      LogHaz[z] = tempH + tempF;
      LogCum[z] = tempL + tempF;
    }
  }
  else {
    int t2;
    double *MyParam = (double *)R_alloc(nbase,sizeof(double));
    for (z=0; z<lx; z++){
      tempF = 0;
      t1 = z*nfix;
      for (i=0; i<nfix; i++){
	tempF += FixObs[i+t1]*ParamF[i];
      }
      t2 = z*nnph;
      for (i=0; i<nbase; i++){
	MyParam[i] = Param[i];
	for (j=1; j<nnph; j++){
	  MyParam[i] += Param[j*nbase+i]*Nph[j+t2];
	}
      }
      tcz0 = TimeCat0[z];
      tcz = TimeCat[z];
      tempL = 0;
      for (i=tcz0; i<tcz; i++){
	tempL += IntNSpl(IntK[i], IntK[i+1], &TotK[i], &MatK[4*i], NsAdj1, NsAdj2, MyBasisB, TempD, Param, N, lW, lleg, leB, nbase, (i+firstK));
      }
      tempL += IntNSpl(IntK[tcz], X[z], &TotK[tcz], &MatK[4*tcz], NsAdj1, NsAdj2, MyBasisB, TempD, Param, N, lW, lleg, leB, nbase, (tcz+firstK));
      tempL -= IntNSpl(IntK[tcz0], X0[z], &TotK[tcz0], &MatK[4*tcz0], NsAdj1, NsAdj2, MyBasisB, TempD, Param, N, lW, lleg, leB, nbase, (tcz0+firstK));
      tempH = NSpl(X[z], &TotK[tcz], &MatK[4*tcz], NsAdj1, NsAdj2, MyBasisB, TempD, Param, leB, nbase, (tcz+firstK));
      tempL = log(tempL);
      Total += tempL + tempH + tempF;
      LogHaz[z] = tempH + tempF;
      LogCum[z] = tempL + tempF;
    }
  }
  LOGICAL(test)[0] = (isinf(fabs(Total)) || isnan(Total));

  /* assemble the return objects as a list */
  PROTECT(rlist= allocVector(VECSXP, 3));
  SET_VECTOR_ELT(rlist, 0, loghaz);
  SET_VECTOR_ELT(rlist, 1, logcum);
  SET_VECTOR_ELT(rlist, 2, test);

  /* add names to the list elements */
  PROTECT(rlistnames = allocVector(STRSXP, 3));
  SET_STRING_ELT(rlistnames, 0, mkChar("LogHaz"));
  SET_STRING_ELT(rlistnames, 1, mkChar("LogCum"));
  SET_STRING_ELT(rlistnames, 2, mkChar("Test"));
  setAttrib(rlist, R_NamesSymbol, rlistnames);

  UNPROTECT(nprotect+2);
  return rlist;
}