/*********************************************/
/* Routine for estimating the variance of    */
/* the log-hazard and log-cumulative hazard  */
/* by the Delta Method                       */
/* (log-hazard described by                  */
/* a quadratic or cubic B-spline)            */
/* Author: H. Charvat                        */
/* Last modified: 2017/02/05                 */
/* Part of the mexhaz 1.4 package            */
/*********************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include "SplineFunc.h"

SEXP DeltaBs23R(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP paramt, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP varcov, SEXP grad)
{
  SEXP varlhaz, varlcum, gradlhaz, gradlcum, rlist, rlistnames;
  int lx = length(x);
  int lnph = length(nph);
  int lfix = length(fixobs);
  int ltotk = length(totk);
  int lleg = length(n);
  int npar = length(paramt);

  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(timecat = coerceVector(timecat,INTSXP));
  PROTECT(fixobs = coerceVector(fixobs,REALSXP));
  PROTECT(paramt = coerceVector(paramt,REALSXP));
  PROTECT(deg = coerceVector(deg,INTSXP));
  PROTECT(n = coerceVector(n,REALSXP));
  PROTECT(lw = coerceVector(lw,REALSXP));
  PROTECT(matk = coerceVector(matk,REALSXP));
  PROTECT(totk = coerceVector(totk,REALSXP));
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
  int nprotect = 16;

  double *X = REAL(x);
  double *Nph = REAL(nph);
  int *TimeCat = INTEGER(timecat);
  double *FixObs = REAL(fixobs);
  double *ParamT = REAL(paramt);
  int Deg = INTEGER(deg)[0];
  double *N = REAL(n);
  double *lW = REAL(lw);
  double *MatK = REAL(matk);
  double *TotK = REAL(totk);
  double *Varcov = REAL(varcov);
  double *VarLHaz = REAL(varlhaz);
  double *VarLCum = REAL(varlcum);

  int nnph = lnph/lx;
  int nfix = lfix/lx;
  int nbase = ltotk-Deg;
  int Cst1 = 2*Deg;
  int Cst2 = 2*Deg-2;
  int i, j, z, t2, t3, tcz;
  double tempL, InvtempL;
  double *TotKPos = &TotK[Deg];

  double *MyGradLH = (double *)R_alloc(npar,sizeof(double));
  double *MyGradLC = (double *)R_alloc(npar,sizeof(double));
  double *MyParam = (double *)R_alloc((nbase+1),sizeof(double));
  double *TempD = (double *)R_alloc(Cst1,sizeof(double));
  double *Res = (double *)R_alloc((nbase+1),sizeof(double));
  double *tempLvec = (double *)R_alloc((nbase+1),sizeof(double));
  
  double **GradLHaz = dmatrix(REAL(gradlhaz), A1, A2);
  double **GradLCum = dmatrix(REAL(gradlcum), A1, A2);

  double (*Fpt)(double, double*, double*, double*, double*, int, int, double*);
  if (Deg==2){
    Fpt = &DeltaSpline2;
  }
  else {
    Fpt = &DeltaSpline3;
  }

  MyParam[0] = 0;

  for (z=0; z<lx; z++){

    t3 = nfix*z;
    for (i=0; i<nfix; i++){
      MyGradLH[i] = FixObs[i+t3];
      MyGradLC[i] = FixObs[i+t3];
    }

    tcz = TimeCat[z];
    VarLHaz[z] = 0;
    VarLCum[z] = 0;
    tempLvec[0] = 0;
    tempL = 0;

    t2 = z*nnph;
    for (i=0; i<nbase; i++){
      MyParam[i+1] = ParamT[i+nfix];
      tempLvec[i+1] = 0;
      for (j=1; j<nnph; j++){
	MyParam[i+1] += ParamT[j*nbase+i+nfix]*Nph[j+t2];
      }
    }

    for (i=0; i<tcz; i++){
      tempL += IntDSpline23((*Fpt), TotKPos[i-1], TotKPos[i], &TotK[i], &MatK[Cst2*i], TempD, &MyParam[i], N, lW, lleg, nbase, i, tempLvec, Res);
    }
    tempL += IntDSpline23((*Fpt), TotKPos[tcz-1], X[z], &TotK[tcz], &MatK[Cst2*tcz], TempD, &MyParam[tcz], N, lW, lleg, nbase, tcz, tempLvec, Res);
    InvtempL = 1/tempL;
    (*Fpt)(X[z], &TotK[tcz], &MatK[Cst2*tcz], TempD, &MyParam[tcz], nbase, tcz, Res);

    for (i=0; i<nnph; i++){
      for (j=0; j<nbase; j++){
	MyGradLH[nfix + i*nbase+j] = Res[j+1]*Nph[i+t2];
	MyGradLC[nfix + i*nbase+j] = tempLvec[j+1]*Nph[i+t2]*InvtempL;
      }
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
