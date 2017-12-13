/*********************************************/
/* Routine for estimating the variance of    */
/* the log-hazard and log-cumulative hazard  */
/* by the Delta Method                       */
/* (log-hazard described by linear B-spline) */
/* Author: H. Charvat                        */
/* Last modified: 2017/02/05                 */
/* Part of the mexhaz 1.4 package            */
/*********************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include "SplineFunc.h"

SEXP DeltaBs1R(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP paramt, SEXP matk, SEXP totk, SEXP varcov, SEXP grad)
{
  SEXP varlhaz, varlcum, gradlhaz, gradlcum, rlist, rlistnames;
  int lx = length(x);
  int lnph = length(nph);
  int lfix = length(fixobs);
  int ltotk = length(totk);
  int npar = length(paramt);

  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(timecat = coerceVector(timecat,INTSXP));
  PROTECT(fixobs = coerceVector(fixobs,REALSXP));
  PROTECT(paramt = coerceVector(paramt,REALSXP));
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
  int nprotect = 13;

  double *X = REAL(x);
  double *Nph = REAL(nph);
  int *TimeCat = INTEGER(timecat);
  double *FixObs = REAL(fixobs);
  double *ParamT = REAL(paramt);
  double *MatK = REAL(matk);
  double *TotK = REAL(totk);
  double *Varcov = REAL(varcov);
  double *VarLHaz = REAL(varlhaz);
  double *VarLCum = REAL(varlcum);

  int nnph = lnph/lx;
  int nfix = lfix/lx;
  int nbase = ltotk-1;
  int i, j, k, z, t2, t3, tcz;
  double Temp, Beta1, Beta2;
  double Cst1 = 0;
  double Cst2 = 0;
  double templ, tempL, InvtempL;

  double *MyGradLH = (double *)R_alloc(npar,sizeof(double));
  double *MyGradLC = (double *)R_alloc(npar,sizeof(double));
  double *MyParam = (double *)R_alloc((nbase+1),sizeof(double));
  double *Res = (double *)R_alloc((nbase+1),sizeof(double));
  double *tempLvec = (double *)R_alloc((nbase+1),sizeof(double));

  double **GradLHaz = dmatrix(REAL(gradlhaz), A1, A2);
  double **GradLCum = dmatrix(REAL(gradlcum), A1, A2);

  Res[0] = 0;
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

    t2 = z*nnph;
    for (i=0; i<nbase; i++){
      MyParam[i+1] = ParamT[i+nfix];
      Res[i+1] = 0;
      tempLvec[i+1] = 0;
      for (j=1; j<nnph; j++){
	MyParam[i+1] += ParamT[j*nbase+i+nfix]*Nph[j+t2];
      }
    }

    // Calculation of lambda, Lambda and necessary integrals //
    Beta1 = MyParam[tcz];
    Beta2 = MyParam[tcz+1];
    Res[tcz] = (TotK[tcz+1]-X[z])/MatK[tcz];
    Res[tcz+1] = (X[z]-TotK[tcz])/MatK[tcz];
    Temp = Beta2-Beta1;
    if (Temp!=0){
      Cst1 = MatK[tcz]/Temp;
      Cst2 = 1/Temp;
      templ = Beta1*Res[tcz]+Beta2*Res[tcz+1];
      tempL = Cst1*(exp(templ)-exp(Beta1));
      tempLvec[tcz] = Cst1*((Res[tcz]+Cst2)*exp(templ)-(1+Cst2)*exp(Beta1));
      tempLvec[tcz+1] = Cst1*((Res[tcz+1]-Cst2)*exp(templ)+Cst2*exp(Beta1));
    }
    else {
      templ = Beta1;
      Cst1 = 0.5*exp(templ)*Res[tcz+1]*MatK[tcz];
      tempL = (X[z]-TotK[tcz])*exp(templ);
      tempLvec[tcz] = Cst1*(1+Res[tcz]);
      tempLvec[tcz+1] = Cst1*Res[tcz+1];
    }

    for (k=tcz; k>0; k--){
      Beta1 = MyParam[k-1];
      Beta2 = MyParam[k];
      Temp = Beta2-Beta1;
      if (Temp!=0) {
	Cst1 = MatK[k-1]/Temp;
	Cst2 = 1/Temp;
	tempL += Cst1*(exp(Beta2)-exp(Beta1));
	tempLvec[k] += Cst1*((1-Cst2)*exp(Beta2)+Cst2*exp(Beta1));
	tempLvec[k-1] += Cst1*(Cst2*exp(Beta2)-(1+Cst2)*exp(Beta1));
      }
      else {
	Cst1 = MatK[k-1]*exp(Beta1);
	tempL += Cst1;
	tempLvec[k] += 0.5*Cst1;
	tempLvec[k-1] += 0.5*Cst1;
      }
    }
    InvtempL = 1/tempL;

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
