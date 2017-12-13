/*********************************************/
/* Routine for estimating the variance of    */
/* the log-hazard and log-cumulative hazard  */
/* by the Delta Method                       */
/* (log-hazard described by                  */
/* a restricted cubic B-spline)              */
/* Author: H. Charvat                        */
/* Last modified: 2017/03/16                 */
/* Part of the mexhaz 1.4 package            */
/*********************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include "SplineFunc.h"

SEXP DeltaNsR(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP paramt, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP intk, SEXP nsadj1, SEXP nsadj2, SEXP varcov, SEXP grad)
{
  SEXP varlhaz, varlcum, gradlhaz, gradlcum, rlist, rlistnames;
  int lx = length(x);
  int lnph = length(nph);
  int lleg = length(n);
  int lfix = length(fixobs);
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
  PROTECT(intk = coerceVector(intk,REALSXP));
  PROTECT(nsadj1 = coerceVector(nsadj1,REALSXP));
  PROTECT(nsadj2 = coerceVector(nsadj2,REALSXP));
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
  int nprotect = 19;

  double *X = REAL(x);
  double *Nph = REAL(nph);
  int *TimeCat = INTEGER(timecat);
  double *FixObs = REAL(fixobs);
  double *ParamT = REAL(paramt);
  int *Deg = INTEGER(deg);
  double *N = REAL(n);
  double *lW = REAL(lw);
  double *MatK = REAL(matk);
  double *TotK = REAL(totk);
  double *IntK = REAL(intk);
  double *NsAdj1 = REAL(nsadj1);
  double *NsAdj2 = REAL(nsadj2);
  double *Varcov = REAL(varcov);
  double *VarLHaz = REAL(varlhaz);
  double *VarLCum = REAL(varlcum);

  int nnph = lnph/lx;
  int nfix = lfix/lx;
  int nbase = Deg[1]-5;
  int leB = Deg[1]-1;
  int firstK = Deg[2];
  int i, j, z, t2, t3, tcz;
  double tempL, InvtempL;
  double TempD[6];

  double *MyGradLH = (double *)R_alloc(npar,sizeof(double));
  double *MyGradLC = (double *)R_alloc(npar,sizeof(double));
  double *MyParam = (double *)R_alloc(nbase,sizeof(double));
  double *MyBasisB = (double *)R_alloc(leB,sizeof(double));
  double *Res = (double *)R_alloc(nbase,sizeof(double));
  double *tempLvec = (double *)R_alloc(nbase,sizeof(double));

  double **GradLHaz = dmatrix(REAL(gradlhaz), A1, A2);
  double **GradLCum = dmatrix(REAL(gradlcum), A1, A2);

  for (z=0; z<lx; z++){

    t3 = nfix*z;
    for (i=0; i<nfix; i++){
      MyGradLH[i] = FixObs[i+t3];
      MyGradLC[i] = FixObs[i+t3];
    }

    tcz = TimeCat[z];
    VarLHaz[z] = 0;
    VarLCum[z] = 0;
    tempL = 0;

    t2 = z*nnph;
    for (i=0; i<nbase; i++){
      MyParam[i] = ParamT[i+nfix];
      tempLvec[i] = 0;
      for (j=1; j<nnph; j++){
	MyParam[i] += ParamT[j*nbase+i+nfix]*Nph[j+t2];
      }
    }

    for (i=0; i<tcz; i++){
      tempL += IntDNSpl(IntK[i], IntK[i+1], &TotK[i], &MatK[4*i], NsAdj1, NsAdj2, MyBasisB, TempD, MyParam, N, lW, lleg, leB, nbase, (i+firstK), tempLvec, Res);
    }
    tempL += IntDNSpl(IntK[tcz], X[z], &TotK[tcz], &MatK[4*tcz], NsAdj1, NsAdj2, MyBasisB, TempD, MyParam, N, lW, lleg, leB, nbase, (tcz+firstK), tempLvec, Res);
    InvtempL = 1/tempL;
    DeltaNSpl(X[z], &TotK[tcz], &MatK[4*tcz], NsAdj1, NsAdj2, MyBasisB, TempD, MyParam, leB, nbase, (tcz+firstK), Res);

    for (i=0; i<nnph; i++){
      for (j=0; j<nbase; j++){
	MyGradLH[nfix + i*nbase+j] = Res[j]*Nph[i+t2];
	MyGradLC[nfix + i*nbase+j] = tempLvec[j]*Nph[i+t2]*InvtempL;
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
