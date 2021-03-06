/************************************************/
/* Routine for estimating the cumulative hazard */
/* (log-hazard described by linear B-spline)    */
/* Author: H. Charvat                           */
/* Last modified: 2019/11/26                    */
/* Part of the mexhaz 1.7 package               */
/************************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>

SEXP HazardBs1R(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP matk, SEXP totk)
{
  SEXP loghaz, hazcum0, hazcum, test, rlist, rlistnames;
  int lx = length(x);
  int lnph = length(nph);
  int ltotk = length(totk);
  int lfix = length(fixobs);

  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(timecat = coerceVector(timecat,INTSXP));
  PROTECT(fixobs = coerceVector(fixobs,REALSXP));
  PROTECT(param = coerceVector(param,REALSXP));
  PROTECT(paramf = coerceVector(paramf,REALSXP));
  PROTECT(matk = coerceVector(matk,REALSXP));
  PROTECT(totk = coerceVector(totk,REALSXP));
  PROTECT(loghaz = allocVector(REALSXP,lx));
  PROTECT(hazcum0 = allocVector(REALSXP,1));
  PROTECT(hazcum = allocVector(REALSXP,lx));
  PROTECT(test = allocVector(LGLSXP,1));
  int nprotect = 12;

  double *X = REAL(x);
  double *Nph = REAL(nph);
  int *TimeCat = INTEGER(timecat);
  double *FixObs = REAL(fixobs);
  double *Param = REAL(param);
  double *ParamF = REAL(paramf);
  double *MatK = REAL(matk);
  double *TotK = REAL(totk);
  double *LogHaz = REAL(loghaz);
  double *HazCum = REAL(hazcum);

  int nnph = lnph/lx;
  int nfix = lfix/lx;
  int nbase = ltotk-1;
  int i, j, k, z, tcz, t1;
  double Temp, tempF, Beta1, Beta2;
  double Total = 0;

  double *MyParam = (double *)R_alloc((nbase+1),sizeof(double));

  MyParam[0] = 0;

  if (nnph==1){

    for (i=0; i<nbase; i++){
      MyParam[i+1] = Param[i];
    }

    for (z=0; z<lx; z++){

      tempF = ParamF[0];
      t1 = z*nfix;
      for (i=1; i<nfix; i++){
	tempF += FixObs[i+t1]*ParamF[i];
      }

      tcz = TimeCat[z];

      Beta1 = MyParam[tcz];
      Beta2 = MyParam[tcz+1];
      Temp = Beta2-Beta1;
      if (Temp!=0){
	LogHaz[z] = (1/MatK[tcz])*(Beta1*(TotK[tcz+1]-X[z])+Beta2*(X[z]-TotK[tcz]));
	HazCum[z] = (MatK[tcz]/Temp)*(exp(LogHaz[z])-exp(Beta1));
      }
      else {
	LogHaz[z] = Beta1;
	HazCum[z] = (X[z]-TotK[tcz])*exp(Beta1);
      }

      if (tcz>0){
	for (k=tcz; k>0; k--){
	  Beta1 = MyParam[k-1];
	  Beta2 = MyParam[k];
	  Temp = Beta2-Beta1;
	  if (Temp!=0) {
	    HazCum[z] += (MatK[k-1]/Temp)*(exp(Beta2)-exp(Beta1));
	  }
	  else {
	    HazCum[z] += MatK[k-1]*exp(Beta1);
	  }
	}
      }
      Total += HazCum[z] + LogHaz[z] + tempF;
      LogHaz[z] += tempF;
      HazCum[z] *= exp(tempF);
    }

  }
  else {

    int t2;
    for (z=0; z<lx; z++){

      tempF = ParamF[0];
      t1 = z*nfix;
      for (i=1; i<nfix; i++){
	tempF += FixObs[i+t1]*ParamF[i];
      }

      tcz = TimeCat[z];

      t2 = z*nnph;
      for (i=0; i<nbase; i++){
	MyParam[i+1] = Param[i];
	for (j=1; j<nnph; j++){
	  MyParam[i+1] += Param[j*nbase+i]*Nph[j+t2];
	}
      }

      Beta1 = MyParam[tcz];
      Beta2 = MyParam[tcz+1];
      Temp = Beta2-Beta1;
      if (Temp!=0){
	LogHaz[z] = (1/MatK[tcz])*(Beta1*(TotK[tcz+1]-X[z])+Beta2*(X[z]-TotK[tcz]));
	HazCum[z] = (MatK[tcz]/Temp)*(exp(LogHaz[z])-exp(Beta1));
      }
      else {
	LogHaz[z] = Beta1;
	HazCum[z] = (X[z]-TotK[tcz])*exp(Beta1);
      }

      for (k=tcz; k>0; k--){
	Beta1 = MyParam[k-1];
	Beta2 = MyParam[k];
	Temp = Beta2-Beta1;
	if (Temp!=0) {
	  HazCum[z] += (MatK[k-1]/Temp)*(exp(Beta2)-exp(Beta1));
	}
	else {
	  HazCum[z] += MatK[k-1]*exp(Beta1);
	}
      }
      Total += HazCum[z] + LogHaz[z] + tempF;
      LogHaz[z] += tempF;
      HazCum[z] *= exp(tempF);
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
