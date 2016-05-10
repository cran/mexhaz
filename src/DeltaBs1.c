/*********************************************/
/* Routine for estimating the variance of    */
/* the log-hazard and log-cumulative hazard  */
/* by the Delta Method                       */
/* (log-hazard described by linear B-spline) */
/* Author: H. Charvat                        */
/* Last modified: 2016/05/10                 */
/* Part of the mexhaz package                */
/*********************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>

SEXP DeltaBs1(SEXP x, SEXP nph, SEXP param, SEXP others, SEXP leint, SEXP whint, SEXP knots, SEXP varcov)
{
  SEXP result;
  int lx = length(x);
  int lnph = length(nph);
  int loth = length(others);
  int npar = length(param);
  int lk = length(knots);
  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(param = coerceVector(param,REALSXP));
  PROTECT(others = coerceVector(others,REALSXP));
  PROTECT(leint = coerceVector(leint,REALSXP));
  PROTECT(whint = coerceVector(whint,INTSXP));
  PROTECT(knots = coerceVector(knots,REALSXP));
  PROTECT(varcov = coerceVector(varcov,REALSXP));
  PROTECT(result = allocVector(REALSXP,(2*lx)));
  double *X = REAL(x);
  double *Nph = REAL(nph);
  double *Param = REAL(param);
  double *Others = REAL(others);
  double *Leint = REAL(leint);
  int *Whint = INTEGER(whint);
  double *Knots = REAL(knots);
  double *Varcov = REAL(varcov);
  double *Result = REAL(result);
  int nnph = lnph/lx;
  int nind = loth/lx;
  int lsdk = lk-1;
  double templ, tempL, InvtempL;
  double *MyGradL = malloc(npar*sizeof(double));
  double *MyGradS = malloc(npar*sizeof(double));
  double *MyParam = malloc((lsdk+1)*sizeof(double));
  double *Res = malloc((lsdk+1)*sizeof(double));
  double *tempLvec = malloc((lsdk+1)*sizeof(double));
  int i, ii, j, k, z, t2, t3, Which;
  double Temp, Beta1, Beta2;
  double Cst1 = 0;
  double Cst2 = 0;

  Res[0] = 0;
  MyParam[0] = 0;

  for (z=0; z<lx; z++){

    t3 = nind*z;
    for (i=0; i<nind; i++){
      MyGradL[i] = Others[i+t3];
      MyGradS[i] = Others[i+t3];
    }

    Which = Whint[z];
    Result[z] = 0;
    Result[lx+z] = 0;
    tempLvec[0] = 0;

    t2 = z*nnph;
    for (i=0; i<lsdk; i++){
      MyParam[i+1] = Param[i+nind];
      Res[i+1] = 0;
      tempLvec[i+1] = 0;
      ii = 1;
      while (ii<nnph){
	MyParam[i+1] += Param[ii*lsdk+i+nind]*Nph[ii+t2];
	ii++;
      }
    }

    // Calculation of lambda, Lambda and necessary integrals //
    Beta1 = MyParam[Which];
    Beta2 = MyParam[Which+1];
    Res[Which] = (Knots[Which+1]-X[z])/Leint[Which];
    Res[Which+1] = (X[z]-Knots[Which])/Leint[Which];
    Temp = Beta2-Beta1;
    if (Temp!=0){
      Cst1 = Leint[Which]/Temp;
      Cst2 = 1/Temp;
      templ = Beta1*Res[Which]+Beta2*Res[Which+1];
      tempL = Cst1*(exp(templ)-exp(Beta1));
      tempLvec[Which] = Cst1*((Res[Which]+Cst2)*exp(templ)-(1+Cst2)*exp(Beta1));
      tempLvec[Which+1] = Cst1*((Res[Which+1]-Cst2)*exp(templ)+Cst2*exp(Beta1));
    }
    else {
      templ = Beta1;
      Cst1 = 0.5*exp(templ)*Res[Which+1]*Leint[Which];
      tempL = (X[z]-Knots[Which])*exp(templ);
      tempLvec[Which] = Cst1*(1+Res[Which]);
      tempLvec[Which+1] = Cst1*Res[Which+1];
    }

    if (Which>0){
      for (k=Which; k>0; k--){
	Beta1 = MyParam[k-1];
	Beta2 = MyParam[k];
	Temp = Beta2-Beta1;
	if (Temp!=0) {
	  Cst1 = Leint[k-1]/Temp;
	  Cst2 = 1/Temp;
	  tempL += Cst1*(exp(Beta2)-exp(Beta1));
	  tempLvec[k] += Cst1*((1-Cst2)*exp(Beta2)+Cst2*exp(Beta1));
	  tempLvec[k-1] += Cst1*(Cst2*exp(Beta2)-(1+Cst2)*exp(Beta1));
	}
	else {
	  Cst1 = Leint[k-1]*exp(Beta1);
	  tempL += Cst1;
	  tempLvec[k] += 0.5*Cst1;
	  tempLvec[k-1] += 0.5*Cst1;
	}
      }
    }
    InvtempL = 1/tempL;

    ii = 0;
    while (ii<nnph){
      for (j=0; j<lsdk; j++){
	MyGradL[nind + ii*lsdk+j] = Res[j+1]*Nph[ii+t2];
	MyGradS[nind + ii*lsdk+j] = tempLvec[j+1]*Nph[ii+t2]*InvtempL;
      }
      ii++;
    }

    for (i=0; i<npar; i++){
      for (j=0; j<npar; j++){
	Result[z] += MyGradL[i]*Varcov[j+npar*i]*MyGradL[j];
	Result[lx+z] += MyGradS[i]*Varcov[j+npar*i]*MyGradS[j];
      }
    }

  }

  UNPROTECT(9);
  free(MyGradL);
  free(MyGradS);
  free(MyParam);
  free(Res);
  free(tempLvec);
  return result;
}
