/*********************************************/
/* Routine for estimating the gradient and   */
/* the Hessian of the log-hazard and of the  */
/* cumulative hazard (aggregated)            */
/* (log-hazard described by a B-spline,      */
/* 1 time, w/ expected hazard)               */
/* Author: H. Charvat                        */
/* Last modified: 2020/09/05                 */
/* Part of the mexhaz 2.0 package            */
/*********************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include "SplineFunc.h"

SEXP HGHAggr_BsRx(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP lambdaobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk)
{
  SEXP loghaz, hazcum, test, gradlhaz, gradcum, hesslhaz, hesscum, rlist, rlistnames;
  int lx = length(x);
  int lnph = length(nph);
  int lfix = length(fixobs);
  int ltotk = length(totk);
  int lleg = length(n);
  int nclust = length(nbyclust);
  int npar = length(param)+length(paramf);
  int nhess = 0.5*npar*(npar+1);
  
  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(timecat = coerceVector(timecat,INTSXP));
  PROTECT(fixobs = coerceVector(fixobs,REALSXP));
  PROTECT(statobs = coerceVector(statobs,INTSXP));
  PROTECT(lambdaobs = coerceVector(lambdaobs,REALSXP));
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
  PROTECT(gradlhaz = allocVector(REALSXP,nclust*npar));
  PROTECT(gradcum = allocVector(REALSXP,nclust*npar));
  PROTECT(hesslhaz = allocVector(REALSXP,nclust*nhess));
  PROTECT(hesscum = allocVector(REALSXP,nclust*nhess));
  PROTECT(test = allocVector(LGLSXP,1));
  int nprotect = 21;

  double *X = REAL(x);
  double *Nph = REAL(nph);
  int *TimeCat = INTEGER(timecat);
  double *FixObs = REAL(fixobs);
  int *StatObs = INTEGER(statobs);
  double *LambdaObs = REAL(lambdaobs);
  int *NByClust = INTEGER(nbyclust);
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
  int Cst2;
  int Cst3 = nfix*npar-0.5*nfix*(nfix-1);
  int i, j, k, l, nc, m, t2, t3, tcz, cc;
  double tempF, tempLH, tempH, tempL, tempDH, tempTotH, tempGH, tempHH, expTF, NexpTF, NiNeTF, NlNeTF;
  double Total = 0;
  int z = 0;
  double *TotKPos = &TotK[Deg];

  double *MyParam = (double *)R_alloc(nbase+1,sizeof(double));
  double *TempD = (double *)R_alloc(Cst1,sizeof(double));
  double *Res = (double *)R_alloc(nbase+1,sizeof(double));
  double *tempLvec = (double *)R_alloc(nbase,sizeof(double));
  double *tempHess = (double *)R_alloc(nbase*nbase,sizeof(double));
  double *tempGLH = (double *)R_alloc(npar,sizeof(double));
  double *tempGC = (double *)R_alloc(npar,sizeof(double));
  
  double **GradLHaz = dmatrix(REAL(gradlhaz), nclust, npar);
  double **GradCum = dmatrix(REAL(gradcum), nclust, npar);
  double **HessLHaz = dmatrix(REAL(hesslhaz), nclust, nhess);
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
      
      tcz = TimeCat[z];
      tempL = 0;
      
      t2 = z*nnph;
      for (i=0; i<nbase; i++){
	MyParam[i+1] = Param[i];
	tempLvec[i] = 0;
	for (j=1; j<nnph; j++){
	  MyParam[i+1] += Param[j*nbase+i]*Nph[j+t2];
	}
	for (k=0; k<nbase; k++){
	  tempHess[k+i*nbase] = 0;
	}
      }
      
      for (i=0; i<tcz; i++){
	tempL += IntDSpline23H((*Fpt), TotKPos[i-1], TotKPos[i], &TotK[i], &MatK[Cst2*i], TempD, &MyParam[i], N, lW, lleg, nbase, i, tempLvec, tempHess, Res);
      }
      tempL += IntDSpline23H((*Fpt), TotKPos[tcz-1], X[z], &TotK[tcz], &MatK[Cst2*tcz], TempD, &MyParam[tcz], N, lW, lleg, nbase, tcz, tempLvec, tempHess, Res);
      tempLH = (*Fpt)(X[z], &TotK[tcz], &MatK[Cst2*tcz], TempD, &MyParam[tcz], nbase, tcz, Res);
          
      Total += tempL + tempLH + tempF;
      tempH = exp(tempLH + tempF);
      tempTotH = tempH+LambdaObs[z];
      tempDH = StatObs[z]*tempH;
      tempGH = tempDH/tempTotH;
      tempHH = tempDH*LambdaObs[z]/pow(tempTotH,2);
      LogHaz[nc] += StatObs[z]*log(tempTotH);
      HazCum[nc] += tempL*expTF;
      
      for (i=0; i<nfix; i++){
	tempGLH[i] = FixObs[i+t3];
	GradLHaz[i][nc] += tempGLH[i]*tempGH;
	tempGC[i] = FixObs[i+t3]*tempL*expTF;
	GradCum[i][nc] += tempGC[i];
      }

      cc = Cst3;
      for (i=0; i<nnph; i++){
	NexpTF = Nph[i+t2]*expTF;
	NiNeTF = Nph[i+t2]*NexpTF;
	for (j=0; j<nbase; j++){
	  tempGLH[nfix + i*nbase+j] = Res[j+1]*Nph[i+t2];
	  GradLHaz[nfix + i*nbase+j][nc] += tempGLH[nfix + i*nbase+j]*tempGH;
	  tempGC[nfix + i*nbase+j] = tempLvec[j]*NexpTF;
	  GradCum[nfix + i*nbase+j][nc] += tempGC[nfix + i*nbase+j];
	  for (k=j; k<nbase; k++){
	    HessCum[cc][nc] += tempHess[k+j*nbase]*NiNeTF;
	    cc++;
	  }
	  for (l=i+1; l<nnph; l++){
	    NlNeTF = Nph[l+t2]*NexpTF;
	    for (k=0; k<nbase; k++){
	      HessCum[cc][nc] += tempHess[k+j*nbase]*NlNeTF;
	      cc++;
	    }
	  }	  
	}
      }
            
      cc = 0;
      for (i=0; i<nfix; i++){
	for (j=i; j<npar; j++){
	  HessCum[cc][nc] += FixObs[i+t3]*tempGC[j];
	  cc++;
	}
      }

      cc = 0;
      for (i=0; i<npar; i++){
	for (j=i; j<npar; j++){
	  HessLHaz[cc][nc] += tempGLH[i]*tempGLH[j]*tempHH;
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
