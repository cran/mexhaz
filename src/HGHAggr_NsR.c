/*********************************************/
/* Routine for estimating the gradient and   */
/* the Hessian of the log-hazard and of the  */
/* cumulative hazard (aggregated)            */
/* (log-hazard described by a natural cubic  */
/* spline, 1 time, no expected hazard)       */
/* Author: H. Charvat                        */
/* Last modified: 2020/09/05                 */
/* Part of the mexhaz 2.0 package            */
/*********************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include "SplineFunc.h"

SEXP HGHAggr_NsR(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP intk, SEXP nsadj1, SEXP nsadj2)
{
  SEXP loghaz, hazcum, test, gradlhaz, gradcum, hesslhaz, hesscum, rlist, rlistnames;
  int lx = length(x);
  int lnph = length(nph);
  int lfix = length(fixobs);
  int lleg = length(n);
  int nclust = length(nbyclust);
  int npar = length(param)+length(paramf);
  int nhess = 0.5*npar*(npar+1);
  
  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
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
  PROTECT(intk = coerceVector(intk,REALSXP));
  PROTECT(nsadj1 = coerceVector(nsadj1,REALSXP));
  PROTECT(nsadj2 = coerceVector(nsadj2,REALSXP));
  PROTECT(loghaz = allocVector(REALSXP,nclust));
  PROTECT(hazcum = allocVector(REALSXP,nclust));
  PROTECT(gradlhaz = allocVector(REALSXP,nclust*npar));
  PROTECT(gradcum = allocVector(REALSXP,nclust*npar));
  PROTECT(hesslhaz = allocVector(REALSXP,1));
  PROTECT(hesscum = allocVector(REALSXP,nclust*nhess));
  PROTECT(test = allocVector(LGLSXP,1));
  int nprotect = 23;

  double *X = REAL(x);
  double *Nph = REAL(nph);
  int *TimeCat = INTEGER(timecat);
  double *FixObs = REAL(fixobs);
  int *StatObs = INTEGER(statobs);
  int *NByClust = INTEGER(nbyclust);
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
  double *HazCum = REAL(hazcum);

  int nnph = lnph/lx;
  int nfix = lfix/lx;
  int nbase = Deg[1]-5;
  int leB = Deg[1]-1;
  int firstK = Deg[2];
  int Cst3 = nfix*npar-0.5*nfix*(nfix-1);
  int i, j, k, l, nc, m, t2, t3, tcz, cc;
  double tempF, tempH, tempL, expTF, NexpTF, NiNeTF, NlNeTF;
  double Total = 0;
  int z = 0;
  double TempD[6];

  double *MyParam = (double *)R_alloc(nbase,sizeof(double));
  double *MyBasisB = (double *)R_alloc(leB,sizeof(double));
  double *Res = (double *)R_alloc(nbase,sizeof(double));
  double *tempLvec = (double *)R_alloc(nbase,sizeof(double));
  double *tempHess = (double *)R_alloc(nbase*nbase,sizeof(double));
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
	MyParam[i] = Param[i];
	tempLvec[i] = 0;
	for (j=1; j<nnph; j++){
	  MyParam[i] += Param[j*nbase+i]*Nph[j+t2];
	}
	for (k=0; k<nbase; k++){
	  tempHess[k+i*nbase] = 0;
	}
      }
  
      for (i=0; i<tcz; i++){
	tempL += IntDNSplH(IntK[i], IntK[i+1], &TotK[i], &MatK[4*i], NsAdj1, NsAdj2, MyBasisB, TempD, MyParam, N, lW, lleg, leB, nbase, (i+firstK), tempLvec, tempHess, Res);
    }
      tempL += IntDNSplH(IntK[tcz], X[z], &TotK[tcz], &MatK[4*tcz], NsAdj1, NsAdj2, MyBasisB, TempD, MyParam, N, lW, lleg, leB, nbase, (tcz+firstK), tempLvec, tempHess, Res);
      tempH = DeltaNSpl(X[z], &TotK[tcz], &MatK[4*tcz], NsAdj1, NsAdj2, MyBasisB, TempD, MyParam, leB, nbase, (tcz+firstK), Res);
            
      Total += tempL + tempH + tempF;
      LogHaz[nc] += StatObs[z]*(tempH + tempF);
      HazCum[nc] += tempL*expTF;
      
      for (i=0; i<nfix; i++){
	GradLHaz[i][nc] += StatObs[z]*FixObs[i+t3];
	tempGC[i] = FixObs[i+t3]*tempL*expTF;
	GradCum[i][nc] += tempGC[i];
      }
      
      cc = Cst3;
      for (i=0; i<nnph; i++){
	NexpTF = Nph[i+t2]*expTF;
	NiNeTF = Nph[i+t2]*NexpTF;
	for (j=0; j<nbase; j++){
	  GradLHaz[nfix + i*nbase+j][nc] += Res[j]*Nph[i+t2]*StatObs[z];
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
