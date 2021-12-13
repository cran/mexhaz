/****************************************************/
/* Part of the NPH test based on counting processes */
/* from Danieli et al., Biostatistics, 2017         */
/* Author: H. Charvat                               */
/* Created: 2017/03/03                              */
/* Last modified: 2021/09/05                        */
/* Part of the mexhaz 2.0 package                   */
/****************************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include "SplineFunc.h"

SEXP GaussProcNPH(SEXP vecnb, SEXP pDtau, SEXP lrowHXS, SEXP maxSXt, SEXP cst, SEXP idxpDx, SEXP keep)
{
  SEXP gproc, gptemp, pval, rlist, rlistnames;

  int npdt = length(pDtau);
  int ntime = length(idxpDx);

  PROTECT(vecnb = coerceVector(vecnb,INTSXP));
  PROTECT(pDtau = coerceVector(pDtau,REALSXP));
  PROTECT(lrowHXS = coerceVector(lrowHXS,VECSXP));
  PROTECT(maxSXt = coerceVector(maxSXt,REALSXP));
  PROTECT(cst = coerceVector(cst,REALSXP));
  PROTECT(idxpDx = coerceVector(idxpDx,INTSXP));
  PROTECT(keep = coerceVector(keep,INTSXP));
  int nprotect = 7;

  int nsim = INTEGER(vecnb)[0];
  int nobs = INTEGER(vecnb)[1];
  int ntest = INTEGER(vecnb)[2];
  double *PDtau = REAL(pDtau);
  double *MaxSXt = REAL(maxSXt);
  double Cst = REAL(cst)[0];
  int *IdxpDx = INTEGER(idxpDx);
  int Keep = INTEGER(keep)[0];

  int nvar = npdt/nobs;
  int nbase = nvar-ntest;

  PROTECT(pval = allocVector(REALSXP,ntest));
  double *Pval = REAL(pval);
  nprotect += 1;

  int i, j, k, s, t, t1, t2;

  double tempGP, ftGP;
  double Matmax = 0;

  double *Dtau = (double *)R_alloc(nvar,sizeof(double));
  double *G = (double *)R_alloc(nobs,sizeof(double));
  double *DXt = (double *)R_alloc(nobs,sizeof(double));
  double *mDXt = (double *)R_alloc(ntime,sizeof(double));

  PROTECT(gproc = allocVector(VECSXP,ntest));
  nprotect += 1;

  for (t=0; t<ntest; t++){
    
    PROTECT(gptemp = allocVector(REALSXP,ntime*(2+Keep)));
    double **Gptemp = dmatrix(REAL(gptemp), ntime, 2+Keep);

    Pval[t] = 0;
    for (k=0; k<ntime; k++){
      Gptemp[0][k] = 0;
      Gptemp[1][k] = 0;
    }
  
    GetRNGstate();
    for (s=0; s<nsim; s++){
      
      for (j=0; j<nvar; j++){
	Dtau[j] = 0;
      }
      
      for (i=0; i<nobs; i++){
	G[i] = norm_rand();
	t1 = i*nvar;
	for (j=0; j<nvar; j++){
	  Dtau[j] += G[i]*PDtau[j+t1];
	}
	DXt[i] = Dtau[nbase+t];
      }
      
      double *RowHXS = REAL(VECTOR_ELT(lrowHXS,t));
      
      Matmax = 0;
      for (k=0; k<ntime; k++){
	mDXt[k] = 0;
	t2 = k*nvar;
	for (j=0; j<nvar; j++){
	  mDXt[k] += RowHXS[j+t2]*Dtau[j];
	}
	tempGP = DXt[IdxpDx[k]] - mDXt[k];
	if ((ftGP = fabs(tempGP))>Matmax){
	  Matmax = ftGP;
	}
	tempGP *= Cst;
	if (s<Keep){
	  Gptemp[2+s][k] = tempGP;
	}
	if (tempGP<=Gptemp[0][k]){
	  Gptemp[0][k] = tempGP;
	}
	if (tempGP>=Gptemp[1][k]){
	  Gptemp[1][k] = tempGP;
	}
      }
      if (Matmax>=MaxSXt[t]){
	Pval[t] += 1;
      }
    }
    PutRNGstate();

    SET_VECTOR_ELT(gproc, t, gptemp);
    Pval[t] = Pval[t]/nsim;
    UNPROTECT(1);
  }
  
  /* assemble the return objects as a list */
  PROTECT(rlist = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(rlist, 0, gproc);
  SET_VECTOR_ELT(rlist, 1, pval);
  
  /* add names to the list elements */
  PROTECT(rlistnames = allocVector(STRSXP, 2));
  SET_STRING_ELT(rlistnames, 0, mkChar("Gproc"));
  SET_STRING_ELT(rlistnames, 1, mkChar("Pval"));

  setAttrib(rlist, R_NamesSymbol, rlistnames);

  UNPROTECT(nprotect+2);
  return rlist;
}
