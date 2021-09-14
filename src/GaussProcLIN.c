/****************************************************/
/* Part of the LIN test based on counting processes */
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

SEXP GaussProcLIN(SEXP vecnb, SEXP event, SEXP pDtau, SEXP lrowHXS, SEXP maxSXt, SEXP cst, SEXP idxpDx, SEXP ordDx, SEXP nbval, SEXP keep)
{
  SEXP gproc, gptemp, pval, rlist, rlistnames;

  int npdt = length(pDtau);
  int nobs = length(event);
  int ntest = length(nbval);
  
  PROTECT(vecnb = coerceVector(vecnb,INTSXP));
  PROTECT(event = coerceVector(event,INTSXP));
  PROTECT(pDtau = coerceVector(pDtau,REALSXP));
  PROTECT(lrowHXS = coerceVector(lrowHXS,VECSXP));
  PROTECT(maxSXt = coerceVector(maxSXt,REALSXP));
  PROTECT(cst = coerceVector(cst,REALSXP));
  PROTECT(idxpDx = coerceVector(idxpDx,VECSXP));
  PROTECT(ordDx = coerceVector(ordDx,VECSXP));
  PROTECT(nbval = coerceVector(nbval,INTSXP));
  PROTECT(keep = coerceVector(keep,INTSXP));
  int nprotect = 10;

  int nsim = INTEGER(vecnb)[0];
  int *Event = INTEGER(event);
  double *PDtau = REAL(pDtau);
  double *MaxSXt = REAL(maxSXt);
  double Cst = REAL(cst)[0];
  int *NbVal = INTEGER(nbval);
  int Keep = INTEGER(keep)[0];

  int nvar = npdt/nobs;

  PROTECT(pval = allocVector(REALSXP,ntest));
  double *Pval = REAL(pval);
  nprotect += 1;

  int i, j, k, s, t, t1, t2, ntime;

  double tempGP, ftGP;
  double Matmax = 0;

  double *Dtau = (double *)R_alloc(nvar,sizeof(double));
  double *G = (double *)R_alloc(nobs,sizeof(double));
  double *DXt = (double *)R_alloc(nobs,sizeof(double));
  double *DXtp = (double *)R_alloc(nobs,sizeof(double));

  PROTECT(gproc = allocVector(VECSXP,2));
  nprotect += 1;
  
  for (t=0; t<ntest; t++){
    
    int *IdxpDx = INTEGER(VECTOR_ELT(idxpDx,t));
    int *OrdDx = INTEGER(VECTOR_ELT(ordDx,t));
    double *RowHXS = REAL(VECTOR_ELT(lrowHXS,t));
    ntime = NbVal[t];
    
    PROTECT(gptemp = allocVector(REALSXP,ntime*(2+Keep)));
    double **Gptemp = dmatrix(REAL(gptemp), ntime, 2+Keep);
    double *mDXt = (double *)R_alloc(ntime,sizeof(double));
    
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
	DXt[i] = G[i]*Event[i];
      }
      
      DXtp[0] = DXt[OrdDx[0]];
      for (i=1; i<nobs; i++){
	DXtp[i] = DXtp[i-1]+DXt[OrdDx[i]];
      }
      
      Matmax = 0;
      for (k=0; k<ntime; k++){
	mDXt[k] = 0;
	t2 = k*nvar;
	for (j=0; j<nvar; j++){
	  mDXt[k] += RowHXS[j+t2]*Dtau[j];
	}
	tempGP = DXtp[IdxpDx[k]] - mDXt[k];
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
    
    Pval[t] = Pval[t]/nsim;
    SET_VECTOR_ELT(gproc, t, gptemp);
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
