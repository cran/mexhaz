/********************************************/
/* Routines for estimating the marginal     */
/* log-likelihood by adaptive Gauss-Hermite */
/* quadrature                               */
/* in the presence of left truncation       */
/* Author: H. Charvat                       */
/* Last modified: 2020/07/06                */
/* Part of the mexhaz 1.8 package           */
/********************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include "OptFunc.h"

SEXP FrailtyAdaptL(SEXP nodes, SEXP nodessquare, SEXP logweights, SEXP clust, SEXP clustd, SEXP expect, SEXP betal, SEXP betaL0, SEXP betaL, SEXP A0, SEXP A, SEXP var, SEXP mh0, SEXP muhatcond)
{
  SEXP vecMuH, vecSigH, vecCstAdj0, vecCstAdj, loglik, rlist, rlistnames;
  int Npoint = length(nodes);
  int Nclust = length(clust);

  PROTECT(nodes = coerceVector(nodes,REALSXP));
  PROTECT(nodessquare = coerceVector(nodessquare,REALSXP));
  PROTECT(logweights = coerceVector(logweights,REALSXP));
  PROTECT(expect = coerceVector(expect,REALSXP));
  PROTECT(betal = coerceVector(betal,REALSXP));
  PROTECT(betaL0 = coerceVector(betaL0,REALSXP));
  PROTECT(betaL = coerceVector(betaL,REALSXP));
  PROTECT(clust = coerceVector(clust,INTSXP));
  PROTECT(clustd = coerceVector(clustd,INTSXP));
  PROTECT(A0 = coerceVector(A0,REALSXP));
  PROTECT(A = coerceVector(A,REALSXP));
  PROTECT(var = coerceVector(var,REALSXP));
  PROTECT(mh0 = coerceVector(mh0,REALSXP));
  PROTECT(muhatcond = coerceVector(muhatcond,INTSXP));
  PROTECT(vecMuH = allocVector(REALSXP,Nclust));
  PROTECT(vecSigH = allocVector(REALSXP,Nclust));
  PROTECT(vecCstAdj0 = allocVector(REALSXP,Nclust));
  PROTECT(vecCstAdj = allocVector(REALSXP,Nclust));
  PROTECT(loglik = allocVector(REALSXP,1));
  int nprotect = 19;

  double *Nodes = REAL(nodes);
  double *Nodessquare = REAL(nodessquare);
  double *Logweights = REAL(logweights);
  double *Expect = REAL(expect);
  double *Betal = REAL(betal);
  double *BetaL0 = REAL(betaL0);
  double *BetaL = REAL(betaL);
  int *Clust = INTEGER(clust);
  int *Clustd = INTEGER(clustd);
  double *AA0 = REAL(A0);
  double *AA = REAL(A);
  double Var = REAL(var)[0];
  double *Muhat0 = REAL(mh0);
  int Muhatcond = INTEGER(muhatcond)[0];
  double *VecMuH = REAL(vecMuH);
  double *VecSigH = REAL(vecSigH);
  double *VecCstAdj0 = REAL(vecCstAdj0);
  double *VecCstAdj = REAL(vecCstAdj);

  int k, z, nclust, lenclust, lenclustd;
  double lbetaL0, lbetaL, muhat, sigmahat0, sigmahat, logsigmahat0, logsigmahat, ddmli, ddmli0, resclust0, resclust, xstar0, xstar;
  int d = 0, dd = 0;
  int essai;
  double casup, cainf;
  double cstadd = 0.5*log(Var) + M_LN_SQRT_PI;
  double LogLik = 0;

  double *logresll = (double *)R_alloc(Npoint,sizeof(double));
  double *logdenom = (double *)R_alloc(Npoint,sizeof(double));

  for (nclust=0; nclust<Nclust; nclust++) {
    lenclust = Clust[nclust];
    lenclustd = Clustd[nclust];
    lbetaL0 = 0;
    lbetaL = 0;
    for (z = 0; z<lenclust; z++){
      lbetaL0 += BetaL0[d+z];
      lbetaL += BetaL[d+z];
    }
    lbetaL = log(lbetaL);
    lbetaL0 = log(lbetaL0);
   
    muhat = ZeroDMLI(-100, 1000, 1e-7, lenclustd, &Expect[dd], &Betal[dd], lbetaL, Var);
    VecMuH[nclust] = muhat;
    if (Muhatcond!=1){
      ddmli = DDMLI(muhat, lenclustd, &Expect[dd], &Betal[dd], lbetaL, Var);
      sigmahat = 1/sqrt(ddmli);
      VecSigH[nclust] = sigmahat;
      if (Muhatcond!=2){

 	/* Denominator */
	casup = 1000.0;
	cainf = -1000.0;
	essai = 0;
	ddmli0 = DDMLI0(Muhat0[nclust], lbetaL0, Var);
	sigmahat0 = 1/sqrt(ddmli0);
	logsigmahat0 = 0.5*log(ddmli0) + cstadd;
	for(k=0; k<Npoint; k++) {
	  xstar0 = Muhat0[nclust] + M_SQRT2*sigmahat0*Nodes[k];
	  logdenom[k] = Logweights[k] - logsigmahat0 + Nodessquare[k] - pow(xstar0,2)/(2*Var) + LogProd0(xstar0, lbetaL0);
	}
	resclust0 = LLGHQClust(Npoint, logdenom, AA0[nclust], lenclust);
	while ((isinf(fabs(resclust0)) || isnan(resclust0)) && (essai<25)){
	  if (isnan(resclust0) || (resclust0==-1.0/0.0)){
	    casup = AA0[nclust];
	  }
	  else {
	    cainf = AA0[nclust];
	  }
	  AA0[nclust] = (cainf+casup)/2;
	  resclust0 = LLGHQClust(Npoint, logdenom, AA0[nclust], lenclust);
	  essai += 1;
	}
	VecCstAdj0[nclust] = AA0[nclust];

	/* Numerator  */
	casup = 1000.0;
	cainf = -1000.0;
	essai = 0;
	logsigmahat = 0.5*log(ddmli) + cstadd;
	for(k=0; k<Npoint; k++) {
	  xstar = muhat + M_SQRT2*sigmahat*Nodes[k];
	  logresll[k] = Logweights[k] - logsigmahat + Nodessquare[k] - pow(xstar,2)/(2*Var) + LogProd(xstar, lenclustd, &Expect[dd], &Betal[dd], lbetaL);
	}
	resclust = LLGHQClust(Npoint, logresll, AA[nclust], lenclust);
	while ((isinf(fabs(resclust)) || isnan(resclust)) && (essai<25)){
	  if (isnan(resclust) || (resclust==-1.0/0.0)){
	    casup = AA[nclust];
	  }
	  else {
	    cainf = AA[nclust];
	  }
	  AA[nclust] = (cainf+casup)/2;
	  resclust = LLGHQClust(Npoint, logresll, AA[nclust], lenclust);
	  essai += 1;
	}
	VecCstAdj[nclust] = AA[nclust];
	
	LogLik += resclust-resclust0;
      }
    }   
    d = d + lenclust;
    dd = dd + lenclustd;
  }
  REAL(loglik)[0] = LogLik;

 /* assemble the return objects as a list */
  PROTECT(rlist= allocVector(VECSXP, 5));
  SET_VECTOR_ELT(rlist, 0, vecMuH);
  SET_VECTOR_ELT(rlist, 1, vecSigH);
  SET_VECTOR_ELT(rlist, 2, vecCstAdj0);
  SET_VECTOR_ELT(rlist, 3, vecCstAdj);
  SET_VECTOR_ELT(rlist, 4, loglik);

  /* add names to the list elements */
  PROTECT(rlistnames = allocVector(STRSXP, 5));
  SET_STRING_ELT(rlistnames, 0, mkChar("MuHat"));
  SET_STRING_ELT(rlistnames, 1, mkChar("SigmaHat"));
  SET_STRING_ELT(rlistnames, 2, mkChar("CstAdj0"));
  SET_STRING_ELT(rlistnames, 3, mkChar("CstAdj"));
  SET_STRING_ELT(rlistnames, 4, mkChar("LogLik"));
  setAttrib(rlist, R_NamesSymbol, rlistnames);

  UNPROTECT(nprotect+2);
  return rlist;
}
