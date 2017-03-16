/********************************************/
/* Routines for estimating the marginal     */
/* log-likelihood by adaptive Gauss-Hermite */
/* quadrature                               */
/* Author: H. Charvat                       */
/* Last modified: 2016/12/27                */
/* Part of the mexhaz 1.3 package           */
/********************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>

double LogProd(double x, int lenclust, double *expect, double *betal, double LSEbetaL);
double DLogProd(double x, int lenclust, double *expect, double *betal, double LSEbetaL);
double DMLI(double x, int lenclust, double *expect, double *betal, double LSEbetaL, double var);
double ZeroDMLI(double xinf, double xmax, double tol, int lenclust, double *expect, double *betal, double LSEbetaL, double var);
double DDLogProd(double x, int lenclust, double *expect, double *betal, double LSEbetaL);
double DDMLI(double x, int lenclust, double *expect, double *betal, double LSEbetaL, double var);
double LLGHQClust(int Npoint, double *logresll, double A, int clust);

SEXP FrailtyAdapt(SEXP nodes, SEXP nodessquare, SEXP logweights, SEXP clust, SEXP clustd, SEXP expect, SEXP betal, SEXP betaL, SEXP A, SEXP var, SEXP muhatcond)
{
  SEXP vecMuH, vecSigH, vecCstAdj, loglik, rlist, rlistnames;
  int Npoint = length(nodes);
  int Nclust = length(clust);

  PROTECT(nodes = coerceVector(nodes,REALSXP));
  PROTECT(nodessquare = coerceVector(nodessquare,REALSXP));
  PROTECT(logweights = coerceVector(logweights,REALSXP));
  PROTECT(expect = coerceVector(expect,REALSXP));
  PROTECT(betal = coerceVector(betal,REALSXP));
  PROTECT(betaL = coerceVector(betaL,REALSXP));
  PROTECT(clust = coerceVector(clust,INTSXP));
  PROTECT(clustd = coerceVector(clustd,INTSXP));
  PROTECT(A = coerceVector(A,REALSXP));
  PROTECT(var = coerceVector(var,REALSXP));
  PROTECT(muhatcond = coerceVector(muhatcond,INTSXP));
  PROTECT(vecMuH = allocVector(REALSXP,Nclust));
  PROTECT(vecSigH = allocVector(REALSXP,Nclust));
  PROTECT(vecCstAdj = allocVector(REALSXP,Nclust));
  PROTECT(loglik = allocVector(REALSXP,1));
  int nprotect = 15;

  double *Nodes = REAL(nodes);
  double *Nodessquare = REAL(nodessquare);
  double *Logweights = REAL(logweights);
  double *Expect = REAL(expect);
  double *Betal = REAL(betal);
  double *BetaL = REAL(betaL);
  int *Clust = INTEGER(clust);
  int *Clustd = INTEGER(clustd);
  double *AA = REAL(A);
  double Var = REAL(var)[0];
  int Muhatcond = INTEGER(muhatcond)[0];
  double *VecMuH = REAL(vecMuH);
  double *VecSigH = REAL(vecSigH);
  double *VecCstAdj = REAL(vecCstAdj);

  int k, z, nclust, lenclust, lenclustd;
  double lbetaL, muhat, sigmahat, logsigmahat, ddmli, resclust, xstar;
  int d = 0, dd = 0;
  int essai;
  double casup, cainf;
  double cstadd = 0.5*log(Var) + M_LN_SQRT_PI;
  double LogLik = 0;

  double *logresll = (double *)R_alloc(Npoint,sizeof(double));

  for (nclust=0; nclust<Nclust; nclust++) {
    casup = 1000.0;
    cainf = -1000.0;
    essai = 0;
    lenclust = Clust[nclust];
    lenclustd = Clustd[nclust];
    lbetaL = 0;
    for (z = 0; z<lenclust; z++){
      lbetaL += exp(BetaL[d+z]);
    }
    lbetaL = log(lbetaL);
    muhat = ZeroDMLI(-100, 1000, 1e-7, lenclustd, &Expect[dd], &Betal[dd], lbetaL, Var);
    VecMuH[nclust] = muhat;
    if (Muhatcond!=1){
      ddmli = DDMLI(muhat, lenclustd, &Expect[dd], &Betal[dd], lbetaL, Var);
      sigmahat = 1/sqrt(ddmli);
      VecSigH[nclust] = sigmahat;
      if (Muhatcond!=2){
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
	LogLik += resclust;
      }
    }
    d = d + lenclust;
    dd = dd + lenclustd;
  }
  REAL(loglik)[0] = LogLik;

 /* assemble the return objects as a list */
  PROTECT(rlist= allocVector(VECSXP, 4));
  SET_VECTOR_ELT(rlist, 0, vecMuH);
  SET_VECTOR_ELT(rlist, 1, vecSigH);
  SET_VECTOR_ELT(rlist, 2, vecCstAdj);
  SET_VECTOR_ELT(rlist, 3, loglik);

  /* add names to the list elements */
  PROTECT(rlistnames = allocVector(STRSXP, 4));
  SET_STRING_ELT(rlistnames, 0, mkChar("MuHat"));
  SET_STRING_ELT(rlistnames, 1, mkChar("SigmaHat"));
  SET_STRING_ELT(rlistnames, 2, mkChar("CstAdj"));
  SET_STRING_ELT(rlistnames, 3, mkChar("LogLik"));
  setAttrib(rlist, R_NamesSymbol, rlistnames);

  UNPROTECT(nprotect+2);
  return rlist;
}

double LogProd(double x, int lenclust, double *expect, double *betal, double LSEbetaL)
{
  int lc;
  double res = 0;
  double temp = 0;
  double temp2 = 0;
  double etemp2 = 0;
  double temp3 = 0;
  temp = LSEbetaL + x;
  for (lc=0; lc<lenclust; lc++){
    temp2 = betal[lc] + x;
    etemp2 = exp(temp2) + expect[lc];
    temp3 = log(etemp2);
    temp3 = (temp3<DBL_MAX) ? temp3 : DBL_MAX;
    res += temp3;
  }
  res += -exp(temp);
  return res;
}

double DLogProd(double x, int lenclust, double *expect, double *betal, double LSEbetaL)
{
  int lc;
  double res = 0;
  double temp = 0;
  double temp1 = 0;
  double temp2 = 0;
  double exptemp = 0;
  double exptemp2 = 0;
  double temp3 = 0;
  temp = x + LSEbetaL;
  for (lc=0; lc<lenclust; lc++){
    temp1 = x + betal[lc];
    exptemp2 = exp(temp1) + expect[lc];
    temp3 = log(exptemp2);
    temp3 = (temp3<DBL_MAX) ? temp3 : DBL_MAX;
    temp2 = temp1 - temp3;
    exptemp = exp(temp2);
    res += exptemp;
  }
  res = (res<DBL_MAX) ? res : DBL_MAX;
  res += -exp(temp);
  return res;
}

double DMLI(double x, int lenclust, double *expect, double *betal, double LSEbetaL, double var)
{
  double res = 0;
  res = x/var - DLogProd(x,lenclust,expect,betal,LSEbetaL);
  return res;
}

/* Root finding algorithm */
double ZeroDMLI(double xmin, double xmax, double tol, int lenclust, double *expect, double *betal, double LSEbetaL, double var)
{
  double x0 = xmin, x1 = xmax;
  double f0, f1, fc, xc;

  f0 = DMLI(x0,lenclust,expect,betal,LSEbetaL,var);
  f1 = DMLI(x1,lenclust,expect,betal,LSEbetaL,var);
  if(f0==0) return x0;
  if(f1==0) return x1;
  if(f0*f1>0){
    return DBL_MAX;
  }
  for(;;) {
    xc = 0.5*(x0+x1);
    if(fabs(x0-x1)<tol) return xc;
    fc = DMLI(xc,lenclust,expect,betal,LSEbetaL,var);
    if(fc==0) {
      return xc;
    }
    if(f0*fc>0) {
      x0 = xc; f0 = fc;
    }
    else {
      x1 = xc; f1 = fc;
    }
  }
}

double DDLogProd(double x, int lenclust, double *expect, double *betal, double LSEbetaL)
{
  int lc;
  double res = 0;
  double temp = 0;
  double temp1 = 0;
  double temp2 = 0;
  double exptemp = 0;
  double exptemp2 = 0;
  temp = x + LSEbetaL;
  for (lc=0; lc<lenclust; lc++){
    temp1 = x + betal[lc];
    exptemp2 = exp(temp1) + expect[lc];
    if (isinf(exptemp2)) exptemp = 0;
    else {
      temp2 = temp1 + log(expect[lc]) - 2*log(exptemp2);
      exptemp = exp(temp2);
    }
    res += exptemp;
  }
  res = (res<DBL_MAX) ? res : DBL_MAX;
  res += -exp(temp);
  return res;
}

double DDMLI(double x, int lenclust, double *expect, double *betal, double LSEbetaL, double var)
{
  double res = 0;
  res = 1/var - DDLogProd(x,lenclust,expect,betal,LSEbetaL);
  return res;
}

double LLGHQClust(int Npoint, double *logresll, double A, int clust)
{
  double res = 0;
  double temp = 0;
  double resll = 0;
  int k = 0;

  for(k=0; k<Npoint; k++) {
    temp = logresll[k] + A*clust;
    resll += exp(temp);
  }
  res += -log(resll) + A*clust;
  return res;
}
