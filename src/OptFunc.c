#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include "OptFunc.h"

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
  res = -log(resll) + A*clust;
  return res;
}

double LogProd0(double x, double LSEbetaL)
{
  double res = 0;
  res = -exp(x + LSEbetaL);
  return res;
}

double DDMLI0(double x, double LSEbetaL, double var)
{
  double res = 0;
  res = 1/var - LogProd0(x,LSEbetaL);
  return res;
}
