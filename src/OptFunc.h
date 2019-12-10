#ifndef _OptFunc_h
#define _OptFunc_h

double LogProd(double x, int lenclust, double *expect, double *betal, double LSEbetaL);
double DLogProd(double x, int lenclust, double *expect, double *betal, double LSEbetaL);
double DMLI(double x, int lenclust, double *expect, double *betal, double LSEbetaL, double var);
double ZeroDMLI(double xinf, double xmax, double tol, int lenclust, double *expect, double *betal, double LSEbetaL, double var);
double DDLogProd(double x, int lenclust, double *expect, double *betal, double LSEbetaL);
double DDMLI(double x, int lenclust, double *expect, double *betal, double LSEbetaL, double var);
double LLGHQClust(int Npoint, double *logresll, double A, int clust);
double LogProd0(double x, double LSEbetaL);
double DDMLI0(double x, double LSEbetaL, double var);

#endif
