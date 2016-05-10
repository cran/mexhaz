/************************************************/
/* Routine for estimating the cumulative hazard */
/* (when the log-hazard is described            */
/* by a quadratic or cubic B-spline)            */
/* Author: H. Charvat                           */
/* Last modified: 2016/05/10                    */
/* Part of the mexhaz package                   */
/************************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>

double Spline2(double x, int ltotk, double *TotK, double *MatK, double *Diff, double *Param, int Pos);
double Spline3(double x, int ltotk, double *TotK, double *MatK, double *Diff, double *Param, int Pos);
double IntSpline23(double (*Spl)(),int pos, double b, int ltotk, double *TotK, double *MatK, double *Diff, double *MyParam, int Pos, double *N, double *lW, int lleg);

SEXP IntBs23(SEXP x, SEXP nph, SEXP timecat, SEXP param, SEXP k, SEXP n, SEXP lw, SEXP matk, SEXP totk)
{
  SEXP result;
  int lx = length(x);
  int lnph = length(nph);
  int lleg = length(n);
  int ltotk = length(totk);
  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(timecat = coerceVector(timecat,INTSXP));
  PROTECT(param = coerceVector(param,REALSXP));
  PROTECT(k = coerceVector(k,INTSXP));
  PROTECT(n = coerceVector(n,REALSXP));
  PROTECT(lw = coerceVector(lw,REALSXP));
  PROTECT(matk = coerceVector(matk,REALSXP));
  PROTECT(totk = coerceVector(totk,REALSXP));
  PROTECT(result = allocVector(REALSXP,(lx*2+1)));
  double *X = REAL(x);
  double *Nph = REAL(nph);
  int *TimeCat = INTEGER(timecat);
  double *Param = REAL(param);
  int *K = INTEGER(k);
  double *N = REAL(n);
  double *lW = REAL(lw);
  double *MatK = REAL(matk);
  double *TotK = REAL(totk);
  double *Result = REAL(result);
  int nnph = lnph/lx;
  double tempL;
  double *MyParam = malloc((ltotk-K[0]+1)*sizeof(double));
  double *Diff = malloc(ltotk*sizeof(double));
  int i, z, ii, t2;
  int lsdk = ltotk-K[0];
  int Pos = -K[0];
  double (*Fpt)(double, int, double*, double*, double*, double*, int);
  if (K[0]==2){
    Fpt = &Spline2;
  }
  else {
    Fpt = &Spline3;
  }

  int Pos2 = K[0]-1;
  int pos1 = 2*lx;
  Result[pos1] = 0;

  for (z=0; z<lx; z++){
    MyParam[0] = 0;
    for (i=0; i<lsdk; i++){
      MyParam[i+1] = Param[i];
      ii = 1;
      t2 = z*nnph;
      while (ii<nnph){
	MyParam[i+1] += Param[ii*lsdk+i]*Nph[ii+t2];
	ii++;
      }
    }
    tempL = 0;
    for (i=0; i<TimeCat[z]; i++){
      tempL += IntSpline23((*Fpt), (Pos2+i), TotK[Pos2+i+1], ltotk, TotK, MatK, Diff, MyParam, Pos, N, lW, lleg);
    }
    tempL += IntSpline23((*Fpt), (Pos2+TimeCat[z]), X[z], ltotk, TotK, MatK, Diff, MyParam, Pos, N, lW, lleg);
    Result[z] = (*Fpt)(X[z], ltotk, TotK, MatK, Diff, MyParam, Pos);
    Result[z+lx] = log(tempL);
    Result[pos1] += Result[z] + Result[z+lx];
  }

  UNPROTECT(10);
  free(MyParam);
  free(Diff);
  return result;
}


double Spline2(double x, int ltotk, double *TotK, double *MatK, double *Diff, double *Param, int Pos){
  int i;
  double res, A, B;
  double *MatKT, *TempD, *ParamT;
  for (i=0; i<ltotk; i++){
    Pos += ((Diff[i] = x - TotK[i]) > 0);
  }
  MatKT = &MatK[2*Pos];
  TempD = &Diff[Pos];
  ParamT = &Param[Pos];
  A = MatKT[0]*TempD[1];
  B = MatKT[1]*TempD[2];
  res = ParamT[2]*A*TempD[1]-ParamT[1]*(B*TempD[0]+A*TempD[3])+ParamT[0]*B*TempD[2];
  return res;
}


double Spline3(double x, int ltotk, double *TotK, double *MatK, double *Diff, double *Param, int Pos){
  int i;
  double res, A, B, C;
  double *MatKT, *TempD, *ParamT;
  for (i=0; i<ltotk; i++){
    Pos += ((Diff[i] = x - TotK[i]) > 0);
  }
  MatKT = &MatK[4*Pos];
  TempD = &Diff[Pos];
  ParamT = &Param[Pos];
  A = MatKT[0]*TempD[2]*TempD[2];
  B = MatKT[1]*TempD[1]*TempD[3]+MatKT[2]*TempD[2]*TempD[4];
  C = MatKT[3]*TempD[3]*TempD[3];
  res = ParamT[3]*A*TempD[2]-ParamT[2]*(B*TempD[1]+A*TempD[5])+ParamT[1]*(C*TempD[0]+B*TempD[4])-ParamT[0]*C*TempD[3];
  return res;
}

double IntSpline23(double (*Spl)(),int pos, double b, int ltotk, double *TotK, double *MatK, double *Diff, double *MyParam, int Pos, double *N, double *lW, int lleg){
  double A = 0.5*(b-TotK[pos]);
  double B = 0.5*(b+TotK[pos]);
  int i;
  double Result = 0;
  double Temp;
  for (i=0; i<lleg; i++){
    Temp = Spl((A*N[i]+B), ltotk, TotK, MatK, Diff, MyParam, Pos);
    Result += exp(lW[i]+Temp);
  }
  Result *= A;
  return Result;
}
