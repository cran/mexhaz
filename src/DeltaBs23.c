/*********************************************/
/* Routine for estimating the variance of    */
/* the log-hazard and log-cumulative hazard  */
/* by the Delta Method                       */
/* (log-hazard described by                  */
/* a quadratic or cubic B-spline)            */
/* Author: H. Charvat                        */
/* Last modified: 2016/05/10                 */
/* Part of the mexhaz package                */
/*********************************************/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>

double DeltaSpline2(double x, int ltotk, double *TotK, double *MatK, double *Diff, double *Param, int Pos, double *Res, int Lsdk);
double DeltaSpline3(double x, int ltotk, double *TotK, double *MatK, double *Diff, double *Param, int Pos, double *Res, int Lsdk);
double IntDSpline23(double (*Spl)(),int pos, double b, int ltotk, double *TotK, double *MatK, double *Diff, double *MyParam, int Pos, double *N, double *lW, int lleg, double *TempV, double *Res, int Lsdk);

SEXP DeltaBs23(SEXP x, SEXP nph, SEXP timecat, SEXP param, SEXP others, SEXP varcov, SEXP k, SEXP n, SEXP lw, SEXP matk, SEXP totk)
{
  SEXP result;
  int lx = length(x);
  int lnph = length(nph);
  int loth = length(others);
  int npar = length(param);
  int lleg = length(n);
  int ltotk = length(totk);
  PROTECT(x = coerceVector(x,REALSXP));
  PROTECT(nph = coerceVector(nph,REALSXP));
  PROTECT(timecat = coerceVector(timecat,INTSXP));
  PROTECT(param = coerceVector(param,REALSXP));
  PROTECT(others = coerceVector(others,REALSXP));
  PROTECT(varcov = coerceVector(varcov,REALSXP));
  PROTECT(k = coerceVector(k,INTSXP));
  PROTECT(n = coerceVector(n,REALSXP));
  PROTECT(lw = coerceVector(lw,REALSXP));
  PROTECT(matk = coerceVector(matk,REALSXP));
  PROTECT(totk = coerceVector(totk,REALSXP));
  PROTECT(result = allocVector(REALSXP,2*lx));
  double *X = REAL(x);
  double *Nph = REAL(nph);
  int *TimeCat = INTEGER(timecat);
  double *Param = REAL(param);
  double *Others = REAL(others);
  double *Varcov = REAL(varcov);
  int *K = INTEGER(k);
  double *N = REAL(n);
  double *lW = REAL(lw);
  double *MatK = REAL(matk);
  double *TotK = REAL(totk);
  double *Result = REAL(result);
  int nnph = lnph/lx;
  int nind = loth/lx;
  int lsdk = ltotk-K[0];
  double tempL, InvtempL;
  double *MyGradL = malloc(npar*sizeof(double));
  double *MyGradS = malloc(npar*sizeof(double));
  double *MyParam = malloc((lsdk+1)*sizeof(double));
  double *Diff = malloc(ltotk*sizeof(double));
  double *Res = malloc((lsdk+1)*sizeof(double));
  double *tempLvec = malloc((lsdk+1)*sizeof(double));
  int i, ii, j, z, t2, t3;
  int Pos = -K[0];
  int Pos2 = K[0]-1;
  double (*Fpt)(double, int, double*, double*, double*, double*, int, double*, int);
  if (K[0]==2){
    Fpt = &DeltaSpline2;
  }
  else {
    Fpt = &DeltaSpline3;
  }

  for (z=0; z<lx; z++){

    t3 = nind*z;
    for (i=0; i<nind; i++){
      MyGradL[i] = Others[i+t3];
      MyGradS[i] = Others[i+t3];
    }

    Result[z] = 0;
    Result[lx+z] = 0;
    MyParam[0] = 0;
    tempLvec[0] = 0;
    tempL = 0;
    t2 = z*nnph;

    for (i=0; i<lsdk; i++){
      MyParam[i+1] = Param[i+nind];
      tempLvec[i+1] = 0;
      ii = 1;
      while (ii<nnph){
	MyParam[i+1] += Param[ii*lsdk+i+nind]*Nph[ii+t2];
	ii++;
      }
    }

    for (j=0; j<TimeCat[z]; j++){
      tempL += IntDSpline23((*Fpt), (Pos2+j), TotK[Pos2+j+1], ltotk, TotK, MatK, Diff, MyParam, Pos, N, lW, lleg, tempLvec, Res, lsdk);
    }
    tempL += IntDSpline23((*Fpt), (Pos2+TimeCat[z]), X[z], ltotk, TotK, MatK, Diff, MyParam, Pos, N, lW, lleg, tempLvec, Res, lsdk);
    InvtempL = 1/tempL;
    (*Fpt)(X[z], ltotk, TotK, MatK, Diff, MyParam, Pos, Res, lsdk);

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

  UNPROTECT(12);
  free(MyGradL);
  free(MyGradS);
  free(MyParam);
  free(Diff);
  free(Res);
  free(tempLvec);
  return result;
}


double DeltaSpline2(double x, int ltotk, double *TotK, double *MatK, double *Diff, double *Param, int Pos, double *Res, int Lsdk){
  int i, j;
  double A, B, res;
  double *MatKT, *TempD, *ParamT, *ResT;
  for (j=0; j<(Lsdk+1); j++){
    Res[j] = 0;
  }
  for (i=0; i<ltotk; i++){
    Pos += ((Diff[i] = x - TotK[i]) > 0);
  }
  MatKT = &MatK[2*Pos];
  TempD = &Diff[Pos];
  ParamT = &Param[Pos];
  ResT = &Res[Pos];
  A = MatKT[0]*TempD[1];
  B = MatKT[1]*TempD[2];
  ResT[0] = B*TempD[2];
  ResT[1] = -(B*TempD[0]+A*TempD[3]);
  ResT[2] = A*TempD[1];
  res = ParamT[2]*ResT[2]+ParamT[1]*ResT[1]+ParamT[0]*ResT[0];
  return res;
}

double DeltaSpline3(double x, int ltotk, double *TotK, double *MatK, double *Diff, double *Param, int Pos, double *Res, int Lsdk){
  int i, j;
  double A, B, C, res;
  double *MatKT, *TempD, *ParamT, *ResT;
  for (j=0; j<(Lsdk+1); j++){
    Res[j] = 0;
  }
  for (i=0; i<ltotk; i++){
    Pos += ((Diff[i] = x - TotK[i]) > 0);
  }
  MatKT = &MatK[4*Pos];
  TempD = &Diff[Pos];
  ParamT = &Param[Pos];
  ResT = &Res[Pos];
  A = MatKT[0]*TempD[2]*TempD[2];
  B = MatKT[1]*TempD[1]*TempD[3]+MatKT[2]*TempD[2]*TempD[4];
  C = MatKT[3]*TempD[3]*TempD[3];
  ResT[0] = -C*TempD[3];
  ResT[1] = (C*TempD[0]+B*TempD[4]);
  ResT[2] = -(B*TempD[1]+A*TempD[5]);
  ResT[3] = A*TempD[2];
  res = ParamT[3]*ResT[3]+ParamT[2]*ResT[2]+ParamT[1]*ResT[1]+ParamT[0]*ResT[0];
  return res;
}

double IntDSpline23(double (*Spl)(),int pos, double b, int ltotk, double *TotK, double *MatK, double *Diff, double *MyParam, int Pos, double *N, double *lW, int lleg, double *TempV, double *Res, int Lsdk){
  double A = 0.5*(b-TotK[pos]);
  double B = 0.5*(b+TotK[pos]);
  int i, j;
  double Result = 0;
  double Temp, Temp2;
  for (i=0; i<lleg; i++){
    Temp = Spl((A*N[i]+B), ltotk, TotK, MatK, Diff, MyParam, Pos, Res, Lsdk);
    Temp2 = exp(lW[i]+Temp);
    Result += Temp2;
    for (j=0; j<(Lsdk+1); j++){
      TempV[j] += A*Res[j]*Temp2;
    }
  }
  Result *= A;
  return Result;
}
