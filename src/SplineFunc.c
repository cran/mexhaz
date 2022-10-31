#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include "SplineFunc.h"

double Spline1(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT){
  int i;
  double res, A;
  for (i=0; i<2; i++){
    TempD[i] = x - TotKT[i];
  }
  A = MatKT[0];
  res = ParamT[1]*A*TempD[0]-ParamT[0]*A*TempD[1];
  return res;
}

double Spline2(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT){
  int i;
  double res, A, B;
  for (i=0; i<4; i++){
    TempD[i] = x - TotKT[i];
  }
  A = MatKT[0]*TempD[1];
  B = MatKT[1]*TempD[2];
  res = ParamT[2]*A*TempD[1]-ParamT[1]*(B*TempD[0]+A*TempD[3])+ParamT[0]*B*TempD[2];
  return res;
}

double Spline3(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT){
  int i;
  double res, A, B, C;
  for (i=0; i<6; i++){
    TempD[i] = x - TotKT[i];
  }
  A = MatKT[0]*TempD[2]*TempD[2];
  B = MatKT[1]*TempD[1]*TempD[3]+MatKT[2]*TempD[2]*TempD[4];
  C = MatKT[3]*TempD[3]*TempD[3];
  res = ParamT[3]*A*TempD[2]-ParamT[2]*(B*TempD[1]+A*TempD[5])+ParamT[1]*(C*TempD[0]+B*TempD[4])-ParamT[0]*C*TempD[3];
  return res;
}

double NSpl(double x, double *TotKT, double *MatKT, double *NsAdj1, double *NsAdj2, double *BasisB, double *TempD, double *Param, int leB, int leN, int Idx){
  int i, j, t1;
  double A, B, C;
  double res = 0;
  double TempN;
  double *BasisT;
  for (i=0; i<leB; i++){
    BasisB[i] = 0;
  }
  if ((TempD[0] = x - NsAdj2[0])<=0){
    BasisB[1] = 1 + TempD[0]*NsAdj2[1];
    BasisB[2] = TempD[0]*NsAdj2[2];
  }
  else if ((TempD[0] = x - NsAdj2[3])>0){
    BasisB[leB-2] = TempD[0]*NsAdj2[4];
    BasisB[leB-1] = 1 + TempD[0]*NsAdj2[5];
  }
  else {
    for (i=0; i<6; i++){
      TempD[i] = x - TotKT[i];
    }
    BasisT = &BasisB[Idx];
    A = MatKT[0]*TempD[2]*TempD[2];
    B = MatKT[1]*TempD[1]*TempD[3]+MatKT[2]*TempD[2]*TempD[4];
    C = MatKT[3]*TempD[3]*TempD[3];
    BasisT[0] = -C*TempD[3];
    BasisT[1] = (C*TempD[0]+B*TempD[4]);
    BasisT[2] = -(B*TempD[1]+A*TempD[5]);
    BasisT[3] = A*TempD[2];
  }
  for (i=0; i<leN; i++){
    TempN = 0;
    t1 = i*(leB-2);
    for (j=2; j<leB; j++){
      TempN += BasisB[j]*NsAdj1[j-2+t1];
    }
    res += Param[i]*TempN;
  }
  return res;
}

double IntSpline23(double (*Spl)(double, double *, double *, double *, double *), double a, double b, double *TotKT, double *MatKT, double *TempD, double *ParamT, double *N, double *lW, int lleg){
  int i;
  double A = 0.5*(b-a);
  double B = 0.5*(b+a);
  double Result = 0;
  double Temp;
  for (i=0; i<lleg; i++){
    Temp = Spl((A*N[i]+B), TotKT, MatKT, TempD, ParamT);
    Result += exp(lW[i]+Temp);
  }
  Result *= A;
  return Result;
}

double IntNSpl(double a, double b, double *TotKT, double *MatKT, double *NsAdj1, double *NsAdj2, double *BasisB, double *TempD, double *Param, double *N, double *lW, int lleg, int leB, int leN, int Idx){
  int i;
  double A = 0.5*(b-a);
  double B = 0.5*(b+a);
  double Result = 0;
  double Temp;
  for (i=0; i<lleg; i++){
    Temp = NSpl((A*N[i]+B), TotKT, MatKT, NsAdj1, NsAdj2, BasisB, TempD, Param, leB, leN, Idx);
    Result += exp(lW[i]+Temp);
  }
  Result *= A;
  return Result;
}

double DeltaSpline1(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT, int Lsdk, int Idx, double *Res){
  int i, j;
  double res, A;
  double *ResT;
  for (j=0; j<(Lsdk+1); j++){
    Res[j] = 0;
  }
  for (i=0; i<2; i++){
    TempD[i] = x - TotKT[i];
  }
  ResT = &Res[Idx];
  A = MatKT[0];
  ResT[0] = -A*TempD[1];
  ResT[1] = A*TempD[0];
  res = ParamT[1]*ResT[1]+ParamT[0]*ResT[0];
  return res;
}

double DeltaSpline2(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT, int Lsdk, int Idx, double *Res){
  int i, j;
  double res, A, B;
  double *ResT;
  for (j=0; j<(Lsdk+1); j++){
    Res[j] = 0;
  }
  for (i=0; i<4; i++){
    TempD[i] = x - TotKT[i];
  }
  ResT = &Res[Idx];
  A = MatKT[0]*TempD[1];
  B = MatKT[1]*TempD[2];
  ResT[0] = B*TempD[2];
  ResT[1] = -(B*TempD[0]+A*TempD[3]);
  ResT[2] = A*TempD[1];
  res = ParamT[2]*ResT[2]+ParamT[1]*ResT[1]+ParamT[0]*ResT[0];
  return res;
}

double DeltaSpline3(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT, int Lsdk, int Idx, double *Res){
  int i, j;
  double res, A, B, C;
  double *ResT;
  for (j=0; j<(Lsdk+1); j++){
    Res[j] = 0;
  }
  for (i=0; i<6; i++){
    TempD[i] = x - TotKT[i];
  }
  ResT = &Res[Idx];
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

double DeltaNSpl(double x, double *TotKT, double *MatKT, double *NsAdj1, double *NsAdj2, double *BasisB, double *TempD, double *ParamT, int leB, int leN, int Idx, double *Res){
  int i, j, t1;
  double A, B, C;
  double res = 0;
  double *BasisT;
  for (i=0; i<leB; i++){
    BasisB[i] = 0;
  }
  if ((TempD[0] = x - NsAdj2[0])<=0){
    BasisB[1] = 1 + TempD[0]*NsAdj2[1];
    BasisB[2] = TempD[0]*NsAdj2[2];
  }
  else if ((TempD[0] = x - NsAdj2[3])>0){
    BasisB[leB-2] = TempD[0]*NsAdj2[4];
    BasisB[leB-1] = 1 + TempD[0]*NsAdj2[5];
  }
  else {
    for (i=0; i<6; i++){
      TempD[i] = x - TotKT[i];
    }
    BasisT = &BasisB[Idx];
    A = MatKT[0]*TempD[2]*TempD[2];
    B = MatKT[1]*TempD[1]*TempD[3]+MatKT[2]*TempD[2]*TempD[4];
    C = MatKT[3]*TempD[3]*TempD[3];
    BasisT[0] = -C*TempD[3];
    BasisT[1] = (C*TempD[0]+B*TempD[4]);
    BasisT[2] = -(B*TempD[1]+A*TempD[5]);
    BasisT[3] = A*TempD[2];
  }
  for (i=0; i<leN; i++){
    Res[i] = 0;
    t1 = i*(leB-2);
    for (j=2; j<leB; j++){
      Res[i] += BasisB[j]*NsAdj1[j-2+t1];
    }
    res += ParamT[i]*Res[i];
  }
  return res;
}

double IntDSpline23(double (*DSpl)(double, double *, double *, double *, double *, int, int, double *), double a, double b, double *TotKT, double *MatKT, double *TempD, double *ParamT, double *N, double *lW, int lleg, int Lsdk, int Idx, double *TempV, double *Res){
  double A = 0.5*(b-a);
  double B = 0.5*(b+a);
  int i, j;
  double Result = 0;
  double Temp, Temp2;
  for (i=0; i<lleg; i++){
    Temp = DSpl((A*N[i]+B), TotKT, MatKT, TempD, ParamT, Lsdk, Idx, Res);
    Temp2 = exp(lW[i]+Temp);
    Result += Temp2;
    for (j=0; j<(Lsdk+1); j++){
      TempV[j] += A*Res[j]*Temp2;
    }
  }
  Result *= A;
  return Result;
}

double IntDNSpl(double a, double b, double *TotKT, double *MatKT, double *NsAdj1, double *NsAdj2, double *BasisB, double *TempD, double *ParamT, double *N, double *lW, int lleg, int leB, int leN, int Idx, double *TempV, double *Res){
  double A = 0.5*(b-a);
  double B = 0.5*(b+a);
  int i, j;
  double Result = 0;
  double Temp, Temp2;
  for (i=0; i<lleg; i++){
    Temp = DeltaNSpl((A*N[i]+B), TotKT, MatKT, NsAdj1, NsAdj2, BasisB, TempD, ParamT, leB, leN, Idx, Res);
    Temp2 = exp(lW[i]+Temp);
    Result += Temp2;
    for (j=0; j<leN; j++){
      TempV[j] += A*Res[j]*Temp2;
    }
  }
  Result *= A;
  return Result;
}

double IntDSpline23H(double (*DSpl)(double, double *, double *, double *, double *, int, int, double *), double a, double b, double *TotKT, double *MatKT, double *TempD, double *ParamT, double *N, double *lW, int lleg, int Lsdk, int Idx, double *TempV, double *Hess, double *Res){
  double A = 0.5*(b-a);
  double B = 0.5*(b+a);
  int i, j, k;
  double Result = 0;
  double Temp, Temp2;
  for (i=0; i<lleg; i++){
    Temp = DSpl((A*N[i]+B), TotKT, MatKT, TempD, ParamT, Lsdk, Idx, Res);
    Temp2 = exp(lW[i]+Temp);
    Result += Temp2;
    for (j=0; j<Lsdk; j++){
      TempV[j] += A*Res[j+1]*Temp2;
      for (k=0; k<Lsdk; k++){
	Hess[k+j*Lsdk] += A*Res[j+1]*Res[k+1]*Temp2;
      }
    }
  }
  Result *= A;
  return Result;
}

double IntDNSplH(double a, double b, double *TotKT, double *MatKT, double *NsAdj1, double *NsAdj2, double *BasisB, double *TempD, double *ParamT, double *N, double *lW, int lleg, int leB, int leN, int Idx, double *TempV, double *Hess, double *Res){
  double A = 0.5*(b-a);
  double B = 0.5*(b+a);
  int i, j, k;
  double Result = 0;
  double Temp, Temp2;
  for (i=0; i<lleg; i++){
    Temp = DeltaNSpl((A*N[i]+B), TotKT, MatKT, NsAdj1, NsAdj2, BasisB, TempD, ParamT, leB, leN, Idx, Res);
    Temp2 = exp(lW[i]+Temp);
    Result += Temp2;
    for (j=0; j<leN; j++){
      TempV[j] += A*Res[j]*Temp2;
      for (k=0; k<leN; k++){
	Hess[k+j*leN] += A*Res[j]*Res[k]*Temp2;
      }
    }
  }
  Result *= A;
  return Result;
}

/* Based on code from the survival package */
double **dmatrix(double *array, int nrow, int ncol)
    {
    int i;
    double **pointer;

    pointer = (double **) R_alloc(ncol, sizeof(double *));
    for (i=0; i<ncol; i++) {
	pointer[i] = array;
	array += nrow;
	}
    return(pointer);
}
