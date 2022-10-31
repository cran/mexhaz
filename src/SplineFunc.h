#ifndef _SplineFunc_h
#define _SplineFunc_h

double Spline1(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT);
double Spline2(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT);
double Spline3(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT);
double NSpl(double x, double *TotKT, double *MatKT, double *NsAdj1, double *NsAdj2, double *BasisB, double *TempD, double *Param, int leB, int leN, int Pos);

double IntSpline23(double (*Spl)(double, double *, double *, double *, double *), double a, double b, double *TotKT, double *MatKT, double *TempD, double *ParamT, double *N, double *lW, int lleg);
double IntNSpl(double a, double b, double *TotKT, double *MatKT, double *NsAdj1, double *NsAdj2, double *BasisB, double *TempD, double *Param, double *N, double *lW, int lleg, int leB, int leN, int Idx);

double DeltaSpline1(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT, int Lsdk, int Idx, double *Res);
double DeltaSpline2(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT, int Lsdk, int Idx, double *Res);
double DeltaSpline3(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT, int Lsdk, int Idx, double *Res);
double DeltaNSpl(double x, double *TotKT, double *MatKT, double *NsAdj1, double *NsAdj2, double *BasisB, double *TempD, double *ParamT, int leB, int leN, int Idx, double *Res);

double IntDSpline23(double (*DSpl)(double, double *, double *, double *, double *, int, int, double *), double a, double b, double *TotKT, double *MatKT, double *TempD, double *ParamT, double *N, double *lW, int lleg, int Lsdk, int Idx, double *TempV, double *Res);
double IntDNSpl(double a, double b, double *TotKT, double *MatKT, double *NsAdj1, double *NsAdj2, double *BasisB, double *TempD, double *ParamT, double *N, double *lW, int lleg, int leB, int leN, int Idx, double *TempV, double *Res);

double IntDSpline23H(double (*DSpl)(double, double *, double *, double *, double *, int, int, double *), double a, double b, double *TotKT, double *MatKT, double *TempD, double *ParamT, double *N, double *lW, int lleg, int Lsdk, int Idx, double *TempV, double *Hess, double *Res);
double IntDNSplH(double a, double b, double *TotKT, double *MatKT, double *NsAdj1, double *NsAdj2, double *BasisB, double *TempD, double *ParamT, double *N, double *lW, int lleg, int leB, int leN, int Idx, double *TempV, double *Hess, double *Res);

double **dmatrix(double *array, int nrow, int ncol);

#endif
