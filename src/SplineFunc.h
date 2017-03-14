#ifndef _SplineFunc_h
#define _SplineFunc_h

double Spline2(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT, int Cst);
double Spline3(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT, int Cst);
double NSpl(double x, double *TotKT, double *MatKT, double *NsAdj1, double *NsAdj2, double *BasisB, double *TempD, double *Param, int leB, int leN, int Pos);

double IntSpline23(double (*Spl)(), double a, double b, double *TotKT, double *MatKT, double *TempD, double *ParamT, double *N, double *lW, int lleg, int Cst);
double IntNSpl(double a, double b, double *TotKT, double *MatKT, double *NsAdj1, double *NsAdj2, double *BasisB, double *TempD, double *Param, double *N, double *lW, int lleg, int leB, int leN, int Idx);

double DeltaSpline2(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT, int Cst, int Lsdk, int Idx, double *Res);
double DeltaSpline3(double x, double *TotKT, double *MatKT, double *TempD, double *ParamT, int Cst, int Lsdk, int Idx, double *Res);
double DeltaNSpl(double x, double *TotKT, double *MatKT, double *NsAdj1, double *NsAdj2, double *BasisB, double *TempD, double *ParamT, int leB, int leN, int Idx, double *Res);

double IntDSpline23(double (*DSpl)(), double a, double b, double *TotKT, double *MatKT, double *TempD, double *ParamT, double *N, double *lW, int lleg, int Cst, int Lsdk, int Idx, double *TempV, double *Res);
double IntDNSpl(double a, double b, double *TotKT, double *MatKT, double *NsAdj1, double *NsAdj2, double *BasisB, double *TempD, double *ParamT, double *N, double *lW, int lleg, int leB, int leN, int Idx, double *TempV, double *Res);

double **dmatrix(double *array, int nrow, int ncol);

#endif
