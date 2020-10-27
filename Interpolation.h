#pragma once

class Interpolation {

public:
	static double F19(double x);
	static double F20(double x);
	static double F24(double x);
	static void MakeDifTable(double a, double b, double step, int variant);
	static double RecF(int start, int finish, double* Arr,int variant);
	static void PrintHead();
	static void NewtonInterpolation(double* difVector, int difVectorN, int variant, double start, double finish, int n,double* Arr,int nArr);
	static double CountPn(double x, int n, double* Arr, int xArr,double* difVector);
	static double FindPolynom(double x,int kol, double* Arr);
	static double dF20(double x, int step);
	static double dF24(double x, int step);
	static double dF19(double x, int step);
	static double d1F19(double x);
	static double d1F20(double x);
	static double d1F24(double x);
	static void PrintBodyNewton(double x, double fx, double pn, double delta, double err);
	static void CubeSpline(double* xArray, double* difVector, int N, int variant);
	static void GetM(double* xArray, double* difVector, int N, int variant,double M5);
	static void PrintSplineHead();
	static void PrintBodySpline(double x, double df, double m, double delta, double err);


};