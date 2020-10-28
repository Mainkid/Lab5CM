#include<iostream>
#include <iomanip>
#include<Windows.h>
#include <fstream>
#include <ctime>
#include "Matrix.h"
#include "RootSearchMethods.h"
#include "Interpolation.h"

using namespace std;

int main()
{
	double a = 1;
	double b = 2;
	double h = 0.2;
	int n = 5;
	double* xArray = new double[n + 1];
	for (int i = 0; i < n + 1; i++)
		xArray[i] = a + i * h;
	

	double** difTable = Interpolation::InitDifTable(n);//del
	double* difVector = new double[n + 1];
	

	Interpolation::MakeDifTable(difTable,xArray,difVector,19);
	Interpolation::NewtonInterpolation(difVector, n + 1, 19, a, b, n, xArray, n);
	Interpolation::CubeSpline(xArray, difVector, n + 1, 19);

	difTable = Interpolation::InitDifTable(n);
	difVector = new double[n + 1];
	Interpolation::ReversedInterpolation(xArray, difVector, difTable, n, 19);




	return 0;
}