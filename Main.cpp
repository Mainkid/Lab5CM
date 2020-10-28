#include<iostream>
#include <iomanip>
#include<Windows.h>
#include <fstream>
#include <ctime>
#include "Matrix.h"
#include "RootSearchMethods.h"
#include "SquareApproximation.h"

using namespace std;

int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	SquareApproximation sqAppr(variants::var20);
	sqAppr.DiscreteMethod();
	printf("\n");
	sqAppr.IntegrationMethod();

	return 0;
}