#define _USE_MATH_DEFINES
#include "Matrix.h"
#include "Interpolation.h"
#include <iostream>
#include<iomanip>
#include<math.h>
#include <cmath>


void Interpolation::MakeDifTable (double a, double b, double step, int variant)
{

	std::cout << "Интерполяционная формула Ньютона" << std::endl;
	std::cout << "Таблица разделенных разностей" << std::endl;
	std::cout<<std::endl;

	int n = (b-a) / step;
	double* xArray = new double[n+1];

	for (int i = 0; i < n+1; i++)
	{
		xArray[i] = a+step * (i);
	}

	for (int i=0; i<n+1;i++)
		for (int j = 0; j < n + 2-i; j++)
		{
			if (j == 0)
				std::cout << xArray[i] << " ";
			else if (j==1)
			{
				
				if (variant == 20)
					std::cout << F20(xArray[i]) << " ";
				else if (variant == 24)
					std::cout << F24(xArray[i]) << " ";
				else if (variant ==19)
					std::cout << F19(xArray[i]) << " ";
			}
			else
			{
				std::cout << RecF(i, j - 1+i, xArray, variant) << " ";


			}

			if (j == n+1-i)
				std::cout << std::endl;
		}



}

double Interpolation::F20(double x)
{
	return pow(M_E, x) - pow(x, 2) + 1;
}

double Interpolation::F24(double x)
{
	return 2 * pow(M_E, x) - 2 * x + 3;
}

double Interpolation::F19(double x)
{
	double new1 = pow(3, x) + 2 * x - 5;
	return pow(3, x) + 2 * x - 5;
}

double Interpolation::RecF(int start, int finish, double* Arr,int variant)
{
	if (finish -start == 1)
	{
		if (variant == 20)
			return (F20(Arr[finish]) - F20(Arr[start])) / (Arr[finish] - Arr[start]);
		else if (variant==24)
			return (F24(Arr[finish]) - F24(Arr[start])) / (Arr[finish] - Arr[start]);
		else if (variant == 19)
		{
			double new1 = (F19(Arr[finish]) - F19(Arr[start])) / (Arr[finish] - Arr[start]);
			return (F19(Arr[finish]) - F19(Arr[start])) / (Arr[finish] - Arr[start]);
		}
	}
	else
	{
		if (variant == 20)
			return (F20(RecF(start + 1, finish, Arr, variant)) - F20(RecF(start, finish - 1, Arr, variant))) / (Arr[finish] - Arr[start]);
		else if (variant==24)
			return (F24(RecF(start + 1, finish, Arr, variant)) - F24(RecF(start, finish - 1, Arr, variant))) / (Arr[finish] - Arr[start]);
		else if (variant==19)
			return (RecF(start + 1, finish, Arr, variant) - RecF(start, finish - 1, Arr, variant)) / (Arr[finish] - Arr[start]);
	}


}