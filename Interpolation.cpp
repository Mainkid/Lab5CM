#define _USE_MATH_DEFINES
#include "Matrix.h"
#include "Interpolation.h"
#include <iostream>
#include<iomanip>
#include<math.h>
#include <cmath>


void Interpolation::MakeDifTable (double a, double b, double step, int variant)
{
	int n = (a - b) / step;
	double* xArray = new double[n];

	for (int i = 0; i < n; i++)
	{
		xArray[i] += step * i;
	}

	for (int i=0; i<n;i++)
		for (int j = 0; j < n + 1; j++)
		{
			if (j == 0)
				std::cout << xArray[i] << " ";
			else if (j==1)
			{
				
				if (variant == 20)
					std::cout << F20(xArray[i]) << " ";
				else if (variant == 24)
					std::cout << F24(xArray[i]) << " ";

			}
			else
			{
				//std::cout<<


			}

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

double Interpolation::RecF(int start, int finish, double* Arr,int variant)
{
	if (start - finish == 1)
	{
		if (variant == 20)
			return (F20(Arr[finish]) - F20(Arr[start])) / (Arr[finish] - Arr[start]);
		else if (variant==24)
			return (F24(Arr[finish]) - F24(Arr[start])) / (Arr[finish] - Arr[start]);
	}
	else
	{



	}


}