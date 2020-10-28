#define _USE_MATH_DEFINES
#include "Matrix.h"
#include "Interpolation.h"
#include <iostream>
#include<iomanip>
#include<math.h>
#include <cmath>

const int width = 15;
const double h = 0.2;
const double a = 1;
const double b = 2;

const double c19=1;
const double c20=0;
const double c24 = 0;



void Interpolation::MakeDifTable (double** difTable, double* xArray,double* difVector, int variant)
{

	std::cout << "Интерполяционная формула Ньютона" << std::endl;
	std::cout << "Таблица разделенных разностей" << std::endl;
	std::cout<<std::endl;

	int n = 5;

	double tmp;

	for (int i=0; i<n+1;i++)
		for (int j = 0; j < n + 2-i; j++)
		{
			if (j == 0)
			{
				std::cout << xArray[i] << " ";
				difTable[i][j] = xArray[i];
			}
			else if (j==1)
			{
				
				if (variant == 20)
				{
					tmp = F20(xArray[i],0);
					std::cout << tmp<< " ";
				}
				else if (variant == 24)
				{
					tmp = F24(xArray[i],0);
					std::cout << tmp << " ";
				}
				else if (variant == 19)
				{
					tmp = F19(xArray[i],0);
					std::cout << tmp << " ";

				}
				difTable[i][j] = tmp;
			}
			else
			{
				tmp = RecF(i, j - 1 + i, xArray, variant, 0);
				std::cout << tmp<< " ";
				difTable[i][j] = tmp;

			}
			if (i == 0&&j>0)
				difVector[j-1] = tmp;


			if (j == n+1-i)
				std::cout << std::endl;
		}
	ShowDifTable(difTable, n + 1);
	
}

void Interpolation::MakeDifTable(double** difTable, double* xArray, double* yArray, double* difVector, int variant)
{
	int n = 5;
	double tmp=0;
	double c = 0;
	if (variant == 19)
		c = c19;
	else if (variant == 20)
		c = c20;
	else if (variant == 24)
		c = c24;

	for (int i = 0; i < n + 1; i++)
		for (int j = 0; j < n + 2 - i; j++)
		{
			if (j == 1)
			{
				
				difTable[i][j] = xArray[i];
			}
			else if (j == 0)
			{

				
				difTable[i][j] = yArray[i];
			}
			else
			{
				tmp = RecF2(i, j - 1 + i, yArray,xArray, variant,c);

				difTable[i][j] = tmp;

			}
			if (i == 0 && j > 0)
				difVector[j - 1] = tmp;

		}
	ShowDifTable(difTable, n + 1);
}

void Interpolation::ShowDifTable(double** difTable, int n)
{

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n+1-i; j++)
		{
			std::cout << " ";
			std::cout << std::setw(12) << std::left << difTable[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

double Interpolation::F20(double x, double c)
{
	if (c == 0)
		return pow(M_E, x) - pow(x, 2) + 1;
	else
		return  pow(M_E, x) - pow(x, 2) + 1 - c;
}

double Interpolation :: dF20(double x,int step)
{
	return pow(M_E, x);
}

double Interpolation::d1F20(double x)
{
	return (pow(M_E, x) - 2 * x);
}


double Interpolation::F24(double x,double c)
{
	if (c == 0)
		return 2 * pow(M_E, x) - 2 * x + 3;
	else
		return 2 * pow(M_E, x) - 2 * x + 3-c;
}

double Interpolation :: dF24(double x, int step)
{
	return 2 * pow(M_E, x);
}

double Interpolation::d1F24(double x)
{

	return 2 * pow(M_E, x) - 2;
}

double Interpolation::F19(double x,double c)
{
if (c==0)
	return pow(3, x) + 2 * x - 5;
else 
	return pow(3, x) + 2 * x - 5-c;
}

double Interpolation::dF19(double x, int step)
{
	return pow(3, x)*pow(log(3), step);
}

double Interpolation::d1F19(double x)
{
	return pow(3, x)*log(3) + 2;
}

double Interpolation::RecF(int start, int finish, double* Arr,int variant,double c)
{
	if (finish -start == 1)
	{
		if (variant == 20)
			return (F20(Arr[finish],c) - F20(Arr[start],c)) / (Arr[finish] - Arr[start]);
		else if (variant==24)
			return (F24(Arr[finish],c) - F24(Arr[start],c)) / (Arr[finish] - Arr[start]);
		else if (variant == 19)
		{
			double new1 = (F19(Arr[finish],c) - F19(Arr[start],c)) / (Arr[finish] - Arr[start]);
			return (F19(Arr[finish],c) - F19(Arr[start],c)) / (Arr[finish] - Arr[start]);
		}
	}
	else
	{//????????????????
		if (variant == 20)
			return (F20(RecF(start + 1, finish, Arr, variant,c),c) - F20(RecF(start, finish - 1, Arr, variant, c),c)) / (Arr[finish] - Arr[start]);
		else if (variant==24)
			return (F24(RecF(start + 1, finish, Arr, variant,c),c) - F24(RecF(start, finish - 1, Arr, variant, c),c)) / (Arr[finish] - Arr[start]);
		else if (variant==19)
			return (RecF(start + 1, finish, Arr, variant, c) - RecF(start, finish - 1, Arr, variant, c)) / (Arr[finish] - Arr[start]);
	}


}

double Interpolation::RecF2(int start, int finish, double* Arr, double* xArr,int variant, double c)
{
	if (finish - start == 1)
	{
			double new1 = (xArr[finish] - xArr[start]) / (Arr[finish] - Arr[start]);
			return ((xArr[finish] - xArr[start]) / (Arr[finish] - Arr[start]));
		
	}
	else
	{//????????????????
		
			return (RecF2(start + 1, finish, Arr,xArr, variant, c) - RecF2(start, finish - 1, Arr,xArr, variant, c)) / (Arr[finish] - Arr[start]);
	}


}

void Interpolation::PrintHead()
{
	std::cout<< std::setw(width) << std::left << "x" << std::setw(width) << std::left << "f(x)" << std::setw(width) << std::left << "Pn(x)" << std::setw(width) << std::left <<
		"Delta" << std::setw(width) << std::left << "Оценка" << std::endl;

}

void Interpolation::PrintSplineHead()
{
std::cout<< std::setw(width) << std::left << "x" << std::setw(width) << std::left << "df/dx(x)" << std::setw(width) << std::left << "m" << std::setw(width) << std::left <<
		"Delta" << std::setw(width) << std::left << "Оценка" << std::endl;

}

void Interpolation :: PrintSplineResultHead()
{
	std::cout << std::setw(width) << std::left << "x" << std::setw(width) << std::left << "f(x)" << std::setw(width) << std::left << "S31(f;x)" << std::setw(width) << std::left <<
		"Delta" << std::setw(width) << std::left << "Оценка" << std::endl;

}

void Interpolation::PrintSplineResultBody(double x, double fx,double S31, double Delta, double Err)
{
	std::cout << std::setw(width) << std::left << x << std::setw(width) << std::left << fx << std::setw(width) << std::left << S31 << std::setw(width) << std::left <<
		Delta << std::setw(width) << std::left << Err << std::endl;

}

void Interpolation::PrintBodySpline(double x, double df, double m, double delta, double err)
{
	std::cout << std::setw(width) << std::left << x << std::setw(width) << std::left <<df << std::setw(width) << std::left << m << std::setw(width) << std::left <<
		delta << std::setw(width) << std::left << err << std::endl;
}

void Interpolation::PrintBodyNewton(double x,double fx, double pn, double delta, double err)
{
	std::cout<<std::setw(width) << std::left << x << std::setw(width) << std::left << fx << std::setw(width) << std::left << pn << std::setw(width) << std::left <<
		delta << std::setw(width) << std::left << err << std::endl;

}

void Interpolation::NewtonInterpolation( double* difVector , int difVectorN, int variant,double start,double finish,int n,double* Arr, int nArr)
{
	double* xArray = new double[5];
	double func;
	double Pn;
	double Rn;
	double M6;

	if (variant == 19)
		M6= dF19(2,6);
	else if (variant==20)
		M6=dF20(2,6);
	else if (variant == 24)
		M6=dF24(2,6);
	std::cout << std::endl;
	std::cout <<"M6: "<<M6<< std::endl;
	std::cout << std::endl;
	for (int i = 0; i < 5; i++)
		xArray[i] = start + (i + 0.5)*0.2;

	PrintHead();

	for (int i = 0; i < 5; i++)
	{
		if (variant == 19)
			func = F19(xArray[i],0);
		else if (variant == 20)
			func = F20(xArray[i],0);
		else if (variant == 24)
			func = F24(xArray[i],0);

		Pn=CountPn(xArray[i], 5, Arr, nArr,difVector);
		Rn = abs(FindPolynom(xArray[i], 6, Arr)*M6/720);
		/*std::cout << func<< std::endl;
		std::cout << Pn - func << std::endl;
		std::cout << Pn << std::endl;*/

		PrintBodyNewton(xArray[i], func, Pn, Pn - func, Rn);
	}

	/*for (int i = 0; i < difVectorN; i++)
		std::cout << difVector[i] << " ";*/



}

double Interpolation::CountPn(double x, int n, double* Arr, int xArr,double* difVector)
{
	double sum=difVector[0];

	for (int i = 1; i < n; i++)
	{
		sum += difVector[i] * FindPolynom(x,i, Arr);
		
	}

	return sum;
}

double Interpolation::FindPolynom(double x,int kol, double* Arr)
	{
		if (kol == 0)
			return 1;
		else
		{
			return (x - Arr[kol - 1])*FindPolynom(x, kol - 1, Arr);
		}

	}

void Interpolation::CubeSpline(double* xArray, double* difVector, int N, int variant)
{
	double M5;
	
	double* mVector = new double[N];
	std::cout << std::endl;
	std::cout << "Интерполяция кубическим сплайном" << std::endl;
	std::cout<<std::endl;
	if (variant == 19)
	{
		M5 = dF19(2, 5);
		mVector[0] = d1F19(xArray[0]);
		mVector[N - 1] = d1F19(xArray[N - 1]);
	}
	else if (variant == 20)
	{
		M5 = dF20(2, 5);
		mVector[0] = d1F20(xArray[0]);
		mVector[N - 1] = d1F20(xArray[N - 1]);
	}
	else if (variant == 24)
	{
		M5 = dF24(2, 5);
		mVector[0] = d1F24(xArray[0]);
		mVector[N - 1] = d1F24(xArray[N - 1]);
	}

	std::cout << std::endl;
	std::cout << "M5: " << M5 << std::endl;
	std::cout << std::endl;

	GetM(xArray, difVector, N, variant,M5);



}

double Interpolation::Phi0(double tau)
{
	return (1 + 2 * tau)*(1 - tau)*(1 - tau);
}

double Interpolation::Phi1(double tau)
{

	return tau * (1 - tau)*(1 - tau);
}

void Interpolation::GetM(double* xArray, double* difVector, int N, int variant,double M5)
{
	double lyambda = 0.5;
	double mu = 0.5;
	double** Matrix = new double*[N - 2];
	double*m = new double[N];
	for (int i = 0; i < N; i++)
		m[i] = 0;
	double f = 0;
	int j = 0;
	for (int i = 0; i < N - 2; i++)
	{
		Matrix[i] = new double[N - 1];

	}

	for (int i = 0; i < N - 2; i++)
		for (int j = 0; j < N - 1; j++)
			Matrix[i][j] = 0;


	for (int i = 1; i < N - 1; i++)
	{
		if (i == 1)
		{
			Matrix[0][0] = 2;
			Matrix[0][1] = 0.5;
			if (variant == 19)
				Matrix[0][N - 2] = 1.5*(F19(xArray[i + 1],0) - F19(xArray[i - 1],0)) / 0.2 - 0.5*d1F19(xArray[i - 1]);
			else if (variant == 20)
				Matrix[0][N - 2] = 1.5*(F20(xArray[i + 1],0) - F20(xArray[i - 1],0)) / 0.2 - 0.5*d1F20(xArray[i - 1]);
			else if (variant == 24)
				Matrix[0][N - 2] = 1.5*(F24(xArray[i + 1],0) - F24(xArray[i - 1],0)) / 0.2 - 0.5*d1F24(xArray[i - 1]);
		}
		else if (i == N - 2)
		{
			Matrix[N - 3][N - 4] = 0.5;
			Matrix[N - 3][N - 3] = 2;

			if (variant == 19)
				Matrix[N - 3][N - 2] = 1.5*(F19(xArray[i + 1],0) - F19(xArray[i - 1],0)) / 0.2 - 0.5*d1F19(xArray[i + 1]);
			else if (variant == 20)
				Matrix[N - 3][N - 2] = 1.5*(F20(xArray[i + 1],0) - F20(xArray[i - 1],0)) / 0.2 - 0.5*d1F20(xArray[i + 1]);
			else if (variant == 24)
				Matrix[N - 3][N - 2] = 1.5*(F24(xArray[i + 1],0) - F24(xArray[i - 1],0)) / 0.2 - 0.5*d1F24(xArray[i + 1]);

		}
		else
		{
			Matrix[i - 1][0 + j] = 0.5;
			Matrix[i - 1][1 + j] = 2;
			Matrix[i - 1][2 + j] = 0.5;

			if (variant == 19)
			{
				double q = xArray[i + 1];
				double m = xArray[i - 1];
				Matrix[i - 1][N - 2] = 1.5*(F19(xArray[i + 1],0) - F19(xArray[i - 1],0)) / 0.2;
				double tt = Matrix[i - 1][N - 2];
				double weq;
			}
			else if (variant == 20)
				Matrix[i - 1][N - 2] = 1.5*(F20(xArray[i + 1],0) - F20(xArray[i - 1],0)) / 0.2;
			else if (variant == 24)
				Matrix[i - 1][N - 2] = 1.5*(F24(xArray[i + 1],0) - F24(xArray[i - 1],0)) / 0.2;

			j++;
		}


	}
	//for (int i = 0; i < N - 2; i++)
	//{
	//	for (int j = 0; j < N - 1; j++)
	//		std::cout << Matrix[i][j] << " ";
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
	int ras = N - 2;


	float  tmp;
	int k;

	for (int i = 0; i < ras; i++)
	{
		tmp = Matrix[i][i];
		for (j = ras; j >= i; j--)
			Matrix[i][j] /= tmp;
		for (j = i + 1; j < ras; j++)
		{
			tmp = Matrix[j][i];
			for (k = ras; k >= i; k--)
			{
				Matrix[j][k] -= tmp * Matrix[i][k];
				if (abs(Matrix[j][k]) < 0.0001)
					Matrix[j][k] = 0;
			}
		}
	}

	//for (int i = 0; i < N - 2; i++)
	//{
	//	for (int j = 0; j < N - 1; j++)
	//		std::cout << Matrix[i][j] << " ";
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;

	if (variant == 19)
	{

		m[0] = d1F19(xArray[0]);
		m[N-1] = d1F19(xArray[N-1]);
	}
	else if (variant == 20)
	{
		m[0] = d1F20(xArray[0]);
		m[N - 1] = d1F20(xArray[N - 1]);

	}
	else if (variant == 24)
	{
		m[0] = d1F24(xArray[0]);
		m[N - 1] = d1F24(xArray[N - 1]);

	}

	/*for (int i = N - 3; i >= 1; i--)
	{
		double d = 0;
		for (int j = i + 1; j <= N - 3; j++)
		{
			d += Matrix[i][j] * Matrix[j][N - 1];
		}

	}*/
	double d;


	for (int i = N-3; i>=0; i--)
	{
		d = 0;
		for (int j = i+1; j < N - 2; j++)
		{
			d += Matrix[i][j]*m[j+1];
		}
		m[i+1] = Matrix[i][N - 2] - d;
		d = m[i+1];
	}

	/*for (int i = 0; i < N; i++)
		std::cout << m[i] << " ";*/
	PrintSplineHead();

	double func;
	for (int i = 0; i < N; i++)
	{
		if (variant == 19)
			func = d1F19(xArray[i]);
		else if (variant == 20)
			func = d1F20(xArray[i]);
		else if (variant == 24)
			func = d1F24(xArray[i]);
		PrintBodySpline(xArray[i], func, m[i], abs(func - m[i]), M5*pow(0.2,4)/60);
	}


	double M4;

	if (variant == 19)
	{
		M4 = dF19(2, 4);
	}
	else if (variant == 20)
	{
		M4 = dF20(2, 4);
	}
	else if (variant == 24)
	{
		M4 = dF24(2, 4);
	}
	std::cout << std::endl;
	std::cout <<"M4: "<< M4<<std::endl;

	double *Array = new double[N - 1];

	for (int i = 0; i < N; i++)
	{

		Array[i] = a + (i + 0.5)*h;
	}

	double tau = 0.5;
	double phi0 = (1 + 2 * tau)*(1 - tau)*(1 - tau);
	double phi1 = tau*(1 - tau)*(1 - tau);
	double S31;
	double fi;
	double fiplus1;

	std::cout << std::endl;

	PrintSplineResultHead();
	double er = (M4 / 384 + M5 * h / 240)*h*h*h*h;
	for (int i = 0; i <N-1; i++)
	{
		
		if (variant == 19)
		{
			fi = F19(xArray[i],0);
			fiplus1 = F19(xArray[i + 1],0);
			f= F19(Array[i],0);
		}
		else if (variant == 20)
		{
			fi = F20(xArray[i],0);
			fiplus1 = F20(xArray[i + 1],0);
			f = F20(Array[i],0);
		}
		else if (variant == 24)
		{
			fi = F24(Array[i],0);
			fiplus1 = F24(xArray[i + 1],0);
			f = F24(Array[i],0);
		}

		S31 = Phi0(tau)*fi + Phi0(1 - tau)*fiplus1 + h * (Phi1(tau)*m[i] - Phi1(1 - tau)*m[i + 1]);

		PrintSplineResultBody(xArray[i], f, S31, abs(S31 - f), er);

	}

}

double** Interpolation::InitDifTable(int n)
{
	double** difTable = new double*[n + 1];
	for (int i = 0; i < n + 1; i++)
	{
		difTable[i] = new double[n + 2];
		
	}

	for (int i = 0; i < n + 1; i++)
		for (int j = 0; j < n + 2; j++)
			difTable[i][j] = 0;
	return difTable;
}

void Interpolation::ReversedInterpolation(double* xArray, double* difVector,double** difTable,int n,int variant)
{
	std::cout << std::endl;

	double* yArray = new double[n + 1];

	for (int i = 0; i < n + 1; i++)
	{
	

		if (variant == 20)
		{
			yArray[i] = F20(xArray[i], c20);
		}
		else if (variant == 24)
		{
			yArray[i] = F24(xArray[i], c24);
		}
		else if (variant == 19)
		{
			yArray[i] = F19(xArray[i], c19);

		}

	}




	MakeDifTable(difTable, xArray,yArray, difVector, variant);

	double x = difTable[0][1];
	double w = 1;
	double c = 1;
	double nevyazka = 0;

	for (int i = 0; i < n + 1; i++)
	{


		if (variant == 20)
		{
			yArray[i] = F20(xArray[i], 0);
		}
		else if (variant == 24)
		{
			yArray[i] = F24(xArray[i], 0);
		}
		else if (variant == 19)
		{
			yArray[i] = F19(xArray[i], 0);

		}

	}

	/*if (variant==19)
		c = F19((b-a)/2,0);
	else if (variant==20)
		c = F20((b - a) / 2, 0);
	else if (variant ==24)
		c = F24((b - a) / 2, 0);*/

	for (int i = 1; i <= n; i++)
	{
		w *= c - yArray[i - 1];
		x += difTable[0][i + 1] * w;
	}
	
	if (variant == 19)
		nevyazka = F19(x,0) - c;
	else if (variant == 20)
		nevyazka = F20(x, 0) - c;
	else if (variant == 24)
		nevyazka = F24(x, 0) - c;

	std::cout << "\nКорень: " << x << std::endl;
	std::cout << "Невязка: " << std::setprecision(8) <<abs(nevyazka) << std::endl;


}