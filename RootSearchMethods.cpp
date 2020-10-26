#include <iostream>
#include<iomanip>
#include "Matrix.h"
#include "RootSearchMethods.h"

using namespace std;

#pragma region Печать

void  RootSearchMethods::PrintTheHead(int i)
{
	int width = 14;

	switch (i)
	{
	case 0:
		cout << setw(width) << left << "Itr" << "|" << setw(width) << left << "x" << "|";
		cout << "|" << setw(width) << left << "y" << "|" << setw(width) << left << "Норма невязки" << "|";
		cout << "|" << setw(width) << left << "F1" << "|" << setw(width) << left << "F2" << endl;
		break;
	case 1:
		cout << setw(width) << left << "Itr" << "|" << setw(width) << left << "x" << "|";
		cout << "|" << setw(width) << left << "y" << "|" << setw(width) << left << "Alfa" << "|";
		cout << "|" << setw(width) << left << "Норма невязки" << "|" << setw(width) << left << "F1";
		cout << "|" << setw(width) << left << "F2" << "|" << setw(width) << left << "FF" << endl;
		break;
	case 2:
		cout << setw(width) << left << "Itr" << "|" << setw(width) << left << "x" << "|";
		cout << "|" << setw(width) << left << "y" << "|" << setw(width) << left << "Норма невязки" << "|";
		cout << "|" << setw(width) << left << "F1" << "|" << setw(width) << left << "F2" << "|" << setw(width) << left << "Норма якобиана" << endl;
		break;
	default:
		break;
	}
}

void RootSearchMethods::PrintIntermediateResults(int Itr, double q, double normaNevyazki, double ErrorEstimate, double x, double y, double F1, double F2)
{
	int width = 14;
	cout << setw(width) << left << Itr << "|" << setw(width) << left << x << "|";
	cout << "|" << setw(width) << left << y << "|" << setw(width) << left << normaNevyazki << "|";
	cout << "|" << setw(width) << left << F1 << "|" << setw(width) << left << F2 << "|" << setw(width) << left << q << endl;
}

void RootSearchMethods::PrintIntermediateResults(int Itr, double normaNevyazki, double x, double y, double F1, double F2)
{
	int width = 14;
	cout << setw(width) << left << Itr << "|" << setw(width) << left << x << "|";
	cout << "|" << setw(width) << left << y << "|" << setw(width) << left << normaNevyazki << "|";
	cout << "|" << setw(width) << left << F1 << "|" << setw(width) << left << F2 << endl;
}

void RootSearchMethods::PrintIntermediateResults1(int itr, double x, double y, double alpha, double normaNevyazki, double f1, double f2, double ff)
{
	int width = 14;
	cout << setw(width) << left << itr << "|" << setw(width) << left << x << "|";
	cout << "|" << setw(width) << left << y << "|" << setw(width) << left << alpha << "|";
	cout << "|" << setw(width) << left << normaNevyazki << "|" << setw(width) << left << f1;
	cout << "|" << setw(width) << left << f2 << "|" << setw(width) << left << ff << endl;
}

void RootSearchMethods::SimpleIterationMethod(int variant, double epsilon, double x, double y)
{
	int Itr = 0;
	double error = 1;
	double xNew = 0;
	double yNew = 0;
	double q = 0;
	double** J = new double* [2];
	double f1;
	double f2;

	cout << "x0: " << x << "; y0: " << y << endl;
	cout << endl;
	if (variant == 24)
	{
		cout << "F1(x,y): " << "cos(x) + y - 1.2" << endl;
		cout << "F2(x,y): " << "2 * x - sin(y - 0.5) - 2;" << endl;
	}
	else if (variant == 20)
	{
		cout << "F1(x,y): " << "sin(y + 2) - x -1.5" << endl;
		cout << "F2(x,y): " << "y + cos(x - 2) - 0.5" << endl;

	}
	cout << endl;

	std::cout << "Якобиан:" << std::endl;
	if (variant == 20)
	{
		std::cout << "0 cos(y+2)" << endl;
		cout << "sin(x-2) 0" << endl;

	}
	else if (variant == 24)
	{
		cout << "sin(x) 0" << endl;
		cout << "0 0.5*cos(y-0.5)" << endl;

	}
	cout << endl;

	for (int i = 0; i < 2; i++)
		J[i] = new double[2];
	std::cout << "Значение якобиана:" << std::endl;
	Jacobian(J, x, y, variant);
	Matrix::Show(J, 2);
	std::cout << std::endl;
	PrintTheHead(2);

	while (error > epsilon)
	{
		Itr++;

		if (variant == 24)
		{
			xNew = 1 + sin(y - 0.5) / 2;
			yNew = 1.2 - cos(x);
		}
		else if (variant == 20)
		{
			xNew = sin(y + 2) - 1.5;
			yNew = 0.5 - cos(x - 2);

		}
		else if (variant == 0)
		{
			xNew = -1 + sin(y + 0.5);
			yNew = -cos(x - 2);

		}

		Jacobian(J, xNew, yNew, variant);
		q = Matrix::Norm3(J, 2);
		f1 = F1(x, y, variant);
		f2 = F2(x, y, variant);
		double* vector = new double[2];
		vector[0] = f1;
		vector[1] = f2;
		error = Matrix::VectorNorm3(vector, 2);

		//error = ErrorNorm(x, xNew, y, yNew)*(1-q)/q;
		RootSearchMethods::PrintIntermediateResults(Itr, q, error, 0, xNew, yNew, f1, f2);
		x = xNew;
		y = yNew;
		delete[] vector;
	}

	Matrix::Clear(J, 2);
}
void RootSearchMethods::Jacobian(double** J, double x, double y, int variant)
{
	if (variant == 24)
	{
		J[0][0] = sin(x);
		J[1][0] = 0;
		J[0][1] = 0;
		J[1][1] = cos(y - 0.5) / 2;
	}
	else if (variant == 20)
	{
		J[0][0] = 0;
		J[1][1] = 0;
		J[0][1] = cos(y + 2);
		J[1][0] = sin(x - 2);
	}
	else if (variant == 0)
	{
		J[0][0] = 0;
		J[0][1] = cos(y + 0.5);
		J[1][0] = sin(x - 2);
		J[1][1] = 0;
	}
}

double RootSearchMethods::ErrorNorm(double x1, double x2, double y1, double y2)
{
	return sqrt(abs(x1 - x2) * abs(x1 - x2) + abs(y1 - y2) * abs(y1 - y2));

}

double RootSearchMethods::F1(double x, double y, double variant)
{
	if (variant == 24)
		return cos(x) + y - 1.2;
	else if (variant == 20)
		return sin(y + 2) - x - 1.5;
	else if (variant == 0)
		return -1 + sin(y + 0.5) - x;

}

double RootSearchMethods::F2(double x, double y, double variant)
{
	if (variant == 24)
		return 2 * x - sin(y - 0.5) - 2;
	else if (variant == 20)
		return y + cos(x - 2) - 0.5;
	else if (variant == 0)
		return cos(x - 2) + y;

}

double** RootSearchMethods::CalculateMatrixF(double x, double y, int var)
{
	double** res = new double* [2];
	for (int i = 0; i < 2; i++)
		res[i] = new double[2];

	switch (var)
	{
	case 20:
		res[0][0] = -1;
		res[0][1] = cos(y + 2);
		res[1][0] = -sin(x - 2);
		res[1][1] = 1;
		break;
	case 24:
		res[0][0] = -sin(x);
		res[0][1] = 1;
		res[1][0] = 2;
		res[1][1] = -cos(y - 0.5);
		break;
	case 0:
		res[0][0] = -1;
		res[0][1] = cos(y + 0.5);
		res[1][0] = -sin(x - 2);
		res[1][1] = 1;
		break;
	default:
		break;
	}

	return res;
}

void RootSearchMethods::CalculateReverseMatrixF(double** F)
{
	double det = F[0][0] * F[1][1] - F[0][1] * F[1][0];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			int temp = i + j + 2;
			F[i][j] /= temp % 2 == 0 ? det : -det;
		}
	}
	swap(F[0][0], F[1][1]);
}

double RootSearchMethods::CalculateTheError(double x, double y, double xNext, double yNext, double q)
{
	double xRes = (xNext - x) * (xNext - x);
	double yRes = (yNext - y) * (yNext - y);
	return sqrt(xRes + yRes) * ((1 - q) / q);
}

double RootSearchMethods::FunctionFF(double x, double y, int var)
{
	return F1(x, y, var) * F1(x, y, var) + F2(x, y, var) * F2(x, y, var);
}

void RootSearchMethods::NewtonMethod(double x, double y, double eps, int var)
{
	cout << endl << "Метод Ньютона: " << endl;
	PrintTheHead(0);
	int itr = 0;
	double xNext, yNext, q, error, normaNevyazki;
	double** J = new double* [2];
	double f1;
	double f2;

	for (int i = 0; i < 2; i++)
		J[i] = new double[2];

	do
	{
		itr++;

		double** F = CalculateMatrixF(x, y, var);
		CalculateReverseMatrixF(F);

		xNext = x - (F[0][0] * F1(x, y, var) + F[0][1] * F2(x, y, var));
		yNext = y - (F[1][0] * F1(x, y, var) + F[1][1] * F2(x, y, var));

		Jacobian(J, xNext, yNext, var);
		q = Matrix::Norm3(J, 2);
		f1 = F1(xNext, yNext, var);
		f2 = F2(xNext, yNext, var);

		double f1 = F1(x, y, var);
		double f2 = F2(x, y, var);
		double* vector = new double[2];
		vector[0] = f1;
		vector[1] = f2;
		normaNevyazki = Matrix::VectorNorm3(vector, 2);
		error = normaNevyazki * ((1 - q) / q);

		PrintIntermediateResults(itr, normaNevyazki, xNext, yNext, f1, f2);

		x = xNext;
		y = yNext;


		Matrix::Clear(F, 2);
		delete[] vector;
	} while (normaNevyazki > eps);
	Matrix::Clear(J, 2);
}

double* RootSearchMethods::CalculateVectorDF(double x, double y, int var)
{
	double* res = new double[2];

	switch (var)
	{
	case 20:
		res[0] = -2 * (sin(y + 2) - x - 1.5) - 2 * sin(x - 2) * (y + cos(x - 2) - 0.5);
		res[1] = 2 * cos(y + 2) * (sin(y + 2) - x - 1.5) + 2 * (y + cos(x - 2) - 0.5);
		break;
	case 24:
		res[0] = -2 * sin(x) * (cos(x) + y - 1.2) + 4 * (2 * x - sin(y - 0.5) - 2);
		res[1] = 2 * (cos(x) + y - 1.2) - 2 * cos(y - 0.5) * (2 * x - sin(y - 0.5) - 2);
		break;
	case 0:
		res[0] = -2 * (-1 + sin(y + 0.5) - x) - 2 * sin(x - 2) * (cos(x - 2) + y);
		res[1] = 2 * cos(y + 0.5) * (-1 + sin(y + 0.5) - x) + 2 * (cos(x - 2) + y);
		break;
	default:
		break;
	}

	return res;
}

void RootSearchMethods::GradientDescentMethods(double x, double y, int var, double eps)
{
	cout << endl << "Метод Градиентного спуска: " << endl;
	PrintTheHead(1);

	int itr = 0;
	double xNext, yNext, q, error, normaNevyazki, alpha, lambda;
	double** J = new double* [2];

	for (int i = 0; i < 2; i++)
		J[i] = new double[2];

	do
	{
		itr++;

		alpha = 1;
		lambda = 0.5;

		double* DF = CalculateVectorDF(x, y, var);
		while (FunctionFF(x - alpha * DF[0], y - alpha * DF[1], var) >= FunctionFF(x, y, var))
			alpha *= lambda;

		xNext = x - alpha * DF[0];
		yNext = y - alpha * DF[1];

		Jacobian(J, xNext, yNext, var);
		q = Matrix::Norm3(J, 2);

		double f1 = F1(x, y, var);
		double f2 = F2(x, y, var);
		double* vector = new double[2];
		vector[0] = f1;
		vector[1] = f2;
		normaNevyazki = Matrix::VectorNorm3(vector, 2);
		error = normaNevyazki * ((1 - q) / q);

		PrintIntermediateResults1(itr, xNext, yNext, alpha, normaNevyazki, F1(xNext, yNext, var), F2(xNext, yNext, var), FunctionFF(xNext, yNext, var));

		x = xNext;
		y = yNext;

		delete[] vector;
		delete[] DF;
	} while (normaNevyazki > eps);
	Matrix::Clear(J, 2);
}

#pragma endregion