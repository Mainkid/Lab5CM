#define _USE_MATH_DEFINES
#include <math.h>
#include<iostream>
#include <string>
#include "SquareApproximation.h"
#include "Matrix.h"

Functions::Functions(variants variant)
{
	this->variant = variant;
}

double Functions::F(double x)
{
	switch (variant)
	{
	case variants::var20:
		return pow(M_E, x) - pow(x, 2) + 1;
	case variants::var24:
		return 2 * pow(M_E, x) - 2 * x + 3;		
	default:
		return 0.;
	}

}

double Functions::FG1(double x) // (f, g1 = 1) - интеграл от f*g1
{
	switch (variant)
	{
	case variants::var20:
		return pow(M_E, x) - pow(x, 3) / 3 + x;
	case variants::var24:
		return 2 * pow(M_E, x) - pow(x, 2) + 3 * x;
	default:
		return 0.;
	}
}
double Functions::FG2(double x) // (f, g2 = x) - интеграл от f*g2
{
	switch (variant)
	{
	case variants::var20:
		return (x - 1) * pow(M_E, x) - pow(x, 4) / 4 + pow(x, 2) / 2;
	case variants::var24:
		return 2 * (x - 1) * pow(M_E, x) - 2 * pow(x, 3) / 3 + 3 * pow(x, 2) / 2;
	default:
		return 0.;
	}
}
double Functions::FG3(double x) // (f, g3 = x^2) - интеграл от f*g3
{
	switch (variant)
	{
	case variants::var20:
		return (pow(x, 2) - 2 * x + 2) * pow(M_E, x) - pow(x, 5) / 5. + pow(x, 3) / 3.;
	case variants::var24:
		return 2 * (pow(x, 2) - 2 * x + 2) * pow(M_E, x) - pow(x, 4) / 2. + pow(x, 3);
	default:
		return 0.;
	}
}

double Functions::FF(double x) // (f, f) - интеграл от f*f
{
	switch (variant)
	{
	case variants::var20:
		return pow(M_E, 2 * x) / 2 + 2 * (pow(M_E, x) - pow(x, 3) / 3) - 2 * (pow(x, 2) - 2 * x + 2) * pow(M_E, x) + pow(x, 5) / 5 + x;
	case variants::var24:
		return 2 * pow(M_E, 2 * x) + 6 * (2 * pow(M_E, x) - pow(x, 2)) - 8 * (x - 1) * pow(M_E, x) + 4 * pow(x, 3) / 3 + 9 * x;
	default:
		return 0.;
	}
}

Table::Table()
{
	m_table = new double* [2];
}

Table::~Table()
{
	for (int i = 0; i < 2; i++)
		if (m_table[i])
			delete[] m_table[i];
	delete[] m_table;
}

void Table::Create(variants variant, double a, double b, int n)
{
	this->n = n;
	for (int i = 0; i < 2; i++)
		m_table[i] = new double[n + 1];
	Functions f(variant);
	double h = (b - a) / n;
	for (int i = 0; i <= this->n; i++)
	{	
		double xi = a + i * h;
		m_table[0][i] = xi;
		m_table[1][i] = f.F(xi);
	}
}

void Table::Print()
{
	printf("%s", "Таблица:\n");
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			if (j == 0 && i == 0)
				printf("%s", " x[i] ");
			if (j == 0 && i == 1)
				printf("%s", " f(x[i]) ");
			printf(" %.5f ", m_table[i][j]);
		}
		printf("\n");
	}
}

double SquareApproximation::ScalarProduct(double* x, double* y, int n)
{
	double res = 0.;
	for (int i = 0; i < n; i++)
	{
		res += x[i] * y[i];
	}
	return res;
}

SquareApproximation::SquareApproximation(variants variant)
{
	this->variant = variant;
	table->Create(variant, 1, 2, 5);
}

SquareApproximation::~SquareApproximation()
{
	delete table;
}

void SquareApproximation::PrintEquation(double* C, int n)
{
	if (n < 1) return;
	printf("%s", "P2(x) = ");
	printf("%f ", C[0]);
	for (int i = 1; i < n; i++)
	{
		std::string ci = std::to_string(C[i]);
		std::string toPrint = C[i] < 0 ? "- " + ci.substr(1, ci.length()) : "+ " + ci;
		std::string xPower = "x" + (i > 1 ? "^" + std::to_string(i) : "");
		printf("%s%s ", toPrint.c_str(), xPower.c_str());
	}
}

void SquareApproximation::DiscreteMethod()
{
	printf("%s", "Дискретный Вариант\n");

	const int n = 6;

	double** G = new double* [3]; // G[0] = g1(x) = 1, G[1] = g2(x) = x, G[2] = g3(x) = x^2
	for (int i = 0; i < 3; i++)
		G[i] = new double[n];

	for (int i = 0; i < n; i++)
	{
		G[0][i] = 1;
		G[1][i] = table->m_table[0][i];
		G[2][i] = pow(table->m_table[0][i], 2);
	}

	double** matrix = new double* [3];
	for (int i = 0; i < 3; i++)
		matrix[i] = new double[3];

	for (int i = 0; i < 3; i++)
	{
		matrix[0][i] = ScalarProduct(G[0], G[i], n);
		matrix[1][i] = ScalarProduct(G[1], G[i], n);
		matrix[2][i] = ScalarProduct(G[2], G[i], n);
	}

	printf("\n");
	table->Print();
	printf("\n");
	printf("%s", "Матрица\n");
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			printf(" %.5f ", matrix[i][j]);
		}
		printf("\n");
	}

	printf("%s", "\nВектор правых частей:\n");
	double* FG = new double[3]; // Вектор правых частей: FG[0] = (f, g1), FG[1] = (f, g2), FG[2] = (f, g3)
	for (int i = 0; i < 3; i++)
	{
		FG[i] = ScalarProduct(table->m_table[1], G[i], n);
		printf(" %f ", FG[i]);
	}
	printf("\n\n");

	
	double* C = Matrix::LU(matrix, FG, 3); // нахождение с1, с2, с3 методом LU разложения
	PrintEquation(C, 3);
	printf("\n");

	double* CG = new double[n]; // вектор - c1*g1+c2*g2+c3*g3
	for (int i = 0; i < n; i++)
	{
		CG[i] = C[0] * G[0][i] + C[1] * G[1][i] + C[2] * G[2][i];
	}
	double errorNorm = sqrt(ScalarProduct(table->m_table[1], table->m_table[1], n) - ScalarProduct(CG, CG, n));
	
	printf("Норма погрешности: %f\n", errorNorm);

	for (int i = 0; i < 3; i++)
	{
		delete[] G[i];
		delete[] matrix[i];
	}
	delete[] matrix;
	delete[] G;
	delete[] FG;
	delete[] C;
	delete[] CG;
}

void SquareApproximation::IntegrationMethod()
{
	printf("%s", "Непрерывный вариант\n\n");
	const int n = 6;
	Functions functions(variant);

	double** matrix = new double* [3]{
		new double[3]{1., 3./2, 7./3},
		new double[3]{3./2, 7./3, 15./4},
		new double[3]{7./3, 15./4, 31./5}
	};

	printf("%s", "Матрица\n");
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			printf(" %.5f ", matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	printf("%s", "\nВектор правых частей:\n");
	
	double* FG = new double[3]{ // Вектор правых частей: FG[0] = (f, g1), FG[1] = (f, g2), FG[2] = (f, g3)
		functions.FG1(2.) - functions.FG1(1.),
		functions.FG2(2.) - functions.FG2(1.),
		functions.FG3(2.) - functions.FG3(1.)
	};
	for (int i = 0; i < 3; i++)
	{
		printf(" %f ", FG[i]);
	}
	printf("\n\n");

	double* C = Matrix::LU(matrix, FG, 3); // нахождение с1, с2, с3 методом LU разложения
	PrintEquation(C, 3);
	printf("\n");

	double f2 = functions.FF(2.) - functions.FF(1.); // ||f||^2 - интеграл от f*f в пределах от 1 до 2
	double g2 = pow(C[0], 2) + pow(C[1], 2) * matrix[0][2] + pow(C[2], 2) * matrix[2][2] + 2 * C[0] * C[1] * matrix[0][1] + 2 * C[0] * C[2] * matrix[0][2] + 2 * C[1] * C[2] * matrix[1][2]; // ||g||^2
	// интеграл от (с1*g1+c2*g2+c3*g3)^2 в пределах от 1 до 2
	double errorNorm = sqrt(f2 - g2);
	printf("Норма погрешности: %f\n", errorNorm);

	for (int i = 0; i < 3; i++)
	{
		delete[] matrix[i];
	}
	delete[] FG;
	delete[] C;
	delete[] matrix;
}